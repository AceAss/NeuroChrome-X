"""Tests: effect tables, catalogue, concordance, matched overlap, meta-analysis."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from neurochromx.integrate.concordance import concordance
from neurochromx.integrate.enrichment import covariate_cells, matched_overlap_test
from neurochromx.integrate.locus_table import (
    EFFECT_COLUMNS,
    build_master_catalogue,
    effects_matrix,
    project_effects,
    read_deseq2_effects,
    read_effects,
    write_effects,
)
from neurochromx.integrate.meta import random_effects
from neurochromx.intervals import Region

DESEQ2 = (
    "element_id\tchrom\tstart\tend\tbaseMean\tlog2FC_mle\tlog2FC\tlfcSE\tstat\tpvalue\tpadj\tcall\n"
    "e1\tchr1\t100\t200\t50\t1.0\t1.0\t0.2\t5.0\t1e-6\t1e-5\tup\n"
    "e2\tchr1\t150\t260\t60\t0.5\t0.5\t0.1\t5.0\t1e-6\t1e-5\tup\n"
    "e3\tchr1\t900\t1000\t10\tNA\tNA\tNA\tNA\tNA\tNA\tns\n"
    "e4\tchr2\t10\t60\t80\t-1.0\t-1.0\t0.25\t-4.0\t1e-4\t1e-3\tdown\n"
)


def test_effect_rows_roundtrip_and_projection(tmp_path: Path) -> None:
    dar = tmp_path / "dar.tsv"
    dar.write_text(DESEQ2)
    rows = read_deseq2_effects(dar, "su2017", "stim_vs_ctrl", n_case=3, n_ctrl=3)
    assert [r["element_id"] for r in rows] == ["e1", "e2", "e3", "e4"]
    assert (
        rows[0]["estimate"] == "1.0"
        and rows[0]["se"] == "0.2"
        and rows[0]["n_case"] == "3"
    )
    out = tmp_path / "effects.tsv"
    write_effects(rows, out)
    assert read_effects(out) == rows
    assert out.read_text().splitlines()[0].split("\t") == EFFECT_COLUMNS

    regions = {
        "e1": Region("chr1", 100, 200), "e2": Region("chr1", 150, 260),
        "e3": Region("chr1", 900, 1000), "e4": Region("chr2", 10, 60),
    }  # fmt: skip
    other = [Region("chr1", 120, 300), Region("chr1", 5000, 5100)]
    master = build_master_catalogue({"su": list(regions.values()), "other": other})
    assert [m.name for m in master] == [
        "master:chr1:100-300",
        "master:chr1:900-1000",
        "master:chr1:5000-5100",
        "master:chr2:10-60",
    ]
    projected = {r["element_id"]: r for r in project_effects(rows, regions, master)}
    assert set(projected) == {
        "master:chr1:100-300",
        "master:chr2:10-60",
    }  # e3 has no estimate
    merged = projected["master:chr1:100-300"]  # e1 (1.0, se 0.2) and e2 (0.5, se 0.1)
    expected = (1.0 / 0.04 + 0.5 / 0.01) / (1 / 0.04 + 1 / 0.01)
    assert float(merged["estimate"]) == pytest.approx(expected, rel=1e-4)
    assert float(merged["se"]) == pytest.approx(
        (1 / (1 / 0.04 + 1 / 0.01)) ** 0.5, rel=1e-4
    )
    assert merged["p"] == "NA"
    est, se = effects_matrix(
        {"a": list(projected.values()), "b": []}, ["master:chr2:10-60", "x"]
    )
    assert (
        est.shape == (2, 2)
        and est[0, 0] == -1.0
        and np.isnan(est[1]).all()
        and np.isnan(se[0, 1])
    )


def test_parquet_roundtrip(tmp_path: Path) -> None:
    pytest.importorskip("pyarrow")
    import pyarrow.parquet as pq

    dar = tmp_path / "dar.tsv"
    dar.write_text(DESEQ2)
    rows = read_deseq2_effects(dar, "s", "c")
    write_effects(rows, tmp_path / "e.parquet")
    table = pq.read_table(tmp_path / "e.parquet")
    assert table.column_names == EFFECT_COLUMNS and table.num_rows == 4


def test_concordance_detects_agreement_and_handles_edge_cases() -> None:
    rng = np.random.default_rng(1)
    a = rng.normal(0, 1, 500)
    b = 0.7 * a + rng.normal(0, 0.7, 500)
    result = concordance(a, b)
    assert result.spearman_rho > 0.5 and result.spearman_p < 1e-10
    assert result.sign_agreement > 0.7 and result.sign_p < 1e-6
    independent = concordance(a, rng.normal(0, 1, 500))
    assert abs(independent.spearman_rho) < 0.15 and independent.sign_p > 0.01
    a2 = a.copy()
    a2[:50] = np.nan
    assert concordance(a2, b).n == 450
    strong = concordance(a, b, min_abs=1.0)
    assert (
        strong.n_sign < result.n_sign and strong.sign_agreement >= result.sign_agreement
    )
    assert concordance(a, b, mask=a > 5).n_sign == 0
    assert np.isnan(
        concordance(np.array([1.0, 2.0]), np.array([1.0, 2.0])).spearman_rho
    )


def _confounded_sets(
    rng: np.random.Generator, n: int = 4000
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    accessibility = rng.normal(0, 1, n)
    prob = 1 / (1 + np.exp(-(1.5 * accessibility - 2.0)))
    return accessibility, rng.random(n) < prob, rng.random(n) < prob


def test_matched_test_removes_covariate_confounding() -> None:
    rng = np.random.default_rng(3)
    accessibility, in_a, in_b = _confounded_sets(rng)  # independent given accessibility
    naive = matched_overlap_test(in_a, in_b, None, np.random.default_rng(4), n_perm=400)
    matched = matched_overlap_test(
        in_a,
        in_b,
        accessibility[:, None],
        np.random.default_rng(4),
        n_perm=400,
        n_bins=10,
    )
    assert (
        naive.pvalue < 0.01 and naive.fold > 1.3
    )  # looks enriched only because of the covariate
    assert matched.pvalue > 0.01 and 0.8 < matched.fold < 1.25


def test_matched_test_finds_real_dependence() -> None:
    rng = np.random.default_rng(5)
    accessibility, in_a, in_b = _confounded_sets(rng)
    boost = in_a & (rng.random(len(in_a)) < 0.4)  # B also takes 40% of A's members
    result = matched_overlap_test(
        in_a,
        in_b | boost,
        accessibility[:, None],
        np.random.default_rng(6),
        n_perm=400,
        n_bins=10,
    )
    assert result.pvalue < 0.01 and result.fold > 1.3
    assert result.observed > result.expected and result.n_perm == 400


def test_covariate_cells_are_balanced_quantile_bins() -> None:
    x = np.linspace(0, 1, 1000)
    cells = covariate_cells(x[:, None], n_bins=5)
    assert set(np.unique(cells)) == {0, 1, 2, 3, 4}
    assert all(abs(int((cells == c).sum()) - 200) <= 2 for c in range(5))
    two = covariate_cells(np.column_stack([x, x[::-1]]), n_bins=4)
    assert (
        len(np.unique(two)) == 4 * 4 - 12 or len(np.unique(two)) > 4
    )  # anti-correlated: fewer than 16 cells


def test_random_effects_recovers_effect_and_heterogeneity() -> None:
    rng = np.random.default_rng(7)
    n_elements, k = 3000, 6
    se = np.full((k, n_elements), 0.1)
    homogeneous = 0.8 + rng.normal(0, 1, (k, n_elements)) * se
    r0 = random_effects(homogeneous, se)
    assert abs(np.nanmean(r0.estimate) - 0.8) < 0.01
    assert np.nanmean(r0.tau2) < 0.003 and np.nanmean(r0.i2) < 0.3
    tau = 0.5
    heterogeneous = (
        0.8
        + rng.normal(0, tau, (k, n_elements))
        + rng.normal(0, 1, (k, n_elements)) * se
    )
    r1 = random_effects(heterogeneous, se)
    assert abs(np.nanmean(r1.tau2) - tau**2) < 0.03 and np.nanmean(r1.i2) > 0.8
    assert np.nanmean(r1.se) > 2 * np.nanmean(
        r0.se
    )  # heterogeneity widens the pooled interval
    assert (r1.k == k).all()


def test_random_effects_missing_data() -> None:
    est = np.array([[1.0, np.nan, np.nan], [1.0, 2.0, np.nan], [np.nan, 4.0, np.nan]])
    se = np.full_like(est, 0.5)
    r = random_effects(est, se)
    assert list(r.k) == [2, 2, 0]
    assert r.estimate[0] == pytest.approx(1.0) and r.estimate[1] == pytest.approx(3.0)
    assert np.isnan(r.estimate[2]) and np.isnan(r.p[2])
    single = random_effects(np.array([[2.0]]), np.array([[0.5]]))
    assert single.estimate[0] == pytest.approx(2.0) and np.isnan(single.tau2[0])
