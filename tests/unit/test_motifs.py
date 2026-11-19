"""Tests for motif scanning, background matching and enrichment (toy truth)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from neurochromx.annotate.promoters import read_genes
from neurochromx.motifs.backgrounds import gc_content, match_background
from neurochromx.motifs.enrichment import benjamini_hochberg, motif_enrichment
from neurochromx.motifs.scan import (
    best_score,
    has_hit,
    motif_from_consensus,
    motif_from_counts,
    relative_score,
)

AP1 = motif_from_consensus("AP1", "TGAGTCA")
CREB = motif_from_consensus("CREB", "TGACGTCA")
MEF2 = motif_from_consensus("MEF2", "CTAAAAATAG")


def revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def test_scan_exact_reverse_strand_mismatch_and_edge_cases() -> None:
    assert relative_score("AAATGAGTCAAA", AP1) == pytest.approx(1.0)
    assert relative_score("AAA" + revcomp("TGAGTCA") + "AAA", AP1) == pytest.approx(1.0)
    assert 0.8 < relative_score("AAATGAGTCTAA", AP1) < 0.9  # one mismatch
    assert not has_hit("AAATGAGTCTAA", AP1, 0.9) and has_hit("AAATGAGTCTAA", AP1, 0.85)
    assert relative_score("NNNNNNNNNNNN", AP1) == pytest.approx(0.0)
    assert best_score("ACG", AP1) == float("-inf")
    assert has_hit("aaatgagtcaaa", AP1)  # lower case is read as bases


def test_motif_from_counts_matches_consensus_logic() -> None:
    counts = np.zeros((4, 4))
    for pos, base in enumerate("ACGT"):
        counts[("ACGT").index(base), pos] = 20
    motif = motif_from_counts("tiny", counts)
    assert relative_score("TTACGTTT", motif) == pytest.approx(1.0)
    assert relative_score("TTTTTTTT", motif) < 0.5


def test_gc_content() -> None:
    assert gc_content("GGCC") == 1.0 and gc_content("ATAT") == 0.0
    assert gc_content("ACGTNN") == 0.5


def test_background_matching_balances_covariates() -> None:
    rng = np.random.default_rng(1)
    fg = np.column_stack([rng.normal(0.50, 0.03, 60), rng.normal(5.0, 1.0, 60)])
    pool = np.column_stack([rng.normal(0.40, 0.06, 1500), rng.normal(4.0, 1.5, 1500)])
    picked = match_background(fg, pool, np.random.default_rng(2), ratio=1)
    assert len(picked) == len(set(picked.tolist())) == 60
    assert abs(pool[:, 0].mean() - fg[:, 0].mean()) > 0.08
    assert abs(pool[picked, 0].mean() - fg[:, 0].mean()) < 0.01
    assert abs(pool[picked, 1].mean() - fg[:, 1].mean()) < 0.15
    with pytest.raises(ValueError):
        match_background(fg, pool[:30], np.random.default_rng(3))
    assert (
        len(match_background(fg, pool[:30], np.random.default_rng(3), replace=True))
        == 60
    )


def test_benjamini_hochberg_known_values() -> None:
    assert benjamini_hochberg([0.01, 0.04, 0.03, 0.005]) == pytest.approx(
        [0.02, 0.04, 0.04, 0.02]
    )


def _toy_sequences(
    out: Path,
) -> tuple[list[str], list[str], np.ndarray, np.ndarray, list[str], list[str]]:
    genome = "".join(
        line for line in (out / "toy.fa").read_text().split(">chrM")[0].splitlines()[1:]
    )
    lines = (out / "truth_regions.tsv").read_text().splitlines()
    head = lines[0].split("\t")
    rows = [dict(zip(head, line.split("\t"), strict=True)) for line in lines[1:]]
    gained = [r for r in rows if r["status"] == "gained"]
    unchanged = [r for r in rows if r["status"] == "unchanged"]

    def seq(r: dict[str, str]) -> str:
        return genome[int(r["start"]) : int(r["end"])]

    def cov(rs: list[dict[str, str]]) -> np.ndarray:
        return np.array(
            [[gc_content(seq(r)), float(int(r["end"]) - int(r["start"]))] for r in rs]
        )

    return (
        [seq(r) for r in gained],
        [seq(r) for r in unchanged],
        cov(gained),
        cov(unchanged),
        gained,
        unchanged,
    )


def test_ap1_is_the_top_enriched_motif_in_gained_peaks(toy_small: Path) -> None:
    fg, pool, fg_cov, pool_cov, _, _ = _toy_sequences(toy_small)
    picked = match_background(fg_cov, pool_cov, np.random.default_rng(5), ratio=2)
    bg = [pool[i] for i in picked]
    results = motif_enrichment(fg, bg, [AP1, CREB, MEF2])
    assert results[0].motif == "AP1"
    assert results[0].padj < 1e-8 and results[0].odds_ratio > 10
    assert results[0].frac_fg >= 0.75 and results[0].frac_bg <= 0.2
    assert all(r.hits_fg > 0 for r in results)


def test_no_enrichment_for_a_random_split_of_unchanged_peaks(toy_small: Path) -> None:
    _, pool, _, _, _, _ = _toy_sequences(toy_small)
    rng = np.random.default_rng(11)
    order = rng.permutation(len(pool))
    fg = [pool[i] for i in order[:48]]
    bg = [pool[i] for i in order[48:]]
    results = motif_enrichment(fg, bg, [AP1, CREB, MEF2])
    assert min(r.padj for r in results) > 0.05


def test_gene_names_still_available_for_ieg_checks(toy_small: Path) -> None:
    names = {g.gene_name for g in read_genes(toy_small / "toy.gtf")}
    assert "Fos" in names
