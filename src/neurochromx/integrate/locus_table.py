"""Element-effect tables and the cross-study master catalogue (PRD sections 8 and 9).

``element_effects`` rows hold one estimate per element, study, contrast and modality.
Studies call peaks independently, so effects from different studies are compared on a
*master catalogue*: the union of all studies' elements, merged where they overlap.
"""

from __future__ import annotations

import math
from collections import defaultdict
from collections.abc import Iterable, Sequence
from pathlib import Path

import numpy as np

from neurochromx.intervals import IntervalIndex, Region

EFFECT_COLUMNS = [
    "element_id", "study_id", "contrast_id", "modality", "estimate", "se", "stat", "p",
    "fdr", "n_case", "n_ctrl", "qc_flag",
]  # fmt: skip


def read_deseq2_effects(
    path: str | Path,
    study_id: str,
    contrast_id: str,
    modality: str = "atac",
    n_case: int | None = None,
    n_ctrl: int | None = None,
    qc_flag: str = "pass",
) -> list[dict[str, str]]:
    """Convert a DESeq2 result table (deseq2_dar.R) to element_effects rows."""
    lines = Path(path).read_text().splitlines()
    head = lines[0].split("\t")
    rows = []
    for line in lines[1:]:
        f = dict(zip(head, line.split("\t"), strict=True))
        rows.append(
            {
                "element_id": f["element_id"],
                "study_id": study_id,
                "contrast_id": contrast_id,
                "modality": modality,
                "estimate": f["log2FC"],
                "se": f["lfcSE"],
                "stat": f["stat"],
                "p": f["pvalue"],
                "fdr": f["padj"],
                "n_case": "NA" if n_case is None else str(n_case),
                "n_ctrl": "NA" if n_ctrl is None else str(n_ctrl),
                "qc_flag": qc_flag,
            }
        )
    return rows


def write_effects(rows: Sequence[dict[str, str]], path: str | Path) -> None:
    """Write ``element_effects`` as TSV, or Parquet when the suffix is ``.parquet``."""
    path = Path(path)
    if path.suffix == ".parquet":
        import pyarrow as pa
        import pyarrow.parquet as pq

        table = pa.table({c: [r[c] for r in rows] for c in EFFECT_COLUMNS})
        pq.write_table(table, path)
        return
    path.write_text(
        "\t".join(EFFECT_COLUMNS) + "\n"
        + "".join("\t".join(r[c] for c in EFFECT_COLUMNS) + "\n" for r in rows)
    )  # fmt: skip


def read_effects(path: str | Path) -> list[dict[str, str]]:
    lines = Path(path).read_text().splitlines()
    head = lines[0].split("\t")
    return [dict(zip(head, line.split("\t"), strict=True)) for line in lines[1:]]


def build_master_catalogue(peaks_by_study: dict[str, Sequence[Region]]) -> list[Region]:
    """Union of all studies' elements, merged where they overlap or touch.

    IDs have the form ``master:<chrom>:<start>-<end>``.
    """
    merged: list[list] = []
    allr = sorted(
        (r for regions in peaks_by_study.values() for r in regions),
        key=lambda r: (r.chrom, r.start, r.end),
    )
    for r in allr:
        if merged and merged[-1][0] == r.chrom and r.start <= merged[-1][2]:
            merged[-1][2] = max(merged[-1][2], r.end)
        else:
            merged.append([r.chrom, r.start, r.end])
    return [Region(c, s, e, f"master:{c}:{s}-{e}") for c, s, e in merged]


def project_effects(
    rows: Iterable[dict[str, str]], regions: dict[str, Region], master: Sequence[Region]
) -> list[dict[str, str]]:
    """Re-key one study's effects onto master elements.

    When several study elements fall in one master element, the estimate is their
    inverse-variance weighted mean (a region-level effect); ``se`` is the combined
    standard error, ``p``/``fdr``/``stat`` are NA, and the largest ``n_case``/``n_ctrl``
    is kept. Study elements with no usable estimate are dropped.
    """
    index = IntervalIndex(master)
    groups: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        region = regions[row["element_id"]]
        hits = index.overlapping(region.chrom, region.start, region.end)
        if not hits:
            continue
        groups[hits[0].name].append(row)
    out = []
    for master_id, members in groups.items():
        usable = []
        for m in members:
            try:
                est, se = float(m["estimate"]), float(m["se"])
            except ValueError:
                continue
            if math.isfinite(est) and math.isfinite(se) and se > 0:
                usable.append((est, se, m))
        if not usable:
            continue
        if len(usable) == 1:
            row = dict(usable[0][2])
            row["element_id"] = master_id
            out.append(row)
            continue
        weights = np.array([1 / se**2 for _, se, _ in usable])
        estimate = float(
            np.sum(weights * np.array([e for e, _, _ in usable])) / weights.sum()
        )
        row = dict(usable[0][2])
        pooled_se = math.sqrt(1 / weights.sum())
        row.update(
            element_id=master_id,
            estimate=f"{estimate:.6g}",
            se=f"{pooled_se:.6g}",
            stat="NA",
            p="NA",
            fdr="NA",
        )
        out.append(row)
    return out


def effects_matrix(
    rows_by_study: dict[str, Sequence[dict[str, str]]], elements: Sequence[str]
) -> tuple[np.ndarray, np.ndarray]:
    """Estimates and standard errors as (n_studies, n_elements); NaN where missing."""
    position = {e: i for i, e in enumerate(elements)}
    est = np.full((len(rows_by_study), len(elements)), np.nan)
    se = np.full_like(est, np.nan)
    for k, rows in enumerate(rows_by_study.values()):
        for row in rows:
            j = position.get(row["element_id"])
            if j is not None:
                est[k, j] = float(row["estimate"])
                se[k, j] = float(row["se"])
    return est, se
