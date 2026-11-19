"""Background selection for motif enrichment (PRD section 7, ADR-008).

The foreground is the set of differential elements; the background is drawn from
non-differential accessible elements, matched on covariates such as GC content and
accessibility, so that generic open-chromatin motifs do not look enriched.
"""

from __future__ import annotations

import numpy as np


def gc_content(seq: str) -> float:
    seq = seq.upper()
    valid = sum(seq.count(b) for b in "ACGT")
    return (seq.count("G") + seq.count("C")) / valid if valid else float("nan")


def match_background(
    foreground: np.ndarray,
    pool: np.ndarray,
    rng: np.random.Generator,
    ratio: int = 1,
    replace: bool = False,
) -> np.ndarray:
    """Greedy nearest-neighbour matching on standardised covariates.

    ``foreground`` is (n, k) and ``pool`` is (m, k). Returns the indices of
    ``n * ratio`` pool rows. Foreground rows are visited in random order; each takes
    the closest unused pool rows. With ``replace=False`` the pool must hold at least
    ``n * ratio`` rows.
    """
    foreground = np.atleast_2d(np.asarray(foreground, dtype=float))
    pool = np.atleast_2d(np.asarray(pool, dtype=float))
    n, m = len(foreground), len(pool)
    if not replace and m < n * ratio:
        raise ValueError(
            f"pool has {m} rows but {n * ratio} are needed; use replace=True"
        )
    mean = pool.mean(axis=0)
    std = pool.std(axis=0)
    std[std == 0] = 1.0
    fg = (foreground - mean) / std
    bg = (pool - mean) / std
    used = np.zeros(m, dtype=bool)
    chosen: list[int] = []
    for i in rng.permutation(n):
        distance = ((bg - fg[i]) ** 2).sum(axis=1)
        if not replace:
            distance = np.where(used, np.inf, distance)
        for j in np.argsort(distance, kind="stable")[:ratio]:
            chosen.append(int(j))
            used[j] = True
    return np.array(chosen, dtype=int)
