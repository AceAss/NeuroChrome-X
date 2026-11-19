"""Overlap of two element sets against a covariate-matched permutation null (ADR-003).

Question: are elements changed in one setting (for example aging) over-represented
among elements changed in another (for example neuronal activation)? A naive test is
confounded when both depend on a covariate such as accessibility: highly accessible
elements are easier to call in either study. The null here redraws the second set
*within covariate cells*, keeping its covariate profile, so only dependence beyond the
covariates counts.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class OverlapResult:
    n_universe: int
    n_a: int
    n_b: int
    observed: int
    expected: float
    fold: float
    pvalue: float
    n_perm: int


def covariate_cells(covariates: np.ndarray, n_bins: int = 5) -> np.ndarray:
    """Integer cell ID per row: quantile bins of each covariate, combined."""
    covariates = np.atleast_2d(np.asarray(covariates, dtype=float))
    if covariates.shape[0] == 1 and covariates.shape[1] > 1:
        covariates = covariates.T
    cell = np.zeros(len(covariates), dtype=np.int64)
    for j in range(covariates.shape[1]):
        edges = np.quantile(covariates[:, j], np.linspace(0, 1, n_bins + 1)[1:-1])
        cell = cell * n_bins + np.searchsorted(edges, covariates[:, j], side="right")
    return cell


def matched_overlap_test(
    in_a: np.ndarray,
    in_b: np.ndarray,
    covariates: np.ndarray | None,
    rng: np.random.Generator,
    n_perm: int = 10_000,
    n_bins: int = 5,
) -> OverlapResult:
    """One-sided permutation test that |A and B| exceeds the matched-null expectation.

    ``in_a`` and ``in_b`` are boolean membership vectors over the same universe.
    ``covariates`` is (n, k) or None for an unmatched test. The p-value is
    ``(1 + #{null >= observed}) / (1 + n_perm)``.
    """
    in_a, in_b = np.asarray(in_a, dtype=bool), np.asarray(in_b, dtype=bool)
    n = len(in_a)
    cells = (
        np.zeros(n, dtype=np.int64)
        if covariates is None
        else covariate_cells(covariates, n_bins)
    )
    members = {c: np.flatnonzero(cells == c) for c in np.unique(cells)}
    b_in_cell = {c: int(in_b[idx].sum()) for c, idx in members.items()}
    observed = int((in_a & in_b).sum())
    null = np.zeros(n_perm, dtype=np.int64)
    for k in range(n_perm):
        count = 0
        for c, idx in members.items():
            take = b_in_cell[c]
            if take:
                drawn = rng.choice(idx, size=take, replace=False)
                count += int(in_a[drawn].sum())
        null[k] = count
    expected = float(null.mean())
    p = (1 + int((null >= observed).sum())) / (1 + n_perm)
    fold = observed / expected if expected > 0 else float("inf")
    return OverlapResult(
        n, int(in_a.sum()), int(in_b.sum()), observed, expected, fold, p, n_perm
    )
