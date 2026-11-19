"""Concordance between two sets of effects measured on the same elements (ADR-003)."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import stats


@dataclass(frozen=True)
class Concordance:
    n: int
    spearman_rho: float
    spearman_p: float
    n_sign: int
    sign_agreement: float
    sign_p: float


def concordance(
    a: np.ndarray, b: np.ndarray, mask: np.ndarray | None = None, min_abs: float = 0.0
) -> Concordance:
    """Spearman correlation on shared elements, and sign agreement among strong effects.

    ``a`` and ``b`` are effect estimates (NaN = not measured). ``mask`` optionally
    restricts the sign test (for example to elements significant in both studies);
    ``min_abs`` drops weak effects from the sign test. The sign test is a one-sided
    binomial test of agreement above 50%.
    """
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    shared = np.isfinite(a) & np.isfinite(b)
    n = int(shared.sum())
    if n < 3:
        return Concordance(n, float("nan"), float("nan"), 0, float("nan"), float("nan"))
    rho, p = stats.spearmanr(a[shared], b[shared])
    keep = shared & (np.abs(a) >= min_abs) & (np.abs(b) >= min_abs)
    if mask is not None:
        keep &= np.asarray(mask, dtype=bool)
    n_sign = int(keep.sum())
    if n_sign == 0:
        return Concordance(n, float(rho), float(p), 0, float("nan"), float("nan"))
    agree = int((np.sign(a[keep]) == np.sign(b[keep])).sum())
    sign_p = float(stats.binomtest(agree, n_sign, 0.5, alternative="greater").pvalue)
    return Concordance(n, float(rho), float(p), n_sign, agree / n_sign, sign_p)
