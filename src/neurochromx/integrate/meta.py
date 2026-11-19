"""Random-effects meta-analysis of per-element effects across studies (ADR-002).

DerSimonian-Laird estimator: pooled = sum(w * y) / sum(w) with w = 1 / (se^2 + tau^2).
Elements measured in fewer than two studies get NaN heterogeneity and a fixed-effect
estimate. Combine only effects from the same kind of contrast.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import stats


@dataclass(frozen=True)
class MetaResult:
    estimate: np.ndarray
    se: np.ndarray
    z: np.ndarray
    p: np.ndarray
    tau2: np.ndarray
    i2: np.ndarray
    k: np.ndarray


def random_effects(estimates: np.ndarray, ses: np.ndarray) -> MetaResult:
    """Pool (n_studies, n_elements) estimates; NaN entries are skipped per element."""
    y = np.asarray(estimates, dtype=float)
    se = np.asarray(ses, dtype=float)
    valid = np.isfinite(y) & np.isfinite(se) & (se > 0)
    k = valid.sum(axis=0)
    yy = np.where(valid, y, 0.0)
    w = np.where(valid, 1.0 / np.where(valid, se, 1.0) ** 2, 0.0)
    sw = w.sum(axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        fixed = (w * yy).sum(axis=0) / sw
        q = (w * (yy - fixed) ** 2).sum(axis=0)
        c = sw - (w**2).sum(axis=0) / sw
        tau2 = np.where(k > 1, np.maximum(0.0, (q - (k - 1)) / c), np.nan)
        w_re = np.where(
            valid, 1.0 / (np.where(valid, se, 1.0) ** 2 + np.nan_to_num(tau2)), 0.0
        )
        sw_re = w_re.sum(axis=0)
        est = (w_re * yy).sum(axis=0) / sw_re
        pooled_se = np.sqrt(1.0 / sw_re)
        z = est / pooled_se
        i2 = np.where(k > 1, np.maximum(0.0, (q - (k - 1)) / q), np.nan)
    est = np.where(k > 0, est, np.nan)
    pooled_se = np.where(k > 0, pooled_se, np.nan)
    p = 2 * stats.norm.sf(np.abs(z))
    return MetaResult(
        est,
        pooled_se,
        np.where(k > 0, z, np.nan),
        np.where(k > 0, p, np.nan),
        tau2,
        i2,
        k,
    )
