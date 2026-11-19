"""Motif scanning with position weight matrices (both strands).

This is a small, dependency-light scanner used for the toy-dataset gates and for
targeted tests. Production enrichment on real data uses HOMER and FIMO/AME with JASPAR
matrices (PRD section 7); this code agrees with them on the logic, not on thresholds.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np

IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT", "S": "CG",
    "W": "AT", "K": "GT", "M": "AC", "N": "ACGT",
}  # fmt: skip
_CODE = np.full(256, 4, dtype=np.int8)
for _i, _b in enumerate(b"ACGT"):
    _CODE[_b] = _i
    _CODE[ord(chr(_b).lower())] = _i
UNIFORM = (0.25, 0.25, 0.25, 0.25)


@dataclass(frozen=True)
class Motif:
    """A motif as log2-odds scores, shape (length, 4), columns A, C, G, T."""

    name: str
    log_odds: np.ndarray

    @property
    def length(self) -> int:
        return int(self.log_odds.shape[0])

    @property
    def max_score(self) -> float:
        return float(self.log_odds.max(axis=1).sum())

    @property
    def min_score(self) -> float:
        return float(self.log_odds.min(axis=1).sum())


def motif_from_probabilities(
    name: str, probs: np.ndarray, background: Sequence[float] = UNIFORM
) -> Motif:
    probs = np.asarray(probs, dtype=float)
    probs = probs / probs.sum(axis=1, keepdims=True)
    return Motif(name, np.log2(probs / np.asarray(background)[None, :]))


def motif_from_consensus(
    name: str,
    consensus: str,
    strength: float = 0.9,
    background: Sequence[float] = UNIFORM,
) -> Motif:
    """PWM from an IUPAC consensus: ``strength`` spread over the allowed bases."""
    rows = []
    for symbol in consensus.upper():
        allowed = IUPAC[symbol]
        row = np.full(
            4, (1 - strength) / max(1, 4 - len(allowed)) if len(allowed) < 4 else 0.25
        )
        for base in allowed:
            row[("ACGT").index(base)] = (
                strength / len(allowed) if len(allowed) < 4 else 0.25
            )
        rows.append(row)
    return motif_from_probabilities(name, np.array(rows), background)


def motif_from_counts(
    name: str,
    counts: np.ndarray,
    pseudocount: float = 0.8,
    background: Sequence[float] = UNIFORM,
) -> Motif:
    """PWM from a JASPAR-style count matrix of shape (4, length), rows A, C, G, T."""
    probs = (
        np.asarray(counts, dtype=float).T + pseudocount * np.asarray(background)
    ) / (np.asarray(counts, dtype=float).T.sum(axis=1, keepdims=True) + pseudocount)
    return motif_from_probabilities(name, probs, background)


def _window_scores(codes: np.ndarray, log_odds: np.ndarray) -> np.ndarray:
    length = log_odds.shape[0]
    extended = np.hstack(
        [log_odds, log_odds.min(axis=1, keepdims=True)]
    )  # N = worst base
    n_windows = len(codes) - length + 1
    scores = np.zeros(n_windows)
    for j in range(length):
        scores += extended[j, codes[j : j + n_windows]]
    return scores


def best_score(seq: str, motif: Motif) -> float:
    """Highest log-odds score over all windows on both strands (-inf if too short)."""
    codes = _CODE[np.frombuffer(seq.encode(), dtype=np.uint8)].astype(np.int64)
    if len(codes) < motif.length:
        return float("-inf")
    forward = _window_scores(codes, motif.log_odds).max()
    reverse = _window_scores(codes, motif.log_odds[::-1, ::-1]).max()
    return float(max(forward, reverse))


def relative_score(seq: str, motif: Motif) -> float:
    """Best score rescaled to 0 (worst possible match) .. 1 (perfect match)."""
    score = best_score(seq, motif)
    return (score - motif.min_score) / (motif.max_score - motif.min_score)


def has_hit(seq: str, motif: Motif, min_rel_score: float = 0.9) -> bool:
    return relative_score(seq, motif) >= min_rel_score
