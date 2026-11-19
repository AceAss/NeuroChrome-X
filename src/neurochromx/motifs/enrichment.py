"""Motif enrichment of foreground vs matched background sequences.

A sequence counts as a hit when its best score on either strand reaches
``min_rel_score`` of the motif's score range. Presence/absence is compared with a
one-sided Fisher exact test and corrected across the tested motifs (BH).
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

from scipy.stats import fisher_exact

from neurochromx.motifs.scan import Motif, has_hit


@dataclass(frozen=True)
class EnrichmentResult:
    motif: str
    n_fg: int
    n_bg: int
    hits_fg: int
    hits_bg: int
    frac_fg: float
    frac_bg: float
    odds_ratio: float
    pvalue: float
    padj: float


def benjamini_hochberg(pvalues: Sequence[float]) -> list[float]:
    m = len(pvalues)
    order = sorted(range(m), key=lambda i: pvalues[i])
    adjusted = [0.0] * m
    running = 1.0
    for rank in range(m, 0, -1):
        i = order[rank - 1]
        running = min(running, pvalues[i] * m / rank)
        adjusted[i] = running
    return adjusted


def motif_enrichment(
    fg_seqs: Sequence[str],
    bg_seqs: Sequence[str],
    motifs: Sequence[Motif],
    min_rel_score: float = 0.9,
) -> list[EnrichmentResult]:
    """Test every motif; results are sorted by p-value."""
    raw = []
    for motif in motifs:
        a = sum(has_hit(s, motif, min_rel_score) for s in fg_seqs)
        c = sum(has_hit(s, motif, min_rel_score) for s in bg_seqs)
        b, d = len(fg_seqs) - a, len(bg_seqs) - c
        _, p = fisher_exact([[a, b], [c, d]], alternative="greater")
        odds = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))  # Haldane correction
        raw.append((motif.name, a, c, odds, float(p)))
    padj = benjamini_hochberg([r[4] for r in raw])
    results = [
        EnrichmentResult(
            name,
            len(fg_seqs),
            len(bg_seqs),
            a,
            c,
            a / len(fg_seqs),
            c / len(bg_seqs),
            odds,
            p,
            q,
        )
        for (name, a, c, odds, p), q in zip(raw, padj, strict=True)
    ]
    return sorted(results, key=lambda r: r.pvalue)
