"""Promoter and enhancer-tier assignment for consensus elements (ADR-007).

Tiers are not exclusive; every result table carries all three flags:

* tier 1: overlaps an ENCODE cCRE whose class is in ``ccre_classes`` (pELS, dELS)
* tier 2: overlaps H3K27ac and H3K4me1 peaks and no H3K4me3 peak (same-tissue ChIP)
* tier 3: a distal element, farther than ``tier3_min_tss_distance`` from any TSS

An element that overlaps a promoter window is ``promoter``; otherwise it is ``enhancer``
if any tier flag is set, else ``other``.
"""

from __future__ import annotations

import bisect
from collections import defaultdict
from collections.abc import Collection, Iterable, Sequence
from dataclasses import dataclass

from neurochromx.annotate.promoters import Gene
from neurochromx.intervals import IntervalIndex, Region


@dataclass(frozen=True)
class ElementAnnotation:
    element_id: str
    chrom: str
    start: int
    end: int
    cls: str  # promoter | enhancer | other
    tier1: bool
    tier2: bool
    tier3: bool
    tss_distance: int


class TssIndex:
    """Distance from an interval to the nearest TSS."""

    def __init__(self, genes: Iterable[Gene]) -> None:
        by_chrom: dict[str, list[int]] = defaultdict(list)
        for gene in genes:
            by_chrom[gene.chrom].append(gene.tss)
        self._tss = {c: sorted(v) for c, v in by_chrom.items()}

    def distance(self, chrom: str, start: int, end: int) -> int | None:
        """0 if a TSS lies inside ``[start, end)``; None if the chromosome has none."""
        positions = self._tss.get(chrom)
        if not positions:
            return None
        i = bisect.bisect_left(positions, start)
        best = None
        for j in (i - 1, i):
            if 0 <= j < len(positions):
                tss = positions[j]
                gap = max(0, start - tss, tss - (end - 1))
                best = gap if best is None else min(best, gap)
        return best


def assign_elements(
    elements: Sequence[Region],
    promoters: Iterable[Region],
    genes: Iterable[Gene],
    ccres: Iterable[Region] = (),
    h3k27ac: Iterable[Region] | None = None,
    h3k4me1: Iterable[Region] | None = None,
    h3k4me3: Iterable[Region] | None = None,
    ccre_classes: Collection[str] = ("pELS", "dELS"),
    tier3_min_tss_distance: int = 2000,
) -> list[ElementAnnotation]:
    """Annotate each element with its class, tier flags and distance to the nearest TSS.

    ``Region.name`` is the element ID. Tier 2 needs all three histone peak sets and is
    False for every element when any of them is None (no same-tissue ChIP available).
    """
    promoter_index = IntervalIndex(promoters)
    ccre_index = IntervalIndex(r for r in ccres if r.label in set(ccre_classes))
    marks = (h3k27ac, h3k4me1, h3k4me3)
    mark_index = (
        [IntervalIndex(m) for m in marks] if all(m is not None for m in marks) else None
    )
    tss_index = TssIndex(genes)

    out = []
    for element in elements:
        c, s, e = element.chrom, element.start, element.end
        distance = tss_index.distance(c, s, e)
        tier1 = ccre_index.any_overlap(c, s, e)
        tier2 = False
        if mark_index is not None:
            ac, me1, me3 = (index.any_overlap(c, s, e) for index in mark_index)
            tier2 = ac and me1 and not me3
        tier3 = distance is not None and distance > tier3_min_tss_distance
        if promoter_index.any_overlap(c, s, e):
            cls = "promoter"
        elif tier1 or tier2 or tier3:
            cls = "enhancer"
        else:
            cls = "other"
        out.append(
            ElementAnnotation(
                element.name,
                c,
                s,
                e,
                cls,
                tier1,
                tier2,
                tier3,
                -1 if distance is None else distance,
            )
        )
    return out
