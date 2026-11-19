"""Enhancer-to-gene links: nearest expressed TSS within a window (ADR-007, v0.1)."""

from __future__ import annotations

import bisect
from collections import defaultdict
from collections.abc import Collection, Sequence
from dataclasses import dataclass

from neurochromx.annotate.promoters import Gene
from neurochromx.intervals import Region


@dataclass(frozen=True)
class Link:
    element_id: str
    gene_id: str
    gene_name: str
    distance: int
    method: str = "nearest_expressed_tss"


def link_nearest_expressed_tss(
    elements: Sequence[Region],
    genes: Sequence[Gene],
    expressed: Collection[str] | None = None,
    window: int = 100_000,
) -> dict[str, Link | None]:
    """Link each element to the nearest expressed gene's TSS within ``window`` bases.

    Distance is measured from the element centre to the TSS. ``expressed`` is a set of
    gene IDs (None means every gene counts). Ties go to the smaller gene ID. Elements
    with no expressed gene in range map to None. ``Region.name`` is the element ID.
    """
    allowed = None if expressed is None else set(expressed)
    by_chrom: dict[str, list[tuple[int, str, Gene]]] = defaultdict(list)
    for gene in genes:
        if allowed is None or gene.gene_id in allowed:
            by_chrom[gene.chrom].append((gene.tss, gene.gene_id, gene))
    for entries in by_chrom.values():
        entries.sort(key=lambda t: (t[0], t[1]))
    keys = {c: [t[0] for t in v] for c, v in by_chrom.items()}

    links: dict[str, Link | None] = {}
    for element in elements:
        centre = element.centre
        entries = by_chrom.get(element.chrom, [])
        best: tuple[int, str, Gene] | None = None
        if entries:
            i = bisect.bisect_left(keys[element.chrom], centre)
            for j in range(max(0, i - 2), min(len(entries), i + 2)):
                distance = abs(entries[j][0] - centre)
                candidate = (distance, entries[j][1], entries[j][2])
                if best is None or candidate[:2] < best[:2]:
                    best = candidate
        if best is not None and best[0] <= window:
            links[element.name] = Link(
                element.name, best[2].gene_id, best[2].gene_name, best[0]
            )
        else:
            links[element.name] = None
    return links
