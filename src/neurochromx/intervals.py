"""Small interval helpers (0-based, half-open) shared by the annotation code."""

from __future__ import annotations

import bisect
import gzip
from collections import defaultdict
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Region:
    """A genomic interval. ``name`` and ``label`` carry optional BED columns 4 and 5."""

    chrom: str
    start: int
    end: int
    name: str = ""
    label: str = ""

    def __post_init__(self) -> None:
        if self.end < self.start:
            raise ValueError(f"end before start in {self}")

    @property
    def centre(self) -> int:
        return (self.start + self.end) // 2


def read_bed(path: str | Path) -> list[Region]:
    """Read a BED file (plain or .gz); columns 4 and 5 become name and label."""
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    regions: list[Region] = []
    with opener(path, "rt") as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            regions.append(
                Region(
                    f[0],
                    int(f[1]),
                    int(f[2]),
                    f[3] if len(f) > 3 else "",
                    f[4] if len(f) > 4 else "",
                )
            )
    return regions


class IntervalIndex:
    """Fast overlap queries over a fixed set of regions."""

    def __init__(self, regions: Iterable[Region]) -> None:
        by_chrom: dict[str, list[Region]] = defaultdict(list)
        for region in regions:
            by_chrom[region.chrom].append(region)
        self._regions = {
            c: sorted(v, key=lambda r: r.start) for c, v in by_chrom.items()
        }
        self._starts = {c: [r.start for r in v] for c, v in self._regions.items()}
        self._max_len = {
            c: max(r.end - r.start for r in v) for c, v in self._regions.items()
        }

    def overlapping(self, chrom: str, start: int, end: int) -> list[Region]:
        """All indexed regions sharing at least one base with ``[start, end)``."""
        if chrom not in self._regions:
            return []
        starts = self._starts[chrom]
        lo = bisect.bisect_left(starts, start - self._max_len[chrom])
        hi = bisect.bisect_left(starts, end)
        return [
            r for r in self._regions[chrom][lo:hi] if r.end > start and r.start < end
        ]

    def any_overlap(self, chrom: str, start: int, end: int) -> bool:
        return bool(self.overlapping(chrom, start, end))
