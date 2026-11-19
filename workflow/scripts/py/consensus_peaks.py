#!/usr/bin/env python3
"""Consensus peaks across replicates (PRD section 7, ADR-008).

A merged interval is kept when peaks from at least ``min_replicates`` samples of the
same condition overlap it; intervals are then merged across conditions. Output is a
BED file whose name column is a stable element ID ``<build>:<chrom>:<start>-<end>``.

    python consensus_peaks.py --samples samples.tsv \
        --peaks-pattern 'results/peaks/{sample}_peaks.narrowPeak' \
        --build mm10 --min-replicates 2 --out consensus.bed

``samples.tsv`` needs the columns ``sample_id`` and ``condition`` (extra columns are
ignored); list only the samples of one assay.
"""

from __future__ import annotations

import argparse
import bisect
import csv
from collections import defaultdict
from pathlib import Path

Interval = tuple[str, int, int]


def read_narrowpeak(path: Path) -> list[Interval]:
    out: list[Interval] = []
    for line in path.read_text().splitlines():
        if line.strip() and not line.startswith(("#", "track")):
            fields = line.split("\t")
            out.append((fields[0], int(fields[1]), int(fields[2])))
    return sorted(out)


def merge(intervals: list[Interval]) -> list[Interval]:
    merged: list[list] = []
    for chrom, start, end in sorted(intervals):
        if merged and merged[-1][0] == chrom and start <= merged[-1][2]:
            merged[-1][2] = max(merged[-1][2], end)
        else:
            merged.append([chrom, start, end])
    return [(c, s, e) for c, s, e in merged]


class PeakIndex:
    """Overlap queries on one sample's peaks."""

    def __init__(self, peaks: list[Interval]) -> None:
        self.by_chrom: dict[str, list[tuple[int, int]]] = defaultdict(list)
        for chrom, start, end in peaks:
            self.by_chrom[chrom].append((start, end))
        self.starts = {c: [s for s, _ in v] for c, v in self.by_chrom.items()}
        self.max_len = {c: max(e - s for s, e in v) for c, v in self.by_chrom.items()}

    def overlaps(self, chrom: str, start: int, end: int) -> bool:
        if chrom not in self.by_chrom:
            return False
        starts = self.starts[chrom]
        lo = bisect.bisect_left(starts, start - self.max_len[chrom])
        hi = bisect.bisect_left(starts, end)
        return any(e > start for _, e in self.by_chrom[chrom][lo:hi])


def consensus(
    peaks: dict[str, list[Interval]],
    condition_of: dict[str, str],
    min_replicates: int = 2,
) -> list[Interval]:
    """Merged intervals supported by >= ``min_replicates`` samples of one condition."""
    union = merge([p for sample in peaks.values() for p in sample])
    index = {sample: PeakIndex(p) for sample, p in peaks.items()}
    groups: dict[str, list[str]] = defaultdict(list)
    for sample, condition in condition_of.items():
        groups[condition].append(sample)
    keep = []
    for chrom, start, end in union:
        for members in groups.values():
            support = sum(index[s].overlaps(chrom, start, end) for s in members)
            if support >= min_replicates:
                keep.append((chrom, start, end))
                break
    return keep


def read_conditions(path: Path) -> dict[str, str]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return {r["sample_id"]: r["condition"] for r in rows}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--samples", type=Path, required=True)
    parser.add_argument("--peaks-pattern", required=True, help="contains {sample}")
    parser.add_argument("--build", required=True, help="genome build for element IDs")
    parser.add_argument("--min-replicates", type=int, default=2)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args(argv)

    condition_of = read_conditions(args.samples)
    peaks = {
        s: read_narrowpeak(Path(args.peaks_pattern.format(sample=s)))
        for s in condition_of
    }
    elements = consensus(peaks, condition_of, args.min_replicates)
    with args.out.open("w") as handle:
        for chrom, start, end in elements:
            handle.write(
                f"{chrom}\t{start}\t{end}\t{args.build}:{chrom}:{start}-{end}\n"
            )
    print(f"{len(elements)} consensus peaks from {len(condition_of)} samples")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
