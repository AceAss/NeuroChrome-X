#!/usr/bin/env python3
"""Count properly paired reads per consensus peak with ``bedtools multicov``.

    python count_peaks.py --bed consensus.bed --bams a.bam b.bam \
        --samples A B --out counts.tsv

Output columns: element_id, chrom, start, end, then one column per sample (in the
order given). BAM files must be coordinate-sorted and indexed.
"""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument(
        "--bed", type=Path, required=True, help="BED with IDs in column 4"
    )
    parser.add_argument("--bams", nargs="+", required=True)
    parser.add_argument("--samples", nargs="+", required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args(argv)
    if len(args.bams) != len(args.samples):
        parser.error("--bams and --samples must have the same length")

    result = subprocess.run(
        ["bedtools", "multicov", "-p", "-bams", *args.bams, "-bed", str(args.bed)],
        capture_output=True,
        text=True,
        check=True,
    )
    with args.out.open("w") as handle:
        handle.write(
            "\t".join(["element_id", "chrom", "start", "end", *args.samples]) + "\n"
        )
        for line in result.stdout.splitlines():
            chrom, start, end, name, *counts = line.split("\t")
            handle.write("\t".join([name, chrom, start, end, *counts]) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
