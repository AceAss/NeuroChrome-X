#!/usr/bin/env python3
"""List expressed genes from a counts matrix (config.rna.gene_filter, PRD section 7).

    python expressed_genes.py --counts rna_counts.tsv --samples s1 s2 s3 \
        --min-count 10 --min-samples 2 --out expressed_genes.txt

The matrix has ``gene_id`` in the first column; sample columns are named by sample ID.
"""

from __future__ import annotations

import argparse
from pathlib import Path


def expressed_genes(
    counts: Path, samples: list[str], min_count: int, min_samples: int
) -> list[str]:
    lines = counts.read_text().splitlines()
    header = lines[0].split("\t")
    missing = [s for s in samples if s not in header]
    if missing:
        raise SystemExit(f"samples not in {counts}: {', '.join(missing)}")
    columns = [header.index(s) for s in samples]
    keep = []
    for line in lines[1:]:
        fields = line.split("\t")
        if sum(int(float(fields[c])) >= min_count for c in columns) >= min_samples:
            keep.append(fields[0])
    return keep


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--counts", type=Path, required=True)
    parser.add_argument("--samples", nargs="+", required=True)
    parser.add_argument("--min-count", type=int, default=10)
    parser.add_argument("--min-samples", type=int, default=2)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args(argv)
    genes = expressed_genes(args.counts, args.samples, args.min_count, args.min_samples)
    args.out.write_text("".join(f"{g}\n" for g in genes))
    print(f"{len(genes)} expressed genes")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
