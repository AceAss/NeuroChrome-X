#!/usr/bin/env python3
"""Motif enrichment of differential elements against a matched background.

Pure-Python engine (PRD section 7: used for the toy gates and targeted JASPAR tests;
HOMER and FIMO/AME remain the engines for genome-wide scans). The foreground is the
set of elements called in ``--direction``; the background is drawn from non-differential
accessible elements (adjusted p at least ``--background-min-padj``), matched on GC
content and accessibility (log2 baseMean) with ``--ratio`` background elements each.

    python motif_enrichment.py --dar dar.tsv --fasta genome.fa --direction up \
        --consensus AP1=TGAGTCA CREB=TGACGTCA --out enrichment.tsv
    python motif_enrichment.py --dar dar.tsv --fasta genome.fa --direction up \
        --meme jaspar.meme --tfs CREB1 MEF2A FOS JUN --out enrichment.tsv
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "src"))

from neurochromx.intervals import Region
from neurochromx.motifs.backgrounds import gc_content, match_background
from neurochromx.motifs.enrichment import motif_enrichment
from neurochromx.motifs.meme import read_meme, select_motifs
from neurochromx.motifs.scan import Motif, motif_from_consensus
from neurochromx.sequences import faidx_sequences

COLUMNS = [
    "motif",
    "direction",
    "n_fg",
    "n_bg",
    "hits_fg",
    "hits_bg",
    "frac_fg",
    "frac_bg",
    "odds_ratio",
    "pvalue",
    "padj",
    "engine",
]


def read_table(path: Path) -> list[dict[str, str]]:
    lines = path.read_text().splitlines()
    head = lines[0].split("\t")
    return [dict(zip(head, line.split("\t"), strict=True)) for line in lines[1:]]


def load_motifs(args: argparse.Namespace) -> list[Motif]:
    motifs: list[Motif] = []
    for item in args.consensus or []:
        name, _, seq = item.partition("=")
        motifs.append(motif_from_consensus(name, seq))
    if args.meme:
        found = read_meme(args.meme)
        motifs.extend(select_motifs(found, args.tfs) if args.tfs else found)
    if not motifs:
        raise SystemExit("no motifs: give --consensus NAME=SEQ or --meme (and --tfs)")
    return motifs


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--dar", type=Path, required=True)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--direction", choices=["up", "down"], required=True)
    parser.add_argument("--consensus", nargs="*", help="NAME=IUPAC consensus")
    parser.add_argument("--meme", type=Path)
    parser.add_argument("--tfs", nargs="*")
    parser.add_argument("--min-rel-score", type=float, default=0.9)
    parser.add_argument("--background-min-padj", type=float, default=0.2)
    parser.add_argument("--ratio", type=int, default=2)
    parser.add_argument("--min-foreground", type=int, default=5)
    parser.add_argument("--seed", type=int, default=20261105)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args(argv)

    rows = read_table(args.dar)
    foreground = [r for r in rows if r["call"] == args.direction]
    pool = [
        r
        for r in rows
        if r["padj"] in ("NA", "") or float(r["padj"]) >= args.background_min_padj
    ]
    header = "\t".join(COLUMNS) + "\n"
    if len(foreground) < args.min_foreground or not pool:
        args.out.write_text(header)
        print(
            f"{len(foreground)} {args.direction} elements: below the minimum; no test"
        )
        return 0

    def regions(rs: list[dict[str, str]]) -> list[Region]:
        return [
            Region(r["chrom"], int(r["start"]), int(r["end"]), r["element_id"])
            for r in rs
        ]

    fg_seq = faidx_sequences(args.fasta, regions(foreground))
    pool_seq = faidx_sequences(args.fasta, regions(pool))

    def covariates(rs: list[dict[str, str]], seqs: list[str]) -> np.ndarray:
        return np.array(
            [
                [gc_content(s), math.log2(float(r["baseMean"]) + 1.0)]
                for r, s in zip(rs, seqs, strict=True)
            ]
        )

    rng = np.random.default_rng(args.seed)
    ratio = max(1, min(args.ratio, len(pool) // len(foreground)))
    picked = match_background(
        covariates(foreground, fg_seq),
        covariates(pool, pool_seq),
        rng,
        ratio=ratio,
        replace=len(pool) < len(foreground) * ratio,
    )
    bg_seq = [pool_seq[i] for i in picked]
    results = motif_enrichment(fg_seq, bg_seq, load_motifs(args), args.min_rel_score)
    with args.out.open("w") as handle:
        handle.write(header)
        for r in results:
            values = [
                r.motif,
                args.direction,
                r.n_fg,
                r.n_bg,
                r.hits_fg,
                r.hits_bg,
                f"{r.frac_fg:.4f}",
                f"{r.frac_bg:.4f}",
                f"{r.odds_ratio:.4f}",
                f"{r.pvalue:.3e}",
                f"{r.padj:.3e}",
                "python",
            ]
            handle.write("\t".join(str(v) for v in values) + "\n")
    print(
        f"{len(foreground)} foreground, {len(bg_seq)} matched background; "
        f"top: {results[0].motif}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
