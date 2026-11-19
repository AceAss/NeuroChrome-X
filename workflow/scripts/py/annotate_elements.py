#!/usr/bin/env python3
"""Annotate consensus elements: class, enhancer tiers, TSS distance, gene link.

PRD sections 8 and 9, ADR-007. Writes ``regulatory_elements.tsv``.

    python annotate_elements.py --elements consensus.bed --gtf genes.gtf \
        --ccres ccres.bed --expressed expressed_genes.txt --out regulatory_elements.tsv
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "src"))

from neurochromx.annotate.enhancers import assign_elements
from neurochromx.annotate.linking import link_nearest_expressed_tss
from neurochromx.annotate.promoters import promoter_regions, read_genes
from neurochromx.intervals import read_bed

COLUMNS = [
    "element_id",
    "chrom",
    "start",
    "end",
    "class",
    "tier1",
    "tier2",
    "tier3",
    "tss_distance",
    "linked_gene_id",
    "linked_gene_name",
    "link_method",
    "link_distance",
    "sources",
]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--elements", type=Path, required=True)
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--ccres", type=Path)
    parser.add_argument("--h3k27ac", type=Path)
    parser.add_argument("--h3k4me1", type=Path)
    parser.add_argument("--h3k4me3", type=Path)
    parser.add_argument("--expressed", type=Path, help="file with one gene_id per line")
    parser.add_argument("--gene-type", nargs="*", default=["protein_coding"])
    parser.add_argument("--promoter-upstream", type=int, default=1000)
    parser.add_argument("--promoter-downstream", type=int, default=500)
    parser.add_argument("--ccre-classes", nargs="*", default=["pELS", "dELS"])
    parser.add_argument("--tier3-min-tss-distance", type=int, default=2000)
    parser.add_argument("--link-window", type=int, default=100_000)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args(argv)

    elements = read_bed(args.elements)
    genes = read_genes(args.gtf, gene_types=args.gene_type or None)
    promoters = promoter_regions(
        genes, args.promoter_upstream, args.promoter_downstream
    )

    def optional(path: Path | None):
        return read_bed(path) if path else None

    annotated = assign_elements(
        elements,
        promoters,
        genes,
        ccres=read_bed(args.ccres) if args.ccres else (),
        h3k27ac=optional(args.h3k27ac),
        h3k4me1=optional(args.h3k4me1),
        h3k4me3=optional(args.h3k4me3),
        ccre_classes=args.ccre_classes,
        tier3_min_tss_distance=args.tier3_min_tss_distance,
    )
    expressed = None
    if args.expressed:
        expressed = {
            x.strip() for x in args.expressed.read_text().splitlines() if x.strip()
        }
    links = link_nearest_expressed_tss(elements, genes, expressed, args.link_window)
    sources = (
        ["gtf"]
        + (["ccres"] if args.ccres else [])
        + (["expression"] if expressed else [])
    )

    with args.out.open("w") as handle:
        handle.write("\t".join(COLUMNS) + "\n")
        for a in annotated:
            link = links[a.element_id]
            row = [
                a.element_id,
                a.chrom,
                a.start,
                a.end,
                a.cls,
                str(a.tier1).lower(),
                str(a.tier2).lower(),
                str(a.tier3).lower(),
                a.tss_distance,
                link.gene_id if link else "NA",
                link.gene_name if link else "NA",
                link.method if link else "NA",
                link.distance if link else "NA",
                ";".join(sources),
            ]
            handle.write("\t".join(str(x) for x in row) + "\n")
    counts = {
        c: sum(a.cls == c for a in annotated) for c in ("promoter", "enhancer", "other")
    }
    print(f"annotated {len(annotated)} elements: {counts}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
