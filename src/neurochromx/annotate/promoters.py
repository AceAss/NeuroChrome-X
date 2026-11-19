"""Genes, TSSs and promoter windows from a GTF (ADR-006, PRD section 7)."""

from __future__ import annotations

import gzip
import re
from collections.abc import Collection, Iterable
from dataclasses import dataclass
from pathlib import Path

from neurochromx.intervals import Region

_ATTRIBUTE = re.compile(r'(\w+) "([^"]*)"')


@dataclass(frozen=True)
class Gene:
    """A gene with 0-based, half-open coordinates."""

    gene_id: str
    gene_name: str
    chrom: str
    start: int
    end: int
    strand: str
    gene_type: str

    @property
    def tss(self) -> int:
        """First transcribed base: ``start`` on +, ``end - 1`` on -."""
        return self.start if self.strand == "+" else self.end - 1


def read_genes(
    gtf_path: str | Path, gene_types: Collection[str] | None = ("protein_coding",)
) -> list[Gene]:
    """Parse the ``gene`` lines of a GTF (1-based, inclusive) into :class:`Gene`.

    ``gene_types=None`` keeps every gene. Genes are returned sorted by position.
    """
    path = Path(gtf_path)
    opener = gzip.open if path.suffix == ".gz" else open
    genes: list[Gene] = []
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            attrs = dict(_ATTRIBUTE.findall(f[8]))
            gene_type = attrs.get("gene_type") or attrs.get("gene_biotype", "")
            if gene_types is not None and gene_type not in gene_types:
                continue
            genes.append(
                Gene(
                    gene_id=attrs["gene_id"],
                    gene_name=attrs.get("gene_name", attrs["gene_id"]),
                    chrom=f[0],
                    start=int(f[3]) - 1,
                    end=int(f[4]),
                    strand=f[6],
                    gene_type=gene_type,
                )
            )
    return sorted(genes, key=lambda g: (g.chrom, g.start))


def promoter_window(
    gene: Gene, upstream: int = 1000, downstream: int = 500
) -> tuple[int, int]:
    """Strand-aware window: ``upstream`` bases before the TSS, ``downstream`` after."""
    if gene.strand == "+":
        return gene.tss - upstream, gene.tss + downstream
    return gene.tss - downstream, gene.tss + upstream


def promoter_regions(
    genes: Iterable[Gene], upstream: int = 1000, downstream: int = 500
) -> list[Region]:
    """Promoter windows as regions named by gene ID (config.annotation defaults)."""
    out = []
    for gene in genes:
        start, end = promoter_window(gene, upstream, downstream)
        out.append(Region(gene.chrom, max(0, start), end, gene.gene_id, gene.gene_name))
    return out
