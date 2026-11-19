"""Tests for intervals, promoters, enhancer tiers and gene linking (toy truth)."""

from __future__ import annotations

from pathlib import Path

import pytest

from neurochromx.annotate.enhancers import TssIndex, assign_elements
from neurochromx.annotate.linking import link_nearest_expressed_tss
from neurochromx.annotate.promoters import (
    Gene,
    promoter_regions,
    promoter_window,
    read_genes,
)
from neurochromx.intervals import IntervalIndex, Region, read_bed


def truth_rows(out: Path) -> list[dict[str, str]]:
    lines = (out / "truth_regions.tsv").read_text().splitlines()
    head = lines[0].split("\t")
    return [dict(zip(head, line.split("\t"), strict=True)) for line in lines[1:]]


def elements(out: Path) -> list[Region]:
    return [
        Region(r["chrom"], int(r["start"]), int(r["end"]), r["region_id"])
        for r in truth_rows(out)
    ]


# --------------------------------------------------------------------------- #
def test_interval_index_overlaps_and_bed_roundtrip(tmp_path: Path) -> None:
    bed = tmp_path / "x.bed"
    bed.write_text("# c\nchr1\t10\t20\ta\tpELS\nchr1\t15\t40\tb\tdELS\nchr2\t0\t5\n")
    regions = read_bed(bed)
    assert [r.name for r in regions] == ["a", "b", ""]
    assert regions[1].label == "dELS"
    index = IntervalIndex(regions)
    assert {r.name for r in index.overlapping("chr1", 19, 21)} == {"a", "b"}
    assert not index.any_overlap("chr1", 40, 50)  # half-open: 40 is outside [15, 40)
    assert index.any_overlap("chr2", 4, 5)
    assert index.overlapping("chrX", 0, 100) == []
    with pytest.raises(ValueError):
        Region("chr1", 5, 4)


def test_read_genes_and_strand_aware_promoters(toy_small: Path) -> None:
    genes = read_genes(toy_small / "toy.gtf")
    assert len(genes) == 60
    names = {g.gene_name for g in genes}
    assert {"Fos", "Arc", "Npas4", "Egr1", "Junb", "Nr4a1"} <= names
    plus = Gene("g1", "G1", "chr19", 1000, 2000, "+", "protein_coding")
    minus = Gene("g2", "G2", "chr19", 1000, 2000, "-", "protein_coding")
    assert plus.tss == 1000 and minus.tss == 1999
    assert promoter_window(plus, 300, 200) == (700, 1200)
    assert promoter_window(minus, 300, 200) == (1799, 2299)
    assert read_genes(toy_small / "toy.gtf", gene_types=["lncRNA"]) == []


def test_toy_promoter_regions_contain_their_tss(toy_small: Path) -> None:
    genes = read_genes(toy_small / "toy.gtf")
    truth = {
        r["linked_gene_id"]: r
        for r in truth_rows(toy_small)
        if r["class"] == "promoter"
    }
    for region, gene in zip(promoter_regions(genes, 300, 200), genes, strict=True):
        assert region.start <= gene.tss < region.end
        row = truth[gene.gene_id]
        assert (
            abs(region.start - int(row["start"])) <= 1
            and abs(region.end - int(row["end"])) <= 1
        )


def test_tss_index_distances() -> None:
    genes = [
        Gene("a", "A", "chr1", 1000, 2000, "+", "p"),  # TSS 1000
        Gene("b", "B", "chr1", 5000, 6000, "-", "p"),  # TSS 5999
    ]
    index = TssIndex(genes)
    assert index.distance("chr1", 990, 1010) == 0
    assert index.distance("chr1", 1100, 1200) == 100
    assert index.distance("chr1", 3000, 3100) == 2000  # nearer to the TSS at 1000
    assert index.distance("chr1", 4000, 4100) == 1900  # nearer to the TSS at 5999
    assert index.distance("chrZ", 0, 10) is None


def test_tiers_match_toy_truth(toy_small: Path) -> None:
    genes = read_genes(toy_small / "toy.gtf")
    ccres = read_bed(toy_small / "toy_ccres.bed")
    els = elements(toy_small)
    rows = {r["region_id"]: r for r in truth_rows(toy_small)}
    annotated = assign_elements(els, promoter_regions(genes, 1000, 500), genes, ccres)
    by_class = {
        c: [a for a in annotated if a.cls == c]
        for c in ("promoter", "enhancer", "other")
    }
    assert len(by_class["promoter"]) == 60
    assert len(by_class["enhancer"]) == 340 and not by_class["other"]
    for a in annotated:
        truth = rows[a.element_id]
        assert (truth["class"] == "promoter") == (a.cls == "promoter")
        assert not a.tier2  # no histone peaks supplied
        if a.cls == "enhancer":
            assert a.tier1 == (truth["in_ccre"] == "true")
            assert a.tier3 and a.tss_distance > 2000
        else:
            assert not a.tier1  # promoters carry the PLS class, not pELS/dELS
    assert sum(a.tier1 for a in by_class["enhancer"]) == round(0.7 * 340)


def test_tier2_needs_k27ac_k4me1_and_no_k4me3() -> None:
    el = [
        Region("c", 100, 200, "e1"),
        Region("c", 300, 400, "e2"),
        Region("c", 500, 600, "e3"),
    ]
    far = [Gene("g", "G", "c", 100_000, 101_000, "+", "p")]
    k27 = [Region("c", 90, 210), Region("c", 290, 410), Region("c", 490, 610)]
    k4me1 = [Region("c", 95, 205), Region("c", 295, 405)]
    k4me3 = [Region("c", 300, 400)]
    out = assign_elements(el, [], far, h3k27ac=k27, h3k4me1=k4me1, h3k4me3=k4me3)
    assert [a.tier2 for a in out] == [True, False, False]
    assert not any(a.tier2 for a in assign_elements(el, [], far, h3k27ac=k27))


def test_links_match_truth_and_respect_expression(toy_small: Path) -> None:
    genes = read_genes(toy_small / "toy.gtf")
    els = elements(toy_small)
    rows = {r["region_id"]: r for r in truth_rows(toy_small)}
    links = link_nearest_expressed_tss(els, genes, window=100_000)
    for element in els:
        link = links[element.name]
        expected = rows[element.name]["linked_gene_id"]
        assert (link.gene_id if link else "NA") == expected
    target = next(e for e in els if links[e.name] is not None)
    nearest = links[target.name].gene_id
    others = {g.gene_id for g in genes} - {nearest}
    relinked = link_nearest_expressed_tss(
        [target], genes, expressed=others, window=100_000
    )[target.name]
    assert relinked is None or relinked.gene_id != nearest
    assert (
        link_nearest_expressed_tss([target], genes, expressed=set())[target.name]
        is None
    )
    assert links[target.name].distance > 0
    assert link_nearest_expressed_tss([target], genes, window=0)[target.name] is None
