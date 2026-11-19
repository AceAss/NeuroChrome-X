"""Unit tests for the helper modules and scripts the workflow calls."""

from __future__ import annotations

import importlib.util
import shutil
from pathlib import Path
from types import ModuleType

import pytest

from neurochromx.intervals import Region
from neurochromx.motifs.meme import read_meme, select_motifs
from neurochromx.motifs.scan import relative_score
from neurochromx.sequences import faidx_sequences

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / "workflow" / "scripts" / "py"
TOY_SHEET = ROOT / "tests" / "data" / "toy" / "samples.tsv"
needs_samtools = pytest.mark.skipif(
    shutil.which("samtools") is None, reason="samtools not on PATH"
)


def load_script(name: str) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / f"{name}.py")
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


MEME = """MEME version 4

ALPHABET= ACGT

MOTIF MA0099.3 FOS::JUN
letter-probability matrix: alength= 4 w= 3 nsites= 20 E= 0
0.9 0.03 0.03 0.04
0.03 0.03 0.9 0.04
0.04 0.9 0.03 0.03

MOTIF MA0018.4 CREB1
letter-probability matrix: alength= 4 w= 2 nsites= 20 E= 0
0.03 0.03 0.04 0.9
0.03 0.9 0.04 0.03
"""


def test_read_and_select_meme(tmp_path: Path) -> None:
    path = tmp_path / "m.meme"
    path.write_text(MEME)
    motifs = read_meme(path)
    assert [m.name for m in motifs] == ["FOS::JUN (MA0099.3)", "CREB1 (MA0018.4)"]
    assert [m.length for m in motifs] == [3, 2]
    assert relative_score("TTAGCTT", motifs[0]) == pytest.approx(1.0)  # AGC
    assert [m.name for m in select_motifs(motifs, ["creb1"])] == ["CREB1 (MA0018.4)"]
    assert select_motifs(motifs, ["JUN"]) == []  # whole name must match, not a part


@needs_samtools
def test_faidx_sequences_match_the_fasta(toy_small: Path) -> None:
    fasta = toy_small / "toy.fa"
    genome = "".join(fasta.read_text().split(">chrM")[0].splitlines()[1:])
    regions = [
        Region("chr19", 1000, 1060),
        Region("chr19", 0, 10),
        Region("chrM", 5, 25),
    ]
    got = faidx_sequences(fasta, regions)
    assert got[0] == genome[1000:1060] and got[1] == genome[0:10]
    assert len(got[2]) == 20
    assert faidx_sequences(fasta, []) == []


def test_expressed_genes_filter(tmp_path: Path) -> None:
    script = load_script("expressed_genes")
    counts = tmp_path / "c.tsv"
    counts.write_text(
        "gene_id\tgene_name\ts1\ts2\ts3\n"
        "g1\tA\t50\t60\t2\n"  # expressed in 2 of 3
        "g2\tB\t50\t1\t2\n"  # only 1
        "g3\tC\t10\t10\t10\n"  # exactly at the threshold
    )
    assert script.expressed_genes(counts, ["s1", "s2", "s3"], 10, 2) == ["g1", "g3"]
    with pytest.raises(SystemExit):
        script.expressed_genes(counts, ["s1", "missing"], 10, 1)


def test_consensus_requires_replicates_within_a_condition() -> None:
    script = load_script("consensus_peaks")
    peaks = {
        "c1": [("chr1", 100, 200), ("chr1", 500, 600)],
        "c2": [("chr1", 150, 250)],
        "s1": [("chr1", 500, 650), ("chr1", 900, 1000)],
        "s2": [("chr1", 520, 640)],
        "s3": [("chr1", 5000, 5100)],
    }
    condition = {"c1": "ctrl", "c2": "ctrl", "s1": "stim", "s2": "stim", "s3": "stim"}
    kept = script.consensus(peaks, condition, min_replicates=2)
    assert kept == [("chr1", 100, 250), ("chr1", 500, 650)]
    assert script.consensus(peaks, condition, min_replicates=1)[-1] == (
        "chr1",
        5000,
        5100,
    )


def test_group_sheet_script(tmp_path: Path) -> None:
    script = load_script("write_group_sheet")
    out = tmp_path / "g.tsv"
    assert (
        script.main(
            [
                "--samples",
                str(TOY_SHEET),
                "--study",
                "toy",
                "--assay",
                "atac",
                "--out",
                str(out),
            ]
        )
        == 0
    )
    lines = out.read_text().splitlines()
    assert lines[0].split("\t")[:3] == ["sample_id", "condition", "batch"]
    assert len(lines) == 7 and lines[1].split("\t")[1] == "control"
    with pytest.raises(SystemExit):
        script.main(
            [
                "--samples",
                str(TOY_SHEET),
                "--study",
                "nope",
                "--assay",
                "atac",
                "--out",
                str(out),
            ]
        )


def test_annotate_elements_script(toy_small: Path, tmp_path: Path) -> None:
    script = load_script("annotate_elements")
    truth = (toy_small / "truth_regions.tsv").read_text().splitlines()[1:]
    bed = tmp_path / "elements.bed"
    bed.write_text(
        "".join(
            f"chr19\t{r.split(chr(9))[2]}\t{r.split(chr(9))[3]}\t{r.split(chr(9))[0]}\n"
            for r in truth
        )
    )
    expressed = tmp_path / "expressed.txt"
    expressed.write_text("TOYG0001\n")
    out = tmp_path / "ann.tsv"
    args = [
        "--elements",
        str(bed),
        "--gtf",
        str(toy_small / "toy.gtf"),
        "--ccres",
        str(toy_small / "toy_ccres.bed"),
    ]
    assert script.main([*args, "--out", str(out)]) == 0
    rows = [line.split("\t") for line in out.read_text().splitlines()[1:]]
    assert len(rows) == 400 and sum(r[4] == "promoter" for r in rows) == 60
    assert all(r[9] != "NA" for r in rows)
    assert script.main([*args, "--expressed", str(expressed), "--out", str(out)]) == 0
    linked = {line.split("\t")[9] for line in out.read_text().splitlines()[1:]}
    assert linked <= {"TOYG0001", "NA"}  # only the one expressed gene can be a target


@needs_samtools
def test_motif_enrichment_script_on_synthetic_dar(
    toy_small: Path, tmp_path: Path
) -> None:
    script = load_script("motif_enrichment")
    rows = (toy_small / "truth_regions.tsv").read_text().splitlines()[1:]
    dar = tmp_path / "dar.tsv"
    lines = ["element_id\tchrom\tstart\tend\tbaseMean\tlog2FC\tpadj\tcall"]
    for row in rows:
        f = row.split("\t")
        status = f[7]
        call = "up" if status == "gained" else "down" if status == "lost" else "ns"
        padj = "1e-6" if call != "ns" else "0.9"
        lines.append(f"{f[0]}\tchr19\t{f[2]}\t{f[3]}\t100\t{f[6]}\t{padj}\t{call}")
    dar.write_text("\n".join(lines) + "\n")
    out = tmp_path / "up.tsv"
    base = [
        "--dar",
        str(dar),
        "--fasta",
        str(toy_small / "toy.fa"),
        "--consensus",
        "AP1=TGAGTCA",
        "CREB=TGACGTCA",
        "MEF2=CTAAAAATAG",
    ]
    assert script.main([*base, "--direction", "up", "--out", str(out)]) == 0
    table = [line.split("\t") for line in out.read_text().splitlines()]
    assert (
        table[0][0] == "motif" and table[1][0] == "AP1" and float(table[1][-2]) < 1e-6
    )
    down = (
        tmp_path / "down.tsv"
    )  # 12 lost regions, none with planted motifs: no enrichment
    assert script.main([*base, "--direction", "down", "--out", str(down)]) == 0
    assert (
        min(float(r.split("\t")[-2]) for r in down.read_text().splitlines()[1:]) > 0.05
    )
    empty = tmp_path / "empty.tsv"
    assert (
        script.main(
            [*base, "--direction", "up", "--min-foreground", "999", "--out", str(empty)]
        )
        == 0
    )
    assert len(empty.read_text().splitlines()) == 1  # header only
