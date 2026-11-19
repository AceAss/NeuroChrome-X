"""Tests for tests/data/toy/make_toy.py (TC-004, built by Claude)."""

from __future__ import annotations

import gzip
import importlib.util
import json
from pathlib import Path
from types import ModuleType

import pytest

from neurochromx.io.samplesheet import load_samplesheet

TOY_DIR = Path(__file__).resolve().parents[1] / "data" / "toy"


def _load_generator() -> ModuleType:
    spec = importlib.util.spec_from_file_location("make_toy", TOY_DIR / "make_toy.py")
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


make_toy = _load_generator()
DEPTH = 6000


@pytest.fixture(scope="module")
def toy(tmp_path_factory: pytest.TempPathFactory) -> tuple[Path, dict]:
    out = tmp_path_factory.mktemp("toy")
    return out, make_toy.generate(out, seed=20261105, depth=DEPTH)


def read_table(path: Path) -> list[dict[str, str]]:
    lines = path.read_text().splitlines()
    header = lines[0].split("\t")
    return [dict(zip(header, line.split("\t"), strict=True)) for line in lines[1:]]


def read_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, list[str]] = {}
    name = ""
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:]
            seqs[name] = []
        else:
            seqs[name].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def read_fastq(path: Path) -> list[tuple[str, str]]:
    with gzip.open(path, "rt") as handle:
        lines = handle.read().splitlines()
    return [(lines[i], lines[i + 1]) for i in range(0, len(lines), 4)]


def has_motif(seq: str, motif: str) -> bool:
    return motif in seq or make_toy.revcomp(motif) in seq


def test_design_counts(toy: tuple[Path, dict]) -> None:
    _, manifest = toy
    assert manifest["truth"] == {
        "regions": 400,
        "promoter": 60,
        "distal": 340,
        "gained": 48,
        "lost": 12,
        "unchanged": 340,
        "genes": 60,
    }


def test_named_iegs_present_and_gained(toy: tuple[Path, dict]) -> None:
    out, _ = toy
    genes = {r["gene_name"]: r for r in read_table(out / "truth_genes.tsv")}
    regions = read_table(out / "truth_regions.tsv")
    for name in make_toy.IEG_GENES:
        assert genes[name]["status"] == "gained"
        promoter = [
            r
            for r in regions
            if r["class"] == "promoter"
            and r["linked_gene_id"] == genes[name]["gene_id"]
        ]
        assert len(promoter) == 1 and promoter[0]["status"] == "gained"
    assert sum(1 for g in genes if g.startswith("ToyGene")) == 54


def test_lost_regions_are_distal_and_have_no_planted_motifs(
    toy: tuple[Path, dict],
) -> None:
    out, _ = toy
    for region in read_table(out / "truth_regions.tsv"):
        if region["status"] == "lost":
            assert region["class"] == "distal"
            assert region["motifs"] == "NA"
            assert float(region["true_log2fc"]) < -1.2 + 1e-9


def test_planted_motif_counts_are_exact(toy: tuple[Path, dict]) -> None:
    out, _ = toy
    regions = read_table(out / "truth_regions.tsv")
    gained = [r for r in regions if r["status"] == "gained"]
    unchanged = [r for r in regions if r["status"] == "unchanged"]
    for name, share in make_toy.PLANT_GAINED.items():
        planted = sum(name in r["motifs"].split(",") for r in gained)
        assert planted == round(share * len(gained))
    for name in make_toy.MOTIFS:
        planted = sum(name in r["motifs"].split(",") for r in unchanged)
        assert planted == round(make_toy.PLANT_UNCHANGED * len(unchanged))


def test_ap1_in_sequence_gained_vs_unchanged(toy: tuple[Path, dict]) -> None:
    out, _ = toy
    seq = read_fasta(out / "toy.fa")["chr19"]
    regions = read_table(out / "truth_regions.tsv")

    def share(status: str) -> float:
        rows = [r for r in regions if r["status"] == status]
        hits = sum(
            has_motif(seq[int(r["start"]) : int(r["end"])], make_toy.MOTIFS["AP1"])
            for r in rows
        )
        return hits / len(rows)

    assert share("gained") >= 0.75
    assert share("unchanged") <= 0.20


def test_fastq_files_pairs_and_read_length(toy: tuple[Path, dict]) -> None:
    out, manifest = toy
    for sample, info in manifest["atac_samples"].items():
        r1 = read_fastq(out / "fastq" / f"{sample}_R1.fastq.gz")
        r2 = read_fastq(out / "fastq" / f"{sample}_R2.fastq.gz")
        assert len(r1) == len(r2) == info["pairs"]
        assert all(len(seq) == make_toy.READ_LEN for _, seq in r1[:200] + r2[:200])
        assert [h.rsplit("/", 1)[0] for h, _ in r1[:50]] == [
            h.rsplit("/", 1)[0] for h, _ in r2[:50]
        ]
        assert set("".join(seq for _, seq in r1[:200])) <= set("ACGT")


def test_read_mix_matches_design(toy: tuple[Path, dict]) -> None:
    _, manifest = toy
    for info in manifest["atac_samples"].values():
        total = info["pairs"]
        assert 0.10 <= info["n_chrM"] / total <= 0.14
        assert 0.08 <= info["n_duplicates"] / total <= 0.12
        assert 0.03 <= info["n_blacklist"] / total <= 0.05
        assert (
            info["n_chrM"]
            + info["n_blacklist"]
            + info["n_region"]
            + info["n_background"]
            == total
        )


def test_rna_counts_and_atac_files_agree_with_toy_sample_sheet(
    toy: tuple[Path, dict],
) -> None:
    out, _ = toy
    sheet = load_samplesheet(TOY_DIR / "samples.tsv")
    rna_ids = [r.sample_id for r in sheet.by_assay("rna")]
    header = (out / "rna_counts.tsv").read_text().splitlines()[0].split("\t")
    assert header[2:] == rna_ids
    for record in sheet.by_assay("atac"):
        assert record.path_or_uri is not None
        for part in record.path_or_uri.split(";"):
            assert (out / "fastq" / Path(part).name).exists()


def test_rna_counts_follow_truth(toy: tuple[Path, dict]) -> None:
    out, _ = toy
    truth = {
        r["gene_id"]: float(r["true_log2fc"])
        for r in read_table(out / "truth_genes.tsv")
    }
    gained = [g for g, v in truth.items() if v > 1.4]
    rows = {r["gene_id"]: r for r in read_table(out / "rna_counts.tsv")}
    for gid in gained:
        ctrl = sum(int(rows[gid][f"toy_rna_ctrl_r{i}"]) for i in (1, 2, 3))
        stim = sum(int(rows[gid][f"toy_rna_stim_r{i}"]) for i in (1, 2, 3))
        assert stim > ctrl


def test_ccre_and_blacklist_files(toy: tuple[Path, dict]) -> None:
    out, _ = toy
    ccres = [
        line.split("\t") for line in (out / "toy_ccres.bed").read_text().splitlines()
    ]
    assert sum(c[4] == "PLS" for c in ccres) == 60
    assert sum(c[4] == "dELS" for c in ccres) == round(make_toy.CCRE_FRACTION * 340)
    blacklist = (out / "toy_blacklist.bed").read_text().splitlines()
    assert len(blacklist) == 2
    assert (out / "toy.chrom.sizes").read_text() == "chr19\t3000000\nchrM\t16299\n"


def test_gtf_has_six_ieg_names(toy: tuple[Path, dict]) -> None:
    out, _ = toy
    text = (out / "toy.gtf").read_text()
    for name in make_toy.IEG_GENES:
        assert f'gene_name "{name}"' in text
    assert text.count("\tgene\t") == 60


def test_determinism_same_seed_same_hashes(tmp_path: Path) -> None:
    a = make_toy.generate(tmp_path / "a", seed=7, depth=1500)
    b = make_toy.generate(tmp_path / "b", seed=7, depth=1500)
    c = make_toy.generate(tmp_path / "c", seed=8, depth=1500)
    assert a["sha256_uncompressed"] == b["sha256_uncompressed"]
    assert a["sha256_uncompressed"] != c["sha256_uncompressed"]
    manifest = json.loads((tmp_path / "a" / "MANIFEST.json").read_text())
    assert manifest["seed"] == 7
    assert (tmp_path / "a" / "fastq" / "toy_atac_ctrl_r1_R1.fastq.gz").read_bytes() == (
        tmp_path / "b" / "fastq" / "toy_atac_ctrl_r1_R1.fastq.gz"
    ).read_bytes()
