"""End-to-end check of the ATAC branch on the synthetic toy dataset (PRD section 11).

Runs the real tools (bowtie2, samtools, bedtools, MACS3, DESeq2 through
workflow/scripts/R/deseq2_dar.R) on the generated toy reads and compares the results
with the planted truth. This is the P1 positive and negative control gate in
miniature. The test is skipped when a tool is missing; the conda environments from
TC-002 provide all of them. It takes a minute or two.
"""

from __future__ import annotations

import importlib.util
import itertools
import shutil
import statistics
import subprocess
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.integration

ROOT = Path(__file__).resolve().parents[2]
TOY = ROOT / "tests" / "data" / "toy"
SCRIPTS_PY = ROOT / "workflow" / "scripts" / "py"
DESEQ2 = ROOT / "workflow" / "scripts" / "R" / "deseq2_dar.R"
TOOLS = ("bowtie2", "bowtie2-build", "samtools", "bedtools", "macs3", "Rscript")
SAMPLES = [f"toy_atac_{tag}_r{rep}" for tag in ("ctrl", "stim") for rep in (1, 2, 3)]
CONDITION = dict(zip(SAMPLES, ["control"] * 3 + ["stimulated"] * 3, strict=True))

missing = [tool for tool in TOOLS if shutil.which(tool) is None]
if missing:
    pytest.skip(f"tools not on PATH: {', '.join(missing)}", allow_module_level=True)
if subprocess.run(["Rscript", "-e", "library(DESeq2)"], capture_output=True).returncode:
    pytest.skip("R package DESeq2 is not installed", allow_module_level=True)


def sh(command: str) -> str:
    done = subprocess.run(
        command,
        shell=True,
        check=True,
        capture_output=True,
        text=True,
        executable="/bin/bash",
    )
    return done.stdout + done.stderr


def load_generator():
    spec = importlib.util.spec_from_file_location("make_toy", TOY / "make_toy.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def count_pairs(args: str) -> int:
    """Read pairs (first mates) in a BAM, optionally restricted to a region."""
    return int(sh(f"samtools view -c -f 64 {args}").split()[0])


def overlap(a: tuple[int, int], b: tuple[int, int]) -> bool:
    return a[0] < b[1] and b[0] < a[1]


@pytest.fixture(scope="module")
def chain(tmp_path_factory: pytest.TempPathFactory) -> dict:
    work = tmp_path_factory.mktemp("toy_chain")
    data = work / "data"
    load_generator().generate(data, seed=20261105, depth=70_000)
    bam = work / "bam"
    peaks = work / "peaks"
    bam.mkdir()
    peaks.mkdir()
    sh(f"bowtie2-build -q {data}/toy.fa {work}/toy")

    qc: dict[str, dict] = {}
    for s in SAMPLES:
        b = bam / s
        log = sh(
            f"bowtie2 -x {work}/toy -1 {data}/fastq/{s}_R1.fastq.gz "
            f"-2 {data}/fastq/{s}_R2.fastq.gz --very-sensitive -X 2000 --no-mixed "
            f"--no-discordant -p 2 2> {b}.log | samtools view -b -o {b}.raw.bam - ; "
            f"grep 'overall alignment rate' {b}.log"
        )
        sh(
            f"samtools view -b -q 30 -f 2 -F 1804 {b}.raw.bam | "
            f"samtools sort -n -o {b}.n.bam - ; "
            f"samtools fixmate -m {b}.n.bam - | samtools sort -o {b}.s.bam - ; "
            f"samtools markdup -r {b}.s.bam {b}.dd.bam ; samtools index {b}.dd.bam ; "
            f"samtools view -b {b}.dd.bam chr19 | "
            f"bedtools intersect -v -abam - -b {data}/toy_blacklist.bed > {b}.bam ; "
            f"samtools index {b}.bam"
        )
        sh(
            f"macs3 callpeak -t {b}.bam -f BAMPE -g 3e6 -q 0.01 --keep-dup all "
            f"-n {s} --outdir {peaks}"
        )
        qc[s] = {
            "rate": float(log.split("%")[0].split()[-1]),
            "q30": count_pairs(f"{b}.n.bam"),
            "dedup": count_pairs(f"{b}.dd.bam"),
            "chrM_dedup": count_pairs(f"{b}.dd.bam chrM"),
            "final": count_pairs(f"{b}.bam"),
            "chrM_final": int(sh(f"samtools view -c {b}.bam chrM").split()[0]),
            "blacklist_final": int(
                sh(
                    f"bedtools intersect -u -abam {b}.bam "
                    f"-b {data}/toy_blacklist.bed | samtools view -c"
                ).split()[0]
            ),
        }

    sheet = work / "atac_samples.tsv"
    sheet.write_text(
        "sample_id\tcondition\tbatch\n"
        + "".join(f"{s}\t{CONDITION[s]}\tb1\n" for s in SAMPLES)
    )
    sh(
        f"{sys.executable} {SCRIPTS_PY}/consensus_peaks.py --samples {sheet} "
        f"--peaks-pattern '{peaks}/{{sample}}_peaks.narrowPeak' --build toy "
        f"--min-replicates 2 --out {work}/consensus.bed"
    )
    sh(
        f"{sys.executable} {SCRIPTS_PY}/count_peaks.py --bed {work}/consensus.bed "
        f"--bams {' '.join(f'{bam}/{s}.bam' for s in SAMPLES)} "
        f"--samples {' '.join(SAMPLES)} --out {work}/counts.tsv"
    )
    return {"work": work, "data": data, "qc": qc, "sheet": sheet}


def read_tsv(path: Path) -> list[dict[str, str]]:
    lines = path.read_text().splitlines()
    head = lines[0].split("\t")
    return [dict(zip(head, line.split("\t"), strict=True)) for line in lines[1:]]


def run_deseq2(chain: dict, sheet: Path, out: Path) -> list[dict[str, str]]:
    sh(
        f"Rscript {DESEQ2} --counts {chain['work']}/counts.tsv --samples {sheet} "
        f"--treatment stimulated --control control --shrink none --out {out}"
    )
    return read_tsv(out)


def classify(
    chain: dict, rows: list[dict[str, str]]
) -> list[tuple[dict[str, str], dict]]:
    truth = [
        {
            "iv": (int(r["start"]), int(r["end"])),
            "status": r["status"],
            "cls": r["class"],
            "lfc": float(r["true_log2fc"]),
        }
        for r in read_tsv(chain["data"] / "truth_regions.tsv")
    ]
    out = []
    for row in rows:
        hits = [
            t for t in truth if overlap((int(row["start"]), int(row["end"])), t["iv"])
        ]
        if len(hits) == 1:
            out.append((row, hits[0]))
    return out


def test_alignment_and_filtering_qc(chain: dict) -> None:
    for s, q in chain["qc"].items():
        assert q["rate"] >= 99.0, s
        assert 0.08 <= 1 - q["dedup"] / q["q30"] <= 0.12, f"{s}: duplicate fraction"
        assert 0.10 <= q["chrM_dedup"] / q["dedup"] <= 0.14, f"{s}: chrM share"
        assert q["chrM_final"] == 0 and q["blacklist_final"] == 0, s


def test_peak_recall_matches_planted_regions(chain: dict) -> None:
    consensus = [
        (int(a), int(b))
        for _, a, b, _ in (
            line.split("\t")
            for line in (chain["work"] / "consensus.bed").read_text().splitlines()
        )
    ]
    by_status: dict[str, list[bool]] = {"gained": [], "lost": [], "unchanged": []}
    for r in read_tsv(chain["data"] / "truth_regions.tsv"):
        region = (int(r["start"]), int(r["end"]))
        by_status[r["status"]].append(any(overlap(region, c) for c in consensus))
    assert sum(by_status["gained"]) >= 46
    assert sum(by_status["lost"]) >= 8
    assert sum(by_status["unchanged"]) >= 323


def test_deseq2_recovers_truth_with_controlled_false_positives(chain: dict) -> None:
    rows = run_deseq2(chain, chain["sheet"], chain["work"] / "dar.tsv")
    paired = classify(chain, rows)
    by_status: dict[str, list[bool]] = {"gained": [], "lost": [], "unchanged": []}
    wrong_direction = 0
    for row, truth in paired:
        called = row["call"] != "ns"
        by_status[truth["status"]].append(called)
        if (
            called
            and truth["status"] != "unchanged"
            and (row["call"] == "up") != (truth["lfc"] > 0)
        ):
            wrong_direction += 1
    n_calls = sum(sum(v) for v in by_status.values())
    false_calls = sum(by_status["unchanged"])
    assert sum(by_status["gained"]) / len(by_status["gained"]) >= 0.90
    assert sum(by_status["lost"]) / len(by_status["lost"]) >= 0.80
    assert false_calls / len(by_status["unchanged"]) <= 0.05
    assert false_calls / n_calls <= 0.10
    assert wrong_direction == 0
    median_unchanged = statistics.median(
        float(row["log2FC_mle"]) for row, t in paired if t["status"] == "unchanged"
    )
    assert (
        abs(median_unchanged) <= 0.15
    ), "normalisation leaves a shift in unchanged peaks"
    ieg_promoters = [
        (r, t) for r, t in paired if t["cls"] == "promoter" and t["status"] == "gained"
    ]
    assert len(ieg_promoters) >= 6 and all(r["call"] == "up" for r, _ in ieg_promoters)


def test_permuted_labels_give_far_fewer_calls(chain: dict) -> None:
    true_calls = (
        sum(r["call"] != "ns" for r in read_tsv(chain["work"] / "dar.tsv"))
        if (chain["work"] / "dar.tsv").exists()
        else None
    )
    if true_calls is None:
        true_calls = sum(
            r["call"] != "ns"
            for r in run_deseq2(chain, chain["sheet"], chain["work"] / "dar.tsv")
        )
    counts = []
    for g1 in itertools.combinations(range(6), 3):
        if 0 not in g1 or list(g1) == [0, 1, 2]:
            continue
        labels = ["control" if i in g1 else "stimulated" for i in range(6)]
        sheet = chain["work"] / f"perm_{''.join(map(str, g1))}.tsv"
        sheet.write_text(
            "sample_id\tcondition\tbatch\n"
            + "".join(f"{s}\t{c}\tb1\n" for s, c in zip(SAMPLES, labels, strict=True))
        )
        rows = run_deseq2(chain, sheet, chain["work"] / "perm_dar.tsv")
        counts.append(sum(r["call"] != "ns" for r in rows))
    assert len(counts) == 9
    assert statistics.median(counts) <= 0.05 * true_calls
