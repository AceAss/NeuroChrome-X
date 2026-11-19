"""The ATAC branch run through Snakemake on the synthetic toy dataset (PRD section 11).

Runs ``workflow/Snakefile`` with the toy configuration (results go to a temporary
directory), then checks the QC table, peak calls, DESeq2 results, annotation and motif
enrichment against the planted truth, and that a second run has nothing left to do.
Needs snakemake, bowtie2, samtools, bedtools, MACS3, deepTools, fastp and R with DESeq2
on PATH together (environment ncx-test). Runs for several minutes.

Set NCX_TEST_WORKDIR to a directory to keep the results between runs: Snakemake then
skips finished jobs, which makes reruns of this test (and debugging) fast.
"""

from __future__ import annotations

import os
import shutil
import statistics
import subprocess
from pathlib import Path

import pytest

pytestmark = pytest.mark.integration

ROOT = Path(__file__).resolve().parents[2]
TOOLS = (
    "snakemake",
    "bowtie2",
    "samtools",
    "bedtools",
    "macs3",
    "bamCoverage",
    "bigwigAverage",
    "fastp",
    "Rscript",
)
missing = [tool for tool in TOOLS if shutil.which(tool) is None]
if missing:
    pytest.skip(f"tools not on PATH: {', '.join(missing)}", allow_module_level=True)
if subprocess.run(["Rscript", "-e", "library(DESeq2)"], capture_output=True).returncode:
    pytest.skip("R package DESeq2 is not installed", allow_module_level=True)

CONFIGS = [
    str(ROOT / "config" / "config.yaml"),
    str(ROOT / "config" / "overrides" / "toy.yaml"),
]


def snakemake(extra: list[str], override: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [
            "snakemake",
            "-s",
            str(ROOT / "workflow" / "Snakefile"),
            "--directory",
            str(ROOT),
            "--cores",
            "2",
            "--rerun-incomplete",
            "--configfile",
            *CONFIGS,
            str(override),
            *extra,
        ],
        capture_output=True,
        text=True,
    )


def read_tsv(path: Path) -> list[dict[str, str]]:
    lines = path.read_text().splitlines()
    head = lines[0].split("\t")
    return [dict(zip(head, line.split("\t"), strict=True)) for line in lines[1:]]


def overlap(a: tuple[int, int], b: tuple[int, int]) -> bool:
    return a[0] < b[1] and b[0] < a[1]


@pytest.fixture(scope="module")
def run(tmp_path_factory: pytest.TempPathFactory) -> dict:
    kept = os.environ.get("NCX_TEST_WORKDIR")
    if kept:
        work = Path(kept)
        work.mkdir(parents=True, exist_ok=True)
    else:
        work = tmp_path_factory.mktemp("toy_snakemake")
    override = work / "override.yaml"
    override.write_text(
        f"paths:\n  results: {work}/results\n"
        f"reference:\n  files:\n    bowtie2_index: {work}/results/index/bowtie2/toy\n"
    )
    first = snakemake([], override)
    assert first.returncode == 0, first.stderr[-3000:]
    return {"results": work / "results" / "toy" / "atac", "override": override}


def truth_by_status() -> list[dict]:
    path = ROOT / "tests" / "data" / "toy" / "generated" / "truth_regions.tsv"
    return [
        {
            "iv": (int(r["start"]), int(r["end"])),
            "status": r["status"],
            "cls": r["class"],
            "lfc": float(r["true_log2fc"]),
        }
        for r in read_tsv(path)
    ]


def test_second_run_has_nothing_to_do(run: dict) -> None:
    again = snakemake(["-n"], run["override"])
    assert again.returncode == 0, again.stderr[-2000:]
    assert "Nothing to be done" in again.stdout + again.stderr


def test_qc_table(run: dict) -> None:
    rows = read_tsv(run["results"] / "qc" / "qc_table.tsv")
    assert len(rows) == 6
    for r in rows:
        assert float(r["mapping_rate"]) >= 0.99
        assert 0.08 <= float(r["dup_fraction"]) <= 0.12
        assert 0.10 <= float(r["chrM_fraction"]) <= 0.14
        assert 0.03 <= float(r["blacklist_fraction"]) <= 0.06
        assert float(r["frip"]) >= 0.5 and float(r["nrf"]) >= 0.85
        assert r["qc_flag"] == "pass"


def test_deseq2_results_match_truth(run: dict) -> None:
    rows = read_tsv(run["results"] / "dar" / "stim_vs_ctrl.tsv")
    truth = truth_by_status()
    by_status: dict[str, list[bool]] = {"gained": [], "lost": [], "unchanged": []}
    wrong_direction = 0
    medians = []
    for row in rows:
        hits = [
            t for t in truth if overlap((int(row["start"]), int(row["end"])), t["iv"])
        ]
        if len(hits) != 1:
            continue
        t = hits[0]
        called = row["call"] != "ns"
        by_status[t["status"]].append(called)
        if (
            called
            and t["status"] != "unchanged"
            and (row["call"] == "up") != (t["lfc"] > 0)
        ):
            wrong_direction += 1
        if t["status"] == "unchanged":
            medians.append(float(row["log2FC_mle"]))
    n_calls = sum(sum(v) for v in by_status.values())
    false_calls = sum(by_status["unchanged"])
    assert sum(by_status["gained"]) / len(by_status["gained"]) >= 0.90
    assert sum(by_status["lost"]) / len(by_status["lost"]) >= 0.80
    assert false_calls / len(by_status["unchanged"]) <= 0.05
    assert false_calls / n_calls <= 0.10
    assert wrong_direction == 0
    assert abs(statistics.median(medians)) <= 0.15


def test_annotation_and_links(run: dict) -> None:
    rows = read_tsv(run["results"] / "elements" / "regulatory_elements.tsv")
    classes = {
        c: sum(r["class"] == c for r in rows) for c in ("promoter", "enhancer", "other")
    }
    assert (
        classes["promoter"] == 60
        and classes["enhancer"] >= 330
        and classes["other"] == 0
    )
    assert all(r["linked_gene_id"] != "NA" for r in rows)
    assert sum(r["tier1"] == "true" for r in rows) >= 220


def test_motif_enrichment_finds_ap1_first(run: dict) -> None:
    up = read_tsv(run["results"] / "motifs" / "stim_vs_ctrl.up.tsv")
    assert up[0]["motif"] == "AP1"
    assert float(up[0]["padj"]) < 1e-8 and float(up[0]["odds_ratio"]) > 10
    assert {r["motif"] for r in up} == {"AP1", "CREB", "MEF2"}


def test_tracks_exist(run: dict) -> None:
    tracks = run["results"] / "tracks"
    assert len(list(tracks.glob("*.bw"))) == 6
    assert {p.name for p in (tracks / "conditions").glob("*.bw")} == {
        "control.mean.bw",
        "stimulated.mean.bw",
    }
