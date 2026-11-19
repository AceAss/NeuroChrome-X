# Shared setup: config, validated sample sheet, wildcard helpers (PRD sections 7 and 9).
import sys
from pathlib import Path

import yaml
from snakemake.exceptions import WorkflowError

ROOT = Path(workflow.basedir).parent
sys.path.insert(0, str(ROOT / "src"))  # lets the workflow run before `pip install -e .`
SCRIPTS = str(ROOT / "workflow" / "scripts")

if not config:  # no --configfile given: fall back to the defaults
    config.update(yaml.safe_load((ROOT / "config" / "config.yaml").read_text()))

from neurochromx.io.samplesheet import load_samplesheet  # noqa: E402

RES = config["paths"]["results"]
REF = config["reference"]["files"]
SHEET = load_samplesheet(config["paths"]["samples"])

ATAC = [r for r in SHEET.records if r.assay.value == "atac"]
RECORD = {r.sample_id: r for r in SHEET.records}
STUDIES = sorted({r.study_id for r in ATAC})
ATAC_BY_STUDY = {s: [r.sample_id for r in ATAC if r.study_id == s] for s in STUDIES}
PE_SAMPLES = [r.sample_id for r in ATAC if r.layout and r.layout.value == "PE"]
SE_SAMPLES = [r.sample_id for r in ATAC if r.layout and r.layout.value == "SE"]
CONDITIONS = {
    (s, c): [r.sample_id for r in ATAC if r.study_id == s and r.condition == c]
    for s in STUDIES
    for c in sorted({r.condition for r in ATAC if r.study_id == s})
}
CONTRASTS = {c["name"]: c for c in config.get("contrasts", [])}

for _r in ATAC:
    if _r.entry_point.value != "fastq":
        raise WorkflowError(
            f"{_r.sample_id}: entry_point {_r.entry_point.value} is not implemented yet "
            "for ATAC (only fastq); see ADR-012 and the README roadmap"
        )
    if not _r.path_or_uri:
        raise WorkflowError(
            f"{_r.sample_id}: no local FASTQ path. SRA download rules are not part of "
            "this workflow yet; download the reads and fill path_or_uri"
        )


def wildcard_alternatives(names):
    """Regex alternation for wildcard_constraints (never empty, so rules stay valid)."""
    return "|".join(sorted(set(names))) or "NO_SAMPLES"


def fastq_paths(sample):
    return [p.strip() for p in RECORD[sample].path_or_uri.split(";") if p.strip()]


def trimming_enabled():
    return config["trimming"].get("enabled", True)


def align_reads(wc):
    mates = 2 if wc.sample in PE_SAMPLES else 1
    if trimming_enabled():
        return [f"{RES}/{wc.study}/atac/trim/{wc.sample}_R{m}.fastq.gz" for m in range(1, mates + 1)]
    return fastq_paths(wc.sample)


def genome_size(chrom_sizes):
    keep = set(config["reference"]["main_chromosomes"])
    total = 0
    with open(chrom_sizes) as handle:
        for line in handle:
            name, length = line.split()[:2]
            if name in keep:
                total += int(length)
    return total


def macs_flags(sample):
    if sample in PE_SAMPLES:
        return "-f BAMPE"
    se = config["peaks"]["single_end_atac"]
    return f"-f BAM --nomodel --shift {se['shift']} --extsize {se['extsize']}"


def motif_args():
    engine = config["motifs"]
    spec = engine.get("python_engine", {})
    consensus = spec.get("consensus") or {}
    if consensus:
        return "--consensus " + " ".join(f"{k}={v}" for k, v in consensus.items())
    tfs = [t for group in engine["targeted"].values() for t in group] + list(engine.get("extended", []))
    return f"--meme {spec['meme']} --tfs " + " ".join(tfs)


def rna_inputs(study):
    """(counts file, sample ids) of the study's RNA samples given as counts, or None."""
    rows = [r for r in SHEET.records if r.study_id == study and r.assay.value == "rna"
            and r.entry_point.value == "counts"]
    if not rows:
        return None
    return rows[0].path_or_uri, [r.sample_id for r in rows]


def all_targets():
    targets = []
    for study in STUDIES:
        base = f"{RES}/{study}/atac"
        targets += [
            f"{base}/qc/qc_table.tsv",
            f"{base}/consensus/counts.tsv",
            f"{base}/elements/regulatory_elements.tsv",
        ]
        for c in CONTRASTS.values():
            if c["study_id"] == study and c["assay"] == "atac":
                targets.append(f"{base}/dar/{c['name']}.tsv")
                targets += [f"{base}/motifs/{c['name']}.{d}.tsv" for d in ("up", "down")]
        if config["tracks"].get("enabled", True):
            targets += [f"{base}/tracks/{s}.bw" for s in ATAC_BY_STUDY[study]]
            targets += [f"{base}/tracks/conditions/{c}.mean.bw" for (s, c) in CONDITIONS if s == study]
        if config["qc"].get("multiqc", False):
            targets.append(f"{base}/qc/multiqc_report.html")
    return targets
