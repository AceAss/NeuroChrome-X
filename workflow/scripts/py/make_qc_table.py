#!/usr/bin/env python3
"""Per-sample QC statistics and the study QC table (PRD section 7, ADR-009).

Two subcommands:

``sample``  counts read pairs at each stage of one ATAC sample and writes one TSV row.
``table``   merges the per-sample rows, adds flags from the thresholds and writes the
            table.

Flags never exclude a sample (policy flag_only): a human signs off in the run report.
TSS enrichment here is a simple version: reads overlapping TSS +/-50 bp per base,
divided by reads overlapping the two flank windows (1900-2000 bp each side) per base. It
is comparable in spirit to the ENCODE score, not identical.
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "src"))

from neurochromx.annotate.promoters import read_genes

SAMPLE_COLUMNS = [
    "sample_id",
    "raw_pairs",
    "mapping_rate",
    "aligned_pairs",
    "q30_pairs",
    "dedup_pairs",
    "chrM_pairs",
    "kept_pairs",
    "final_pairs",
    "n_peaks",
    "peak_read_pairs",
    "nrf",
    "pbc1",
    "pbc2",
    "tss_center_density",
    "tss_flank_density",
]
TABLE_COLUMNS = [
    "sample_id",
    "raw_pairs",
    "mapping_rate",
    "q30_pairs",
    "mapq_removed_fraction",
    "dup_fraction",
    "chrM_fraction",
    "blacklist_fraction",
    "final_pairs",
    "n_peaks",
    "frip",
    "nrf",
    "pbc1",
    "pbc2",
    "tss_enrichment",
    "qc_flag",
]


def run(command: str) -> str:
    done = subprocess.run(
        command,
        shell=True,
        check=True,
        capture_output=True,
        text=True,
        executable="/bin/bash",
    )
    return done.stdout


def count_pairs(bam: str, layout: str, region: str = "", extra_flags: str = "") -> int:
    flags = f"-f 64 {extra_flags}" if layout == "PE" else f"-F 2308 {extra_flags}"
    return int(run(f"samtools view -c {flags} {bam} {region}").split()[0])


def bowtie2_summary(log: Path) -> tuple[int, float]:
    text = log.read_text()
    reads = int(re.match(r"\s*(\d+) reads", text).group(1))
    rate = float(re.search(r"([\d.]+)% overall alignment rate", text).group(1)) / 100
    return reads, rate


def complexity(q30_bam: str) -> tuple[float, float, float]:
    """NRF, PBC1, PBC2 on unique fragments of proper pairs (ENCODE definitions)."""
    awk_fragment = (
        'awk \'BEGIN{OFS="\\t"} '
        "{t=$9; if (t<0) t=-t; s=($9>0)?$4:$8; print $3,s,s+t}'"
    )
    awk_counts = "awk '{c[$1]++; n+=$1; d++} END{print n+0, d+0, c[1]+0, c[2]+0}'"
    out = run(
        f"samtools view -f 66 {q30_bam} | {awk_fragment} | sort -S 1G | uniq -c "
        f"| {awk_counts}"
    ).split()
    total, distinct, one, two = (int(x) for x in out)
    nan = float("nan")
    return (
        distinct / total if total else nan,
        one / distinct if distinct else nan,
        one / two if two else nan,
    )


def tss_densities(
    final_bam: str, gtf: Path, gene_types: list[str]
) -> tuple[float, float]:
    genes = read_genes(gtf, gene_types=gene_types or None)
    windows = []
    for gene in genes:
        tss = gene.tss
        for kind, lo, hi in (
            ("center", tss - 50, tss + 50),
            ("flank", tss - 2000, tss - 1900),
            ("flank", tss + 1900, tss + 2000),
        ):
            if lo >= 0:
                windows.append((gene.chrom, lo, hi, kind))
    with tempfile.TemporaryDirectory() as tmp:
        bed = Path(tmp) / "tss_windows.bed"
        bed.write_text("".join(f"{c}\t{s}\t{e}\t{k}\n" for c, s, e, k in windows))
        out = run(f"bedtools multicov -p -bams {final_bam} -bed {bed}")
    totals = {"center": [0, 0], "flank": [0, 0]}
    for line in out.splitlines():
        _, start, end, kind, count = line.split("\t")
        totals[kind][0] += int(count)
        totals[kind][1] += int(end) - int(start)
    density = {k: (v[0] / v[1] if v[1] else float("nan")) for k, v in totals.items()}
    return density["center"], density["flank"]


def cmd_sample(a: argparse.Namespace) -> int:
    reads, rate = bowtie2_summary(a.bowtie2_log)
    pe = a.layout == "PE"
    contigs = " ".join(a.keep_contigs)
    aligned = count_pairs(a.raw_bam, a.layout, extra_flags="-f 2" if pe else "")
    n_peaks = sum(1 for line in a.peaks.read_text().splitlines() if line.strip())
    peak_reads = int(
        run(
            f"bedtools intersect -u -abam {a.final_bam} -b {a.peaks} | "
            f"samtools view -c {'-f 64' if pe else '-F 2308'} -"
        ).split()[0]
    )
    nrf, pbc1, pbc2 = complexity(a.q30_bam) if pe else (float("nan"),) * 3
    center, flank = tss_densities(a.final_bam, a.gtf, a.gene_type)
    row = [
        a.sample,
        reads,
        rate,
        aligned,
        count_pairs(a.q30_bam, a.layout),
        count_pairs(a.dedup_bam, a.layout),
        count_pairs(a.dedup_bam, a.layout, a.mito),
        count_pairs(a.dedup_bam, a.layout, contigs),
        count_pairs(a.final_bam, a.layout),
        n_peaks,
        peak_reads,
        nrf,
        pbc1,
        pbc2,
        center,
        flank,
    ]
    a.out.write_text(
        "\t".join(SAMPLE_COLUMNS) + "\n" + "\t".join(str(x) for x in row) + "\n"
    )
    return 0


def ratio(num: float, den: float) -> float:
    return num / den if den else float("nan")


def cmd_table(a: argparse.Namespace) -> int:
    lines = ["\t".join(TABLE_COLUMNS)]
    for path in a.stats:
        head, values = path.read_text().splitlines()[:2]
        s = dict(zip(head.split("\t"), values.split("\t"), strict=True))
        f = {k: float(v) for k, v in s.items() if k != "sample_id"}
        frip = ratio(f["peak_read_pairs"], f["final_pairs"])
        tss = ratio(f["tss_center_density"], f["tss_flank_density"])
        flags = []
        if frip < a.frip_min:
            flags.append("low_frip")
        if tss < a.tss_min:
            flags.append("low_tss")
        if f["nrf"] < a.nrf_min:
            flags.append("low_nrf")
        if f["pbc1"] < a.pbc1_min:
            flags.append("low_pbc1")
        row = [
            s["sample_id"],
            int(f["raw_pairs"]),
            f"{f['mapping_rate']:.4f}",
            int(f["q30_pairs"]),
            f"{1 - ratio(f['q30_pairs'], f['aligned_pairs']):.4f}",
            f"{1 - ratio(f['dedup_pairs'], f['q30_pairs']):.4f}",
            f"{ratio(f['chrM_pairs'], f['dedup_pairs']):.4f}",
            f"{1 - ratio(f['final_pairs'], f['kept_pairs']):.4f}",
            int(f["final_pairs"]),
            int(f["n_peaks"]),
            f"{frip:.4f}",
            f"{f['nrf']:.4f}",
            f"{f['pbc1']:.4f}",
            f"{f['pbc2']:.4f}",
            f"{tss:.3f}",
            ",".join(flags) or "pass",
        ]
        lines.append("\t".join(str(x) for x in row))
    a.out.write_text("\n".join(lines) + "\n")
    flagged = sum(not line.endswith("\tpass") for line in lines[1:])
    print(f"QC table for {len(a.stats)} samples; flagged: {flagged}")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = parser.add_subparsers(dest="command", required=True)
    s = sub.add_parser("sample", help="count pairs at each stage for one sample")
    s.add_argument("--sample", required=True)
    s.add_argument("--layout", choices=["PE", "SE"], default="PE")
    for name in ("raw-bam", "q30-bam", "dedup-bam", "final-bam"):
        s.add_argument(f"--{name}", required=True)
    s.add_argument("--bowtie2-log", type=Path, required=True)
    s.add_argument("--peaks", type=Path, required=True)
    s.add_argument("--gtf", type=Path, required=True)
    s.add_argument("--gene-type", nargs="*", default=["protein_coding"])
    s.add_argument("--keep-contigs", nargs="+", required=True)
    s.add_argument("--mito", default="chrM")
    s.add_argument("--out", type=Path, required=True)
    s.set_defaults(func=cmd_sample)
    t = sub.add_parser("table", help="merge per-sample rows and add flags")
    t.add_argument("--stats", nargs="+", type=Path, required=True)
    t.add_argument("--frip-min", type=float, default=0.15)
    t.add_argument("--tss-min", type=float, default=5.0)
    t.add_argument("--nrf-min", type=float, default=0.8)
    t.add_argument("--pbc1-min", type=float, default=0.8)
    t.add_argument("--out", type=Path, required=True)
    t.set_defaults(func=cmd_table)
    args = parser.parse_args(argv)
    for attr in ("raw_bam", "q30_bam", "dedup_bam", "final_bam"):
        if hasattr(args, attr):
            setattr(args, attr, str(getattr(args, attr)))
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
