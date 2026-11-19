# ADR-001: Tool stack beyond the original brief

- **Status:** Accepted (PRD v1.0, 2026-11-05)
- **PRD reference:** §3 D1; §7
- **Deciders:** project owner; drafted by Claude

## Context

The original brief lists FastQC, Bowtie2, SAMtools, bedtools and deepTools. Bowtie2 is neither splice-aware nor bisulfite-aware, so it cannot align RNA-seq or bisulfite reads. The brief also leaves out trimming, duplicate handling, peak calling and report aggregation.

## Decision

- ATAC and ChIP: fastp (adapter trimming only), then Bowtie2, samtools filtering, Picard MarkDuplicates and MACS3.
- RNA-seq: STAR with featureCounts when the machine has about 32 GB of RAM or more; otherwise Salmon, quantification only. The machine audit (TC-002) decides and records it in `compute.rna_path`.
- Bisulfite: Bismark with its Bowtie2 back end.
- QC and tracks: FastQC, MultiQC, a custom `qc_table.tsv`, and deepTools bamCoverage. samtools and bedtools stay as general glue.

## Consequences

- More tools mean more conda environments (ADR-010).
- STAR needs roughly 30 GB of RAM for a mouse index (verify on this machine).
- Bismark genome preparation is large, so bisulfite from FASTQ is gated by disk and RAM (ADR-012).

## Alternatives considered

- HISAT2 or kallisto for RNA-seq: lower memory; the fallback if neither STAR nor Salmon fits.
- bwa-mem2 for DNA reads: no benefit over Bowtie2 here.
- MACS2: superseded by MACS3.

## Verification

- `docs/tool_versions.md` lists every tool with its version (TC-002).
- The toy ATAC run passes end to end.
- `docs/machine_audit.md` records the chosen RNA path.

## Revisit when

- The machine audit reports less than about 32 GB of RAM.
- A chosen dataset's protocol needs a different aligner.
