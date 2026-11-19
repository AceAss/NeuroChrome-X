# ADR-011: Python package plus R scripts, no rpy2

- **Status:** Accepted (PRD v1.0, 2026-11-05)
- **PRD reference:** §3 D11
- **Deciders:** project owner; drafted by Claude

## Context

The best statistical tools for count data (DESeq2, csaw, DSS) are R packages; orchestration, annotation, integration and export are easier to test in Python.

## Decision

- Python package `neurochromx` for orchestration, annotation, integration and export.
- R/Bioconductor scripts for DESeq2, csaw and DSS, called by Snakemake with command-line arguments.
- Data pass between languages as TSV or Parquet files. No rpy2 in v0.1.

## Consequences

- Two toolchains to maintain, but a clean file-based interface.

## Alternatives considered

- pydeseq2 and Python-only: less validated for this use.
- rpy2: tighter coupling and harder environment management.

## Verification

- Every R script has a test on toy inputs.
- Input and output files follow PRD §9.

## Revisit when

- A Python implementation matches the R results on the toy and real data.
