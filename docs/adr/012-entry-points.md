# ADR-012: Every stage accepts processed inputs

- **Status:** Accepted (PRD v1.0, 2026-11-05)
- **PRD reference:** §3 D12; §9
- **Deciders:** project owner; drafted by Claude

## Context

Re-aligning raw reads (especially WGBS) is unrealistic on a laptop, and many studies publish usable processed files.

## Decision

- `samples.tsv` has an `entry_point` column: fastq, bam, peaks, counts or cov.
- Each workflow stage declares which entry points it accepts, and processed inputs are normalised into canonical formats (narrowPeak, counts TSV, Bismark coverage).
- The validator enforces the allowed assay and entry-point combinations and requires `genome_native` for every processed input.

## Consequences

- More code paths to test, but projects can start from processed data.
- Processed inputs inherit the authors' pipeline choices, which must be recorded in the registry.

## Alternatives considered

- FASTQ only: rejected, as it would block methylation work on modest hardware.

## Verification

- The sample-sheet tests (CL-003) cover the combinations.
- The toy dataset uses fastq (ATAC) and counts (RNA).

## Revisit when

- Cluster access removes the need for processed inputs.
