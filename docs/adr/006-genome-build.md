# ADR-006: Genome build mm10 with GENCODE vM25 and UCSC chromosome names

- **Status:** Accepted (PRD v1.0, 2026-11-05)
- **PRD reference:** §3 D6
- **Deciders:** project owner; drafted by Claude

## Context

Most legacy data, the ENCODE mm10 blacklist and the cCRE registry are mm10. Mixing chromosome-naming styles is a common source of silent errors.

## Decision

- All results use mm10 (GRCm38) and GENCODE vM25.
- Re-align from FASTQ whenever possible, which avoids liftOver.
- Processed inputs on another build (mm9, mm39) are flagged by the sample-sheet validator and are either lifted over with a documented chain file or excluded.
- Chromosome names use the UCSC style (chr1, chrX, chrM) everywhere. TC-006 checks consistency across FASTA, GTF, blacklist and cCREs.

## Consequences

- mm39 resources are not used in v0.x.
- Any model or resource available only on another build needs a documented coordinate conversion.

## Alternatives considered

- mm39: newer, but fewer matching resources.
- Ensembl-style names without chr: rejected, to avoid mixed naming.

## Verification

- TC-006 naming-consistency report.
- The validator warns on processed inputs that are not on mm10.

## Revisit when

- A key resource is only available on mm39.
