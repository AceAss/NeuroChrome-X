# ADR-009: Unique mapping is primary; repeat accessibility is analysed separately

- **Status:** Accepted (PRD v1.0, 2026-11-05)
- **PRD reference:** §3 D9
- **Deciders:** project owner; drafted by Claude

## Context

The aging dataset B1 reports higher accessibility with age concentrated in repetitive elements, which a MAPQ filter of 30 largely removes.

## Decision

- The primary analysis keeps MAPQ at least 30 (unique mapping).
- A separate multi-mapper and repeat-class analysis is a stretch. Whenever it is run, both results are reported.
- QC records the fraction of reads removed by the MAPQ filter per sample.

## Consequences

- Primary results are blind to most repeat-derived signal; say so in Methods.

## Alternatives considered

- Keeping multi-mappers in the primary analysis: rejected, because peak calls and counts become unreliable.

## Verification

- The MAPQ-removed fraction appears in `qc_table.tsv`.

## Revisit when

- Repeat-class analysis is scheduled.
