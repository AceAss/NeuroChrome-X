# ADR-008: Normalization of accessibility counts

- **Status:** Accepted (PRD v1.0, 2026-11-05)
- **PRD reference:** §3 D8; §7
- **Deciders:** project owner; drafted by Claude

## Context

Peak-based normalization assumes most peaks do not change, which can hide global shifts such as strong activation or aging.

## Decision

- Primary: DESeq2 on consensus-peak counts, design `~ batch + condition`, apeglm shrinkage, FDR below 0.05 and absolute log2FC at least 0.585.
- Sensitivity: csaw-style normalization on large background bins (10 kb).
- Report both. Flag any dataset where the DAR counts differ by more than 2-fold or where effect directions flip.

## Consequences

- Two runs per differential analysis, and one extra R environment package (csaw).

## Alternatives considered

- TMM only.
- Spike-in normalization: not available in public data.

## Verification

- Both normalizations run on the toy data.
- The sensitivity table exists for every dataset.
- Evidence from the toy dataset (Claude's sandbox, `tests/data/toy/README.md`): library-size (CPM) normalisation with a Welch t-test called 38 of 337 unchanged peaks, while DESeq2's median-of-ratios called 3, because strong gains take a larger share of the reads.

## Revisit when

- A dataset with spike-ins is added.
