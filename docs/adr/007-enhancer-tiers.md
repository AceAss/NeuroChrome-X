# ADR-007: Enhancer tiers and gene linking

- **Status:** Accepted (PRD v1.0, 2026-11-05)
- **PRD reference:** §3 D7; §8 item 1
- **Deciders:** project owner; drafted by Claude

## Context

There is no gold-standard list of mouse hippocampal enhancers, and enhancer-to-gene links are models rather than facts.

## Decision

- T1: overlap with ENCODE mm10 enhancer-like cCREs (pELS, dELS).
- T2: same-tissue chromatin-state evidence, H3K27ac and H3K4me1 without H3K4me3.
- T3: a distal ATAC peak, more than 2 kb from a TSS.
- Gene link in v0.1: the nearest expressed TSS within 100 kb. TAD- or ABC-based linking is a stretch; the Hi-C data in A2 can validate links.
- Every result table carries an enhancer-tier column, and results are reported per tier.

## Consequences

- Three definitions multiply result tables, but conclusions that hold in all tiers are stronger.

## Alternatives considered

- A single definition: rejected, because there is no reference standard.
- Nearest gene with no distance limit: rejected.

## Verification

- Tier-assignment unit tests on the toy truth table.
- Tier overlap is reported for every dataset.

## Revisit when

- Hi-C or ABC-based links are added (stretch).
