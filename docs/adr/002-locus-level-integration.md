# ADR-002: Integrate at locus level, not by pooling samples

- **Status:** Accepted (PRD v1.0, 2026-11-05)
- **PRD reference:** §3 D2; §8
- **Deciders:** project owner; drafted by Claude

## Context

Public datasets are rarely paired across modalities, and studies differ in protocol, age, sex and cell population.

## Decision

- Differential analysis runs within each study and contrast. Results are stored as effect sizes with standard errors (`element_effects`, PRD §9).
- Integration happens at element or gene level: concordance, enrichment against matched backgrounds, and random-effects meta-analysis.
- Joint sample-level models such as MOFA2 are used only where several modalities come from the same animals. B1 may qualify; check its sample metadata.

## Consequences

- Raw counts are never pooled across studies.
- Sample-level features for BioAge-X come from per-sample scores (X1, X2), not from merged matrices.

## Alternatives considered

- A batch-corrected pooled matrix: rejected, because study and protocol are confounded with condition.
- MOFA+ across unmatched samples: rejected for the same reason.

## Verification

- Every integrated result lists its contributing study and contrast IDs.
- A unit test with shuffled study labels gives null concordance.

## Revisit when

- A truly paired multi-omic dataset is added.
