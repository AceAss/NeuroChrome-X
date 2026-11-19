# ADR-003: Age and activation are compared by overlap, not interaction

- **Status:** Accepted (PRD v1.0, 2026-11-05)
- **PRD reference:** §3 D3; §8 item 3
- **Deciders:** project owner; drafted by Claude

## Context

No public hippocampal dataset found so far has both age and neuronal activation as factors.

## Decision

- RQ3 is answered by concordance and overlap analyses with matched-background permutation tests: at least 10,000 permutations, matched on accessibility, GC content and distance to the TSS.
- An interaction model runs only if a dataset's design supports it.
- Outputs say overlap and concordance, never attenuation, unless an interaction test supports it.

## Consequences

- Claims stay at the level of association between independent experiments.
- The overlap definitions are frozen before aging data are examined (ADR-013).

## Alternatives considered

- Subtracting effects across studies to mimic an interaction: rejected, because protocols differ.

## Verification

- Each overlap result reports an effect size and a permutation p-value.
- The permutation null is tested on the toy dataset.

## Revisit when

- A dataset with both factors is found.
