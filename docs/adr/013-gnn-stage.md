# ADR-013: The GNN stage is exploratory, with frozen definitions

- **Status:** Accepted (PRD v1.0, 2026-11-05)
- **PRD reference:** §3 D13; §8 items 5-6; §15
- **Deciders:** project owner; drafted by Claude

## Context

Most aging datasets have two age levels and a handful of replicates, so flexible models overfit easily and age regression is poorly identified.

## Decision

- Use classification or contrast-level targets, leave-one-study-out cross-validation, and ridge/ElasticNet baselines.
- A GNN claim requires beating the baselines; otherwise it is reported as a negative result.
- Element and feature definitions are frozen on activity data, recorded in an ADR with a file hash, before any aging data are examined.

## Consequences

- Some planned analyses may end as negative results, and that is acceptable.

## Alternatives considered

- Age regression across studies: rejected, with so few age levels.
- Tuning the feature definitions on aging data: rejected as circular.

## Verification

- The baseline gate in PRD §11.
- The run manifest includes the hash of the frozen definition.

## Revisit when

- A dataset with a dense age series is added.
