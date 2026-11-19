# ADR-010: Snakemake with per-rule conda environments

- **Status:** Accepted (PRD v1.0, 2026-11-05)
- **PRD reference:** §3 D10; §10
- **Deciders:** project owner; drafted by Claude

## Context

The project needs a reproducible, inspectable workflow that a single developer can run locally and later on a cluster.

## Decision

- Snakemake version 8 or later with per-rule conda environments. Containers (Apptainer) come later.
- nf-core/atacseq is an external benchmark on one dataset (peak overlap), not the engine.

## Consequences

- Environments are pinned with lockfiles (ADR-014).
- `snakemake -n` is part of CI.

## Alternatives considered

- Nextflow with nf-core pipelines: robust, but a heavier stack and less control for a portfolio project.
- Plain scripts or Make: not reproducible enough.

## Verification

- `snakemake -n` passes on every push.
- Lockfiles are committed for every environment.

## Revisit when

- A cluster or cloud run needs a different executor.
