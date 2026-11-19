# ADR-014: Reproducibility baseline

- **Status:** Accepted (PRD v1.0, 2026-11-05)
- **PRD reference:** §3 D14
- **Deciders:** project owner; drafted by Claude

## Context

Results must be re-creatable by someone else, and every number must be traceable to inputs, code and settings.

## Decision

- Pinned conda environments with committed lockfiles.
- Config-driven runs: `config/config.yaml` plus per-study overrides.
- Provenance for every download (accession or URL, date, checksum), reusing the BioAge-X layer where TC-003 shows it fits.
- A fixed seed (`project.seed`) for every stochastic step.
- A run manifest next to every result directory: git SHA, config hash, tool versions and input checksums.

## Consequences

- Small overhead per run, in exchange for traceability.

## Alternatives considered

- Notebook-only analysis: rejected, because it is hard to reproduce.

## Verification

- The determinism gate in PRD §11: a re-run reproduces tables and peak sets.
- CI runs the toy end-to-end test on every push.

## Revisit when

- A move to containers changes how tool versions are pinned.
