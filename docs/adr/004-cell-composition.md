# ADR-004: Cell composition is a first-class confound

- **Status:** Accepted (PRD v1.0, 2026-11-05)
- **PRD reference:** §3 D4; §12
- **Deciders:** project owner; drafted by Claude

## Context

Aging shifts cell-type proportions in bulk hippocampus, which can mimic a change in accessibility or expression within cells.

## Decision

- Prefer sorted-nuclei datasets.
- For bulk tissue, add composition proxies (marker-gene accessibility or expression for neurons, astrocytes, microglia and oligodendrocytes) as covariates, or run them as a sensitivity analysis.
- Label every bulk-derived result as bulk, in outputs and in the README.
- Module X6b adds pseudobulk single-nucleus data as a check.

## Consequences

- More covariates and some loss of power, but defensible results.

## Alternatives considered

- Ignoring composition: rejected.
- Full deconvolution with a reference atlas: deferred until a matching reference is chosen.

## Verification

- Each differential result has a paired sensitivity result.
- Composition proxies are stored in the run manifest.

## Revisit when

- A cell-type-resolved aging dataset is processed (X6b).
