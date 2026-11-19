# ADR-005: Methylation is reported as mC/hmC, CpG only in v0.1

- **Status:** Accepted (PRD v1.0, 2026-11-05)
- **PRD reference:** §3 D5; §11 coverage gate
- **Deciders:** project owner; drafted by Claude

## Context

Bisulfite sequencing reads 5mC and 5hmC together. Neurons also carry abundant non-CpG methylation. RRBS covers CpG-dense regions, so many distal enhancers get no coverage.

## Decision

- Bisulfite results are labelled mC/hmC, never 5mC.
- v0.1 analyses CpG context only; non-CpG methylation is a later stretch.
- A coverage gate (fraction of enhancers with at least 3 covered CpGs) runs per dataset before any enhancer-level claim. Its threshold is set in a later ADR after the first dataset.
- MeDIP data (A3) use an optional adapter and are never mixed with bisulfite effect sizes.
- "Likely DNMT target sites" stays undefined until OD3 is settled in P2.

## Consequences

- Enhancer-methylation claims may be dropped for RRBS datasets.
- WGBS data would remove that limit but are heavy to process (ADR-012).

## Alternatives considered

- Calling it 5mC: rejected as scientifically wrong for bisulfite data.
- oxBS or TAB-seq to separate 5mC from 5hmC: only if public data exist.

## Verification

- A coverage-gate table exists for every methylation dataset.
- A test checks that methylation outputs carry the mC/hmC label.

## Revisit when

- OD3 is settled.
- A suitable WGBS dataset becomes available.
