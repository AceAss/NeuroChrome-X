# Architecture decision records

Accepted with PRD v1.0 on 2026-11-05. ADR-001 to ADR-014 record the design decisions D1 to
D14 of `docs/PRD.md` section 3. Agents must not change an accepted ADR; propose a new
one instead (see `000-template.md`).

| ADR | Title | PRD | Status |
|---|---|---|---|
| [ADR-001](001-tool-stack.md) | Tool stack beyond the original brief | D1 | Accepted |
| [ADR-002](002-locus-level-integration.md) | Integrate at locus level, not by pooling samples | D2 | Accepted |
| [ADR-003](003-age-activity-overlap.md) | Age and activation are compared by overlap, not interaction | D3 | Accepted |
| [ADR-004](004-cell-composition.md) | Cell composition is a first-class confound | D4 | Accepted |
| [ADR-005](005-methylation-readout.md) | Methylation is reported as mC/hmC, CpG only in v0.1 | D5 | Accepted |
| [ADR-006](006-genome-build.md) | Genome build mm10 with GENCODE vM25 and UCSC chromosome names | D6 | Accepted |
| [ADR-007](007-enhancer-tiers.md) | Enhancer tiers and gene linking | D7 | Accepted |
| [ADR-008](008-count-normalization.md) | Normalization of accessibility counts | D8 | Accepted |
| [ADR-009](009-repeats.md) | Unique mapping is primary; repeat accessibility is analysed separately | D9 | Accepted |
| [ADR-010](010-workflow-engine.md) | Snakemake with per-rule conda environments | D10 | Accepted |
| [ADR-011](011-languages.md) | Python package plus R scripts, no rpy2 | D11 | Accepted |
| [ADR-012](012-entry-points.md) | Every stage accepts processed inputs | D12 | Accepted |
| [ADR-013](013-gnn-stage.md) | The GNN stage is exploratory, with frozen definitions | D13 | Accepted |
| [ADR-014](014-reproducibility.md) | Reproducibility baseline | D14 | Accepted |

Planned ADRs (not yet written): enhancer-gene link thresholds, methylation coverage-gate
threshold, nf-core peak-overlap threshold, definition of "likely DNMT target sites" (OD3),
BioAge-X species handling (OD6), freeze of the plasticity-index definition (X2).
