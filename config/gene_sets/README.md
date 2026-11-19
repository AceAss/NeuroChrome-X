# Gene sets (CL-004)

| File | Genes | Purpose |
|---|---|---|
| `iegs_v1.tsv` | 49 (25 tier 1, 24 tier 2) | Immediate-early and activity-regulated genes: positive controls, IGV bookmarks, locus-level tests |
| `plasticity_v1.tsv` | 37 (22 tier 1, 15 tier 2) | Synaptic-plasticity, memory and epigenetic-regulator loci |

These replace the `*_seed.tsv` lists created by TC-001. Delete the seeds with `git rm` once TC-006 reports that the v1 lists map cleanly. All 49 symbols from the PRD Appendix D seeds are kept; `Jund` and `Dusp5` moved from the IEG seeds to tier 2 because their induction is variable.

## Columns

| Column | Meaning |
|---|---|
| `symbol` | Mouse gene symbol as used by MGI and GENCODE vM25 |
| `aliases` | Alternative names, `\|`-separated, used only to map symbols across annotation releases |
| `category` | Functional class (tf, effector, regulator, kinase, phosphatase, nmdar, ampar, scaffold, epigenetic_regulator, ...) |
| `tier` | 1 = canonical, use for positive controls and headline tests; 2 = reported or plausible, use for sensitivity analyses |
| `source` | Where the entry comes from and what still needs checking |
| `curated_on` | Date of this version |
| `status` | `curated_v1` (tier 1), `candidate_v1` (tier 2), later `frozen_vN` |
| `notes` | Free text |

## Rules

1. Default filters: `tier == 1` for positive controls (PRD §11) and `tier <= 2` for sensitivity analyses.
2. Mapping: match `symbol`, then any alias, to the GENCODE vM25 gene name. Never map by position. TC-006 writes `resources/annotation/geneset_mapping.tsv` and lists every unmapped symbol.
3. Honest status: these lists come from the curator's knowledge of the canonical literature. They have **not** yet been checked against the Tyssowski 2018 and Hrvatin 2018 supplementary tables. That cross-check is required before freezing.
4. Freezing: after the cross-check, set `status` to `frozen_v1`, update `curated_on`, record the SHA-256 of each file in an ADR, and do not edit the files afterwards. The freeze must happen before any aging analysis (PRD D13, ADR-013).
5. `Bdnf` is multi-promoter: track promoters and exons separately when analysing accessibility.
6. `Homer1`, `Arc` and `Bdnf` are kept only in the IEG file; the plasticity file has no overlap with it.
