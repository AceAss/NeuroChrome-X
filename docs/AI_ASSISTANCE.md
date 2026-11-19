# AI assistance log

Required by the explain-it gate (PRD section 4). One row per file or module. "What I verified" is for you to fill in: running the tests yourself and being able to explain the module count as verification; "Claude says it works" does not.

| Date | File or module | What AI did | What I verified |
|---|---|---|---|
| 2026-11-09 | `docs/PRD.md`, `docs/adr/001-014` | Claude drafted the PRD and the 14 decision records | |
| 2026-11-09 | `config/config.yaml`, `config/overrides/toy.yaml` | Claude wrote the configuration and checked it against PRD section 7 by script | |
| 2026-11-09 | `config/gene_sets/*_v1.tsv` | Claude curated the lists from its knowledge of the literature; not yet cross-checked against the published tables | |
| 2026-11-09 | `src/neurochromx/io/samplesheet.py`, `tests/unit/test_samplesheet.py` | Claude wrote the validator and 26 tests; tests run in Claude's sandbox | |
| 2026-11-09 | `tests/data/toy/make_toy.py`, `tests/unit/test_toy_generator.py` | Claude wrote the synthetic dataset generator and 12 tests | |
| 2026-11-09 | `src/neurochromx/intervals.py`, `annotate/`, `motifs/`, `tests/unit/test_annotate.py`, `tests/unit/test_motifs.py` | Claude wrote the annotation, linking, motif and background-matching modules and 15 tests | |
| 2026-11-09 | `workflow/scripts/py/consensus_peaks.py`, `count_peaks.py`, `workflow/scripts/R/deseq2_dar.R`, `tests/integration/test_toy_atac_chain.py` | Claude wrote the ATAC-chain scripts and ran them with real tools on the toy data in its sandbox | |
| 2026-11-09 | `docs/task-cards/TC-002` to `TC-007` | Claude wrote the task cards for Antigravity | |
| 2026-11-15 | `workflow/Snakefile`, `workflow/rules/*.smk`, `workflow/envs/*.yaml` | Claude wrote the ATAC-branch workflow and ran all 68 jobs on the toy data in its sandbox (without conda); the real environments are untested | |
| 2026-11-15 | `workflow/scripts/py/{annotate_elements,expressed_genes,motif_enrichment,make_qc_table,write_group_sheet}.py`, `src/neurochromx/{sequences.py,motifs/meme.py}` | Claude wrote the workflow helper scripts and their unit tests | |
| 2026-11-15 | `src/neurochromx/integrate/*.py`, `tests/unit/test_integrate.py` | Claude wrote the effect-table, master-catalogue, concordance, matched-overlap and meta-analysis code with tests on synthetic data | |
| 2026-11-15 | `tests/integration/test_toy_snakemake.py` | Claude wrote the end-to-end workflow test | |
