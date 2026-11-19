# Explain-it notes (PRD section 4 gate)

A module is not done until you can explain it in your own words. Read each section, close the file, then answer the questions aloud or in writing. The answer hints are short on purpose: use your own words. Log what you checked in `docs/AI_ASSISTANCE.md`.

## 1. Sample-sheet validator (`io/samplesheet.py`)

1. `samples.tsv` is the only input the workflow trusts, so every row is checked before anything runs.
2. Each row is validated twice: controlled vocabularies (assay, tissue, entry point, genome build) and cross-field rules (chip needs a target; paired-end FASTQ needs two paths; processed inputs must name their genome build).
3. The whole sheet is checked too: unique sample IDs, unique replicate slots, and soft warnings (a group with one sample, a stimulated sample without a time point, processed data on mm9).
4. Errors are collected and reported together, with row numbers, so you fix a sheet in one pass.
5. Missing values are a small fixed set of tokens (empty, NA, none, ...), which is why a control stimulus is written `home_cage`, not `none`.

Questions:
- Why are small groups warnings and not errors? (A run with n = 1 is possible but statistically weak; the researcher must decide, not the code.)
- Why must processed inputs state `genome_native`? (Coordinates mean nothing without the build; mm9 peaks on an mm10 pipeline silently give wrong overlaps.)
- Why pydantic for rows but plain Python for the whole-sheet checks? (Rows are independent and typed; sheet-level rules need all rows at once.)
- What bug would a duplicated replicate slot cause later? (The same biological group counted twice, inflating confidence.)

## 2. Toy dataset generator (`tests/data/toy/make_toy.py`)

1. It builds a tiny fake genome, genes, accessible regions and reads in which the truth is known: which regions gain or lose accessibility, and which carry which motif.
2. Because the truth is known, you can measure the pipeline: did it find the gained regions, how many false calls, did AP-1 come out on top.
3. Motif counts are planted as exact numbers (for example 38 of 48 gained regions), not probabilities, so tests do not fail by bad luck.
4. It adds realistic nuisances: duplicates (10%), chrM reads (12%), blacklist pile-ups (4%), uniform background (18%), nucleosome-like fragment sizes.
5. All randomness comes from one seed, and the MANIFEST holds hashes, so a rerun is provably identical.

Questions:
- What can the toy not tell you? (Anything about real biology: cell-type mixtures, batch effects, real noise structure, enhancer biology.)
- Why synthetic instead of a small real dataset? (No download, exact truth, fast CI, no licence or privacy questions.)
- Why do the control and stimulated libraries differ in composition, and why does that matter? (The gained regions absorb reads, which biases library-size normalisation; see section 3.)

## 3. ATAC chain and the control gates (`workflow/scripts/`, `tests/integration/`)

1. Reads are aligned with Bowtie2 (up to 2 kb fragments), then filtered to MAPQ 30 and proper pairs (`-f 2 -F 1804` also drops unmapped, secondary, QC-fail and duplicate-flagged reads).
2. Duplicates are removed, chrM and blacklist reads are dropped, and MACS3 calls peaks per sample.
3. Consensus peaks are the merged peaks seen in at least 2 replicates of one condition; reads are counted in them and DESeq2 tests the contrast.
4. Gates compare the output with the planted truth: recovery of gained regions, false calls among unchanged regions, wrong-direction calls, the median log2FC of unchanged regions, and the number of calls under permuted labels.
5. In the sandbox run, library-size (CPM) normalisation made 38 unchanged peaks look different, while DESeq2's median-of-ratios called 3. That gap is composition bias: strong gains take a larger share of the reads.

Questions:
- Why remove chrM and blacklist regions? (chrM is a huge open-chromatin-like artefact in ATAC; blacklisted regions pile up reads in every assay and create false peaks.)
- Why require 2 replicates for a consensus peak? (Single-replicate peaks are mostly noise; replicate support is a cheap reproducibility filter.)
- Why is a permuted-label gate used, and why does a 3 v 3 design give only 9 other splits? (Permutation should remove the real signal; with 6 samples there are 10 distinct 3 v 3 splits, one is the truth, so 9 remain, and they still share part of the truth.)
- Why not use MAPQ 30 for repeat-derived signal? (Multi-mappers get low MAPQ and are removed; ADR-009 treats repeats as a separate analysis.)

## 4. Annotation and gene linking (`annotate/`)

1. A GTF gives genes and TSSs (TSS is the first transcribed base, so it depends on strand). Promoters are windows around the TSS (default 1000 up, 500 down).
2. Elements that overlap a promoter are promoters. The rest are enhancers if they pass any of three tiers: cCRE overlap (T1), histone evidence H3K27ac + H3K4me1 without H3K4me3 (T2), or a distal ATAC peak (T3).
3. The tiers are not exclusive, and every result is reported per tier, because no single enhancer definition is a gold standard.
4. A link assigns an element to the nearest expressed TSS within 100 kb, measured from the element centre.
5. A link is a model, not a measured contact: it is not evidence that the gene is the target.

Questions:
- Why 0-based half-open coordinates? (BED convention; avoids off-by-one errors when intervals touch.)
- Why only expressed genes for linking? (A silent gene is an unlikely target; it removes many false links.)
- How could Hi-C or ABC change the result? (They use 3D contact or activity-by-contact to link elements to distant genes, beyond nearest-gene.)

## 5. Motif scanning, background matching and enrichment (`motifs/`)

1. A motif is a position weight matrix of log-odds scores; a sequence scores its best window on either strand.
2. The score is rescaled to 0-1, and a hit means at least 0.9 of the way from the worst to the best possible match (an exact 7-mer match scores 1.0; one mismatch about 0.86).
3. The foreground (for example gained peaks) is compared with a background of non-differential peaks matched on GC and accessibility, so that generic promoter or open-chromatin motifs do not look enriched.
4. Presence versus absence is compared with a one-sided Fisher test, then corrected across motifs (Benjamini-Hochberg).
5. This code is for the toy gates and targeted tests. On real data, HOMER and FIMO/AME with JASPAR matrices are the tools (PRD section 7).

Questions:
- Why match the background on GC? (GC-rich sequences contain more CpG- and GC-box-like motifs by chance.)
- Why is a random split of unchanged peaks a useful negative control? (It must give no enrichment; if it does, the test or the background is biased.)
- What does an odds ratio of 20 with n = 48 tell you, and what does it not tell you? (A strong association in these sequences; not that the TF binds, nor that it causes the change.)

## 6. The Snakemake workflow (`workflow/`)

1. A Snakemake rule says how to make one kind of file from others; Snakemake builds the dependency graph (DAG) from the files you ask for and runs only what is missing or out of date.
2. Every parameter comes from `config/config.yaml`; a second `--configfile` (the toy override) changes keys one by one. Nothing is hard-coded, so the same rules run on the toy data and on real data.
3. Wildcards like `{sample}` are patterns, so they need constraints. Without them, the name `x.q30` could match `{sample}` of the BAM rule and the DAG search never ends. That is why sample IDs may not contain dots.
4. Big intermediate BAMs are `temp()` and are deleted once every consumer has run. The QC step needs them, so it is a consumer too; `--notemp` keeps them for debugging.
5. Each rule names a conda environment, so every tool version is pinned and recorded. Always dry-run (`-n`) first.

Questions:
- Why does Snakemake rerun a step even though its output exists? (An input, the code or the parameters changed; `--rerun-triggers` controls which of these count.)
- Why is the first alignment slow and the second run instant? (Outputs exist and nothing changed: "Nothing to be done".)
- What is the risk of a rule that writes to a path it does not declare as output? (Snakemake cannot track it, so reruns and cleanup go wrong.)

## 7. QC metrics (`make_qc_table.py`)

1. Mapping rate, duplicate fraction, chrM share and blacklist share say how clean the library is; the MAPQ-removed fraction records how much multi-mapping signal the unique-mapping filter discards (ADR-009).
2. NRF is distinct fragments divided by all fragments; PBC1 is fragments seen once divided by distinct fragments; PBC2 is seen-once divided by seen-twice. They measure library complexity: low values mean few molecules were sequenced many times.
3. FRiP is the share of reads inside called peaks: a signal-to-noise measure.
4. TSS enrichment here is read density at TSSs divided by density in flanks 2 kb away. Good ATAC is much denser at promoters. It is a simple version, not the ENCODE score.
5. A flag never excludes a sample: a person looks and signs off. Thresholds sit in the config because public hippocampus data vary.

Questions:
- Why can duplicates be high even when the experiment is fine? (Low input or over-sequencing of a small library.)
- Why do ATAC libraries have so many chrM reads, and why remove them? (Mitochondrial DNA is accessible and abundant; it wastes depth and says nothing about nuclear chromatin.)
- Why a TSS ratio and not a raw count? (Counts scale with depth; the ratio compares two regions in the same library.)

## 8. Integration statistics (`integrate/`)

1. Studies call their own peaks, so effects are compared on a master catalogue: all elements merged where they overlap. When several study peaks fall in one master element, their effects are averaged with inverse-variance weights.
2. Random-effects meta-analysis (DerSimonian-Laird) pools one element's effect across studies and estimates between-study variation (tau squared, I squared). Heterogeneity widens the interval.
3. Concordance asks whether two experiments rank and sign effects alike: Spearman correlation, and the share of same-sign effects with a binomial test.
4. The overlap test asks whether elements changed in aging are over-represented among elements changed by activation, against a null that keeps accessibility (and other covariates) matched.
5. The matched null matters: both studies find highly accessible elements more easily, so a naive overlap looks enriched even when nothing links the two processes. The tests show that effect on synthetic data.

Questions:
- Why not pool raw counts across studies? (Protocol and cell mixture differ; study would be confounded with condition.)
- What does it mean that I squared is 0.9? (About 90% of the variation between study estimates is real heterogeneity, not sampling noise.)
- What would convince you that aging and activation overlap is real and not accessibility? (A significant fold over the matched null, stable across bin choices, and replicated in a second dataset.)
