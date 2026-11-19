# NeuroChrom-X Claude pack v3

Built 15 Nov 2026 for PRD v1.1. This is everything Claude could build and test without your laptop, plus the exact steps for when you have it. Tick the boxes as you go.

## 1. Where things stand

**Written and tested in Claude's sandbox** (Python 3.12, bowtie2 2.5.2, samtools 1.19.2, bedtools 2.31.1, MACS3 3.1.0, fastp 0.23.4, deepTools 4.0.0, Snakemake 9.27, R 4.3.3, DESeq2 1.42.0):

| Piece | Tests |
|---|---|
| Sample-sheet validator, toy sample sheet | 26 |
| Synthetic toy dataset generator (3 Mb genome, 6 ATAC samples, RNA counts, known truth) | 12 |
| Intervals, promoters, enhancer tiers, gene linking, motif scanning, background matching, motif enrichment | 15 |
| Sequence helper, MEME reader, workflow scripts (consensus peaks, expressed genes, group sheet, annotation, motif CLI) | 7 |
| Integration statistics: effect tables, master catalogue, concordance, matched overlap test, random-effects meta-analysis | 8 |
| **Fast tests in total** | **68** |
| ATAC chain with the real tools, called directly | 4 integration |
| **The same chain run through Snakemake (68 jobs: toy data, fastp, Bowtie2, filtering, duplicates, MACS3, consensus, counts, DESeq2, annotation, motifs, tracks, QC table)** | 6 integration |

Result on the toy data, through Snakemake: all 48 planted gained regions recovered and called, 3 false calls among 337 unchanged regions, no wrong-direction calls, AP-1 the top enriched motif (odds ratio about 22), all samples pass QC, and a second run does nothing. `ruff` and `black` are clean. Reference numbers are in `tests/data/toy/README.md`.

**Written but not run anywhere:** the 14 ADRs, the gene sets (not yet cross-checked against the published tables), the conda environment files, the Picard, FastQC, MultiQC, apeglm and single-end code paths, and the task cards.

**Not built yet (known gaps):** ChIP, RNA-alignment and methylation branches; SRA/ENA download rules; BAM, peak and counts entry points for ATAC; HOMER as motif engine; the background-bin normalization run (ADR-008 sensitivity); IGV session export; the BioAge-X export; the advanced track. They come after the first real-data run.

## 2. Calendar: the repo is reviewed in the last week of December or the first week of January

Plan for about 10 focused hours a week. Quality comes before speed, and the repo should be review-ready by about 20 Dec, then frozen for the review window.

| Week | Goal | Cards |
|---|---|---|
| 16-22 Nov | Laptop set up; repo with all tests green; machine audit; environments | TC-001, TC-002, TC-004 |
| 23-29 Nov | BioAge-X audit; real metadata for A1 and B1 confirmed; reference files | TC-003, TC-005, TC-006 |
| 30 Nov-6 Dec | Workflow verified with conda and on every code path; first real ATAC run on A1 (Claude writes TC-008 after TC-005 and TC-006) | TC-007, TC-008 |
| 7-13 Dec | A1 analysis and gates on real data; RNA branch for A1; figures | TC-009 onward |
| 14-20 Dec | B1 aging data; the activity-versus-aging overlap test; README with results; release v0.1.0; CI badge | |
| 21 Dec-7 Jan | Review window: main is frozen, only bug fixes. Extra work (histone, methylation, BioAge-X bridge) goes on branches | |

If you have fewer hours, cut from the bottom: the toy benchmark, the A1 real-data run, the README and the PRD already make a solid repository. The application text due on 15 Nov may only claim what is true on the day it is sent (PRD section 5, claims ladder).

## 3. Do now, from your phone (30-60 minutes)

- [ ] Create the Claude Project "NeuroChrom-X". Export the PRD (v1.1) as Markdown or PDF and add it to the Project knowledge. Paste Appendix A of the PRD into the Project instructions.
- [ ] Optional: create an empty private GitHub repository named `neurochrom-x` (no README, no licence). It is your backup. The repo stays private until you decide to publish.
- [ ] Send Claude the GitHub link of BioAge-X. If it is public, Claude can audit it without waiting for the laptop.
- [ ] Read PRD sections 1 to 5, then `docs/explain-it/EXPLAIN_IT.md` sections 1 to 3. Answer each question in a notes file in your own words.
- [ ] Decide your laptop's operating system. On Windows, plan to use WSL2 with Ubuntu: Bowtie2, Bismark and MACS3 are Linux and macOS tools. Everything below assumes a Linux or macOS shell. The GTX 1650 can be used from WSL2 once the NVIDIA driver is installed, but nothing in the plan needs it yet.

## 4. Laptop day 1 (about 3 hours)

### Step 0: prerequisites (30-45 minutes)

- [ ] Windows only: in an administrator PowerShell run `wsl --install`, restart, open Ubuntu, and work inside it. Keep the repo in the Linux home folder (`~/neurochrom-x`), not under `/mnt/c`.
- [ ] Install Git and Python 3.11 or later. Check: `git --version` and `python3 --version`.
- [ ] Install Miniforge (conda with the conda-forge channel, plus mamba) following the instructions on the Miniforge GitHub page. Check `mamba --version`.
- [ ] Install Antigravity and sign in. VS Code is optional.
- [ ] Note your laptop's cores, RAM and free disk (Linux: `nproc`, `free -g`, `df -h .`).

### Step 1: create the repo with the pack in it (10 minutes)

```
mkdir -p ~/neurochrom-x && cd ~/neurochrom-x
git init -b main
unzip ~/Downloads/neurochrom-x-claude-pack.zip -d /tmp/ncx
cp -r /tmp/ncx/neurochrom-x-claude-pack/. .
```

- [ ] Export the PRD from the Claude document as Markdown and save it as `docs/PRD.md`.
- [ ] Check: `ls` shows `AGENTS.md`, `config`, `docs`, `src`, `tests`, `workflow`.
- [ ] Commit: `git add -A && git commit -m "chore: add Claude pack and PRD"`

### Step 2: TC-001, the scaffold (about 15 minutes of agent time)

- [ ] Open `~/neurochrom-x` in Antigravity. Start a new task and paste the full text of `docs/task-cards/TC-001.md`.
- [ ] Check the Run Report yourself: acceptance test 3 must show **68 passed**, test 4 must print `OK: 12 samples, 1 study(ies), 4 group(s) (schema 0.1)`, and test 1 must show no modified pack files.
- [ ] Paste the Run Report to Claude. Do not merge before Claude says Accept.
- [ ] If tests fail, paste the first error with 20 lines of context. Do not let the agent edit tests or pack files to make them pass.

### Step 3: TC-002, machine audit and environments (60-90 minutes, mostly waiting)

- [ ] Paste `docs/task-cards/TC-002.md`. The environment files are already written; the agent creates them and records the versions. Plan for roughly 10-15 GB of free disk.
- [ ] Paste the Run Report and the audit JSON to Claude. Claude reads the recommendations (RNA path STAR or Salmon, FASTQ or processed inputs) and tells you what they mean.

### Step 4: TC-004, verify the toy dataset and the real-tool tests (30-45 minutes)

- [ ] Paste `docs/task-cards/TC-004.md`. It creates the `ncx-test` environment and runs everything on your machine.
- [ ] Check the Run Report: 68 fast tests pass, 10 integration tests pass, the two MANIFEST hash maps match, and the numbers match `tests/data/toy/README.md`.
- [ ] If a number differs materially, paste both to Claude. Never change a threshold yourself.

## 5. Laptop days 2 and 3

### Step 5: TC-003, the BioAge-X audit (30-60 minutes)

- [ ] `git clone` your BioAge-X repository next to this one. Put the absolute path on the `BIOAGEX_PATH:` line of the card text when you paste it.
- [ ] Check that BioAge-X was not changed: `git -C <path> status --porcelain` prints nothing.
- [ ] Paste the Run Report and `docs/bioagex_interface.md` to Claude. It settles the species question (OD6) and the export adapter.

### Step 6: TC-005, dataset registry for A1 and B1 (30-45 minutes, needs internet)

- [ ] Paste `docs/task-cards/TC-005.md`. It reads GEO and ENA metadata only and downloads no reads.
- [ ] Paste the Run Report and `docs/registry_review_A1_B1.md` to Claude. Claude checks every claim in PRD section 6 before any data are downloaded.

### Step 7: TC-006, reference resources (30-60 minutes plus downloads)

- [ ] Paste `docs/task-cards/TC-006.md`. Expect roughly 1-2 GB of downloads and a Bowtie2 index build of 10-60 minutes.
- [ ] Paste the Run Report and `docs/resources_report.md` to Claude. Claude checks unmapped gene symbols and freezes the gene sets after the cross-check against the published tables.

### Step 8: TC-007, verify the workflow with conda (1-2 hours)

- [ ] Paste `docs/task-cards/TC-007.md`. Claude already ran the whole toy workflow in a sandbox; this card runs it with the real conda environments and tests the paths Claude could not (Picard, apeglm, FastQC and MultiQC, single-end reads).
- [ ] Paste the Run Report, the QC table and the DESeq2 summary to Claude.

## 6. What is in the zip

```
AGENTS.md                       rules for Antigravity (PRD Appendix B)
README.md                       this file
config/config.yaml              every PRD section 7 default; overrides/toy.yaml for CI
config/gene_sets/               iegs_v1.tsv (49), plasticity_v1.tsv (37), README
docs/adr/                       ADR-001 to ADR-014, template, index
docs/explain-it/EXPLAIN_IT.md   study notes and interview questions per module (8 sections)
docs/task-cards/                TC-001 to TC-007 and the card template
docs/AI_ASSISTANCE.md           pre-filled log; you fill in "what I verified"
src/neurochromx/                intervals, sequences, io/samplesheet, annotate/, motifs/, integrate/
tests/                          68 fast tests, 10 integration tests, toy dataset generator
workflow/Snakefile, rules/      the ATAC branch (nine rule files)
workflow/envs/                  ten conda environment files (no version pins; TC-002 records them)
workflow/scripts/               consensus peaks, counts, annotation, motifs, QC table, DESeq2
```

## 7. PRD changes since v1.0

Items 1 and 2 and the status log are already applied in the Claude document (PRD v1.1; export it again before adding it to the Project knowledge). The rest are not applied yet; do them when you next edit the PRD.

1. **Negative-control gate (section 11), refined.** Permuted labels use all distinct label splits (at most 20; a 3 v 3 design has 9 besides the true one) and must give a median DAR count of at most 5% of the true-label count. Where the truth is known, at most 5% of unchanged elements may be called. Reason: with 3 to 4 replicates per group, "20 shuffles" cannot exist, and shuffled splits still share part of the true signal.
2. **Normalization sanity (same row).** The median log2FC of unchanged elements is within 0.15 under DESeq2's median-of-ratios. Evidence: on the toy data library-size normalization called 38 of 337 unchanged peaks, DESeq2 called 3.
3. **Section 9, `regulatory_elements`.** Enhancer tiers are not exclusive, so the table has `class` (promoter, enhancer, other) plus the flags `tier1`, `tier2` and `tier3`, instead of the classes enhancer_T1 to enhancer_T3.
4. **Section 9, sample IDs** contain letters, digits, `_` and `-` only. Dots clash with file suffixes in the rule patterns.
5. **Section 7, QC table.** It reports mapping rate, MAPQ-removed fraction (ADR-009), duplicate fraction, chrM and blacklist fractions, FRiP, NRF, PBC1, PBC2 and a simple TSS enrichment (read density at TSSs over density in flanks 2 kb away; not the ENCODE score).
6. **Section 5, roadmap.** P1 is shorter: the workflow, validator, toy benchmark and integration statistics exist and are tested. TC-007 verifies them in the real environments.
7. **Section 10, layout.** Adds `workflow/scripts/{py,R}/`, `tests/integration/`, `tests/conftest.py`, `docs/explain-it/`, `config/overrides/`, `src/neurochromx/{intervals,sequences}.py`, `motifs/{scan,meme}.py` and `integrate/`.

## 8. If something goes wrong

- A test or command fails: paste the exact command, the first error and about 20 lines of context. Do not edit tests, thresholds or pack files to make it pass.
- Antigravity edits a protected file: `git checkout -- <file>`, then tell it that the Task Card forbids it.
- conda is slow to solve: use `mamba`, and run TC-002 when you can leave the laptop alone.
- Snakemake looks stuck on a dry run: the cause is usually an ambiguous wildcard pattern. Paste the Snakefile edit you made.
- Out of disk: stop, run `df -h .`, and send Claude the audit. Real data are far larger than the toy set (about 22 MB).
- Integration tests fail by a small margin on your machine: paste both sets of numbers. Different tool versions can shift counts slightly; Claude decides whether a threshold needs review.

## 9. What to paste back to Claude each time

1. The card id and the Run Report (changed files, commands with exit codes, test results, metrics, deviations, errors, open questions, disk and time used).
2. Any file the card names as an output for review (for example `docs/bioagex_interface.md`).
3. Your own one-line verdict: did the acceptance tests you checked pass?

Claude answers Accept, Changes or Re-spec, plus a short explain-it note for each module.
