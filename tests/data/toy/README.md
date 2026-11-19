# Toy dataset (synthetic)

A fully synthetic mouse-like dataset with known truth, built by `make_toy.py` (seed 20261105). CI uses it to test the pipeline and the PRD section 11 control gates without any download. Nothing here is real data.

```
python tests/data/toy/make_toy.py --out tests/data/toy/generated --seed 20261105
```

Generated files go to `tests/data/toy/generated/` (git-ignored, about 22 MB at the default depth of 70,000 read pairs per sample). Only the generator, this README, `samples.tsv` and the tests are committed.

## Design

- Genome: `chr19` (3 Mb, random sequence, GC drifting by 50 kb window) and `chrM` (16.3 kb). Build name `toy`.
- 60 protein-coding genes on chr19. Six carry real IEG names (Fos, Arc, Npas4, Egr1, Junb, Nr4a1); the rest are ToyGene01 to ToyGene54.
- 400 accessible regions: 60 promoters (one per gene) and 340 distal regions of 300-600 bp.
- Truth, stimulated vs control: 48 gained (the six IEG promoters and 42 distal, log2FC 1.2 to 3.0), 12 lost (distal, log2FC -2.5 to -1.2), 340 unchanged.
- Planted motifs (exact counts): AP-1 `TGAGTCA` in 38 of 48 gained regions, CREB `TGACGTCA` in 19, MEF2 `CTAAAAATAG` in 14; each of the three in 27 of the 340 unchanged regions; none in lost regions.
- ATAC: 2 conditions x 3 replicates, paired-end 2 x 50 bp, about 70,000 pairs per sample. Read origin: 66% regions, 4% blacklist pile-ups, 18% uniform background, 12% chrM. 10% exact duplicates.
- RNA: counts for 60 genes x 6 samples (negative binomial, dispersion 0.1); genes next to gained or lost regions change in the same direction.
- Also written: `toy.gtf`, `toy.fa`, `toy.chrom.sizes`, `toy_blacklist.bed` (2 regions), `toy_ccres.bed` (60 PLS and 238 dELS), `truth_regions.tsv`, `truth_genes.tsv`, `MANIFEST.json` (SHA-256 of every file's uncompressed content).

Hashes are stable for a given numpy version; pin numpy in the environment lockfile.

## Reference results (Claude's sandbox: bowtie2 2.5.2, samtools 1.19.2, bedtools 2.31.1, MACS3 3.1.0, R 4.3.3, DESeq2 1.42.0)

| Check | Result |
|---|---|
| Alignment rate | at least 99.98% in every sample |
| Duplicates removed | about 10% of q30 proper pairs |
| chrM share after dedup | 12.0% |
| Peaks per sample | 363 to 395 (design has 400 regions) |
| Consensus peaks (at least 2 replicates) | 396 |
| Planted regions recovered | 48/48 gained, 11/12 lost, 337/340 unchanged |
| DESeq2 (median-of-ratios, no shrinkage) | 48/48 gained called, 11/11 detected lost called, 3/337 unchanged called, 0 wrong direction, median log2FC of unchanged -0.05 |
| DESeq2 under the 9 permuted 3 v 3 label splits | 0 calls in every split (true labels: 59) |
| Contrast: library-size (CPM) normalisation + Welch t-test | 38/337 unchanged called: composition bias, which is why ADR-008 uses median-of-ratios |

If your numbers differ materially, report them to Claude before changing any threshold.

## Tests

```
pytest -q -m "not integration"      # seconds: generator design, truth, determinism
pytest -q -m integration            # a few minutes: real tools on the toy reads
```

The integration tests need bowtie2, samtools, bedtools, MACS3 and R with DESeq2 on PATH together (environment `ncx-test`, TC-004).
