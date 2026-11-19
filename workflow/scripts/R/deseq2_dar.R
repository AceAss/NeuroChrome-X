#!/usr/bin/env Rscript
# Differential accessibility (or any count-based contrast) with DESeq2.
# PRD section 7 and ADR-008: median-of-ratios normalisation, design ~ batch + condition,
# apeglm shrinkage, FDR 0.05 and |log2FC| >= 0.585 as defaults (all configurable).
#
# Rscript deseq2_dar.R --counts counts.tsv --samples samples.tsv \
#   --treatment stimulated --control control --out dar.tsv [--size_factors sf.tsv]
#
# counts.tsv : element_id, chrom, start, end, <one column per sample>
# samples.tsv: sample_id, condition, batch (other columns ignored; batch optional)
# A batch column with a single level is dropped from the design (it would be collinear).
suppressPackageStartupMessages({
  library(optparse)
  library(DESeq2)
})

option_list <- list(
  make_option("--counts", type = "character"),
  make_option("--samples", type = "character"),
  make_option("--treatment", type = "character"),
  make_option("--control", type = "character"),
  make_option("--out", type = "character"),
  make_option("--shrink", type = "character", default = "apeglm",
              help = "apeglm or none [default %default]"),
  make_option("--fdr", type = "double", default = 0.05),
  make_option("--min_abs_log2fc", type = "double", default = 0.585),
  make_option("--size_factors", type = "character", default = NULL,
              help = "optional TSV (sample_id, size_factor), e.g. from background bins"),
  make_option("--seed", type = "integer", default = 20261105L)
)
opt <- parse_args(OptionParser(option_list = option_list))
for (name in c("counts", "samples", "treatment", "control", "out")) {
  if (is.null(opt[[name]])) stop("missing required option --", name)
}
set.seed(opt$seed)

counts <- read.delim(opt$counts, check.names = FALSE, stringsAsFactors = FALSE)
meta <- read.delim(opt$samples, check.names = FALSE, stringsAsFactors = FALSE)
anno <- counts[, c("element_id", "chrom", "start", "end")]
mat <- as.matrix(counts[, setdiff(colnames(counts), colnames(anno)), drop = FALSE])
rownames(mat) <- counts$element_id

meta <- meta[meta$sample_id %in% colnames(mat) &
               meta$condition %in% c(opt$treatment, opt$control), , drop = FALSE]
mat <- mat[, meta$sample_id, drop = FALSE]
storage.mode(mat) <- "integer"
meta$condition <- relevel(factor(meta$condition), ref = opt$control)
terms <- "condition"
if ("batch" %in% colnames(meta) && length(unique(meta$batch)) > 1) {
  meta$batch <- factor(meta$batch)
  terms <- c("batch", terms)
}
design <- as.formula(paste("~", paste(terms, collapse = " + ")))
message("design: ", deparse(design), "; samples: ", nrow(meta))

dds <- DESeqDataSetFromMatrix(countData = mat, colData = meta, design = design)
if (!is.null(opt$size_factors)) {
  sf <- read.delim(opt$size_factors, stringsAsFactors = FALSE)
  sizeFactors(dds) <- sf$size_factor[match(meta$sample_id, sf$sample_id)]
} else {
  dds <- estimateSizeFactors(dds)
}
dds <- DESeq(dds, quiet = TRUE)
coef <- paste0("condition_", opt$treatment, "_vs_", opt$control)
res <- results(dds, contrast = c("condition", opt$treatment, opt$control))
lfc_mle <- res$log2FoldChange
lfc <- lfc_mle
if (opt$shrink == "apeglm") {
  if (!requireNamespace("apeglm", quietly = TRUE)) {
    stop("apeglm is not installed; install it or pass --shrink none")
  }
  lfc <- as.data.frame(lfcShrink(dds, coef = coef, type = "apeglm", quiet = TRUE))$log2FoldChange
} else if (opt$shrink != "none") {
  stop("--shrink must be apeglm or none")
}

out <- data.frame(
  anno, baseMean = res$baseMean, log2FC_mle = lfc_mle, log2FC = lfc,
  lfcSE = res$lfcSE, stat = res$stat, pvalue = res$pvalue, padj = res$padj
)
out$call <- ifelse(!is.na(out$padj) & out$padj < opt$fdr & abs(out$log2FC) >= opt$min_abs_log2fc,
                   ifelse(out$log2FC > 0, "up", "down"), "ns")
write.table(out, opt$out, sep = "\t", quote = FALSE, row.names = FALSE)
message("wrote ", nrow(out), " rows; calls: ", paste(names(table(out$call)), table(out$call), collapse = ", "))
writeLines(c(paste("R", getRversion()), paste("DESeq2", packageVersion("DESeq2"))),
           paste0(opt$out, ".versions"))
