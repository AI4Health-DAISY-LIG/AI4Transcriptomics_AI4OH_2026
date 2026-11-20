# Reproduction of the Figure 3 rlog matrix from the GEO counts, compared with the authors'
# saved matrix and with the replay (r/pca_replay.R).
#
# Inputs
#   data/raw/GSE115650_biopsy_genecounts.csv.gz   GEO raw counts (27,506 genes x 43 samples)
#   data/raw/authors/sanitized.rlg.rda            authors' saved rlog matrix (reference)
#   results/step2_qc/pca_replay_scores.csv        replay: GEO title <-> library ID, PC scores
#   results/step2_qc/pca_replay_variance.csv      replay: variance explained
#   results/step2_qc/pca_replay_genes.txt         replay: the 500 genes behind Figure 3
# Method, as the authors (makeSE.R L164-L181): median-of-ratios size factors, then
# rlog(blind = TRUE) with the default parametric dispersion trend; PCA as plotPCA
# (top 500 genes by variance, centred, not scaled).
# Outputs
#   results/step2_qc/rlog_geo_comparison.txt      size factors, rlog values and PCA vs the authors
#   results/step2_qc/rlog_geo_pca_variance.csv    variance explained by PC1 to PC10 (%)
#   results/step2_qc/rlog_geo_pca_scores.csv      PC1 to PC10 per sample, signs aligned to the replay
#   data/processed/rlog_geo.csv.gz                recomputed rlog matrix (git-ignored)
# Run from practical/:  pixi run --manifest-path r/pixi.toml Rscript r/rlog_from_geo.R

suppressMessages(library(DESeq2))
out <- "results/step2_qc"

counts <- as.matrix(read.csv(gzfile("data/raw/GSE115650_biopsy_genecounts.csv.gz"),
                             row.names = 1, check.names = FALSE))
storage.mode(counts) <- "integer"
replay <- read.csv(file.path(out, "pca_replay_scores.csv"), check.names = FALSE)
replay_genes <- readLines(file.path(out, "pca_replay_genes.txt"))
ref_pct <- read.csv(file.path(out, "pca_replay_variance.csv"))$variance_explained_pct

# Authors' matrix, renamed to GEO titles with the link established by the replay.
rlg <- get(load("data/raw/authors/sanitized.rlg.rda"))
ref_bpv <- attr(rlg, "betaPriorVar")
colnames(rlg) <- replay$geo_title[match(colnames(rlg), replay$library)]
stopifnot(!anyNA(colnames(rlg)), setequal(colnames(rlg), colnames(counts)),
          setequal(rownames(rlg), rownames(counts)))
rlg <- rlg[rownames(counts), colnames(counts)]

# Size factors and rlog from the GEO counts, design ~ 1 as rlog(blind = TRUE) uses.
dds <- DESeqDataSetFromMatrix(counts, colData = data.frame(row.names = colnames(counts)),
                              design = ~ 1)
dds <- estimateSizeFactors(dds)
sf_diff <- max(abs(sizeFactors(dds) - rlg$sizeFactor))
rld <- rlog(dds, blind = TRUE)
new <- assay(rld)
ref <- assay(rlg)
d <- abs(new - ref)

# PCA as plotPCA computes it, with the top 500 genes chosen afresh.
top <- order(matrixStats::rowVars(new), decreasing = TRUE)[1:500]
p <- prcomp(t(new[top, ]))
pct <- 100 * p$sdev^2 / sum(p$sdev^2)
k <- 1:10
ref_x <- as.matrix(replay[match(colnames(new), replay$geo_title), paste0("PC", k)])
signs <- sign(diag(cor(p$x[, k], ref_x)))
x_aligned <- sweep(p$x[, k], 2, signs, "*")
score_cor <- diag(cor(x_aligned, ref_x))

lines <- c(
  sprintf("R %s, DESeq2 %s", getRversion(), packageVersion("DESeq2")),
  sprintf("size factors, max abs difference vs authors: %.3g", sf_diff),
  sprintf("betaPriorVar: recomputed %.6g, authors %s", attr(rld, "betaPriorVar"),
          if (is.null(ref_bpv)) "not stored" else sprintf("%.6g", ref_bpv)),
  sprintf("rlog values, abs difference: max %.4g, median %.4g, 99th percentile %.4g",
          max(d), median(d), quantile(d, 0.99)),
  sprintf("rlog values, correlation over all genes x samples: %.6f", cor(as.vector(new), as.vector(ref))),
  sprintf("top 500 genes shared with the replay: %d / 500", length(intersect(rownames(new)[top], replay_genes))),
  sprintf("PC1 %.2f%% (replay %.2f%%), PC2 %.2f%% (replay %.2f%%)", pct[1], ref_pct[1], pct[2], ref_pct[2]),
  sprintf("score correlation with replay (sign-aligned): PC1 %.6f, PC2 %.6f", score_cor[1], score_cor[2]),
  sprintf("PC2 extremes: %s", paste(colnames(new)[order(abs(p$x[, 2]), decreasing = TRUE)[1:2]], collapse = ", "))
)
writeLines(lines, file.path(out, "rlog_geo_comparison.txt"))
cat(lines, sep = "\n")

write.csv(data.frame(component = paste0("PC", k), variance_explained_pct = round(pct[k], 4)),
          file.path(out, "rlog_geo_pca_variance.csv"), row.names = FALSE)
write.csv(data.frame(geo_title = rownames(x_aligned), round(x_aligned, 6)),
          file.path(out, "rlog_geo_pca_scores.csv"), row.names = FALSE)
dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)
write.csv(round(new, 6), gzfile("data/processed/rlog_geo.csv.gz"))

# Genes whose rlog values differ (max abs difference > 0.1), and both PCAs without them.
gene_max <- apply(d, 1, max)
diff_genes <- rownames(new)[gene_max > 0.1]
dc <- counts[diff_genes, , drop = FALSE]
write.csv(data.frame(gene = diff_genes, max_abs_diff = round(gene_max[diff_genes], 4),
                     total_count = rowSums(dc), zero_samples = rowSums(dc == 0),
                     max_count = apply(dc, 1, max), in_replay_top500 = diff_genes %in% replay_genes),
          file.path(out, "rlog_geo_differing_genes.csv"), row.names = FALSE)
pct_without <- function(m) {
  m <- m[setdiff(rownames(m), diff_genes), ]
  top <- order(matrixStats::rowVars(m), decreasing = TRUE)[1:500]
  s <- prcomp(t(m[top, ]))$sdev^2
  100 * s / sum(s)
}
p_new <- pct_without(new)
p_ref <- pct_without(ref)
extra <- c(
  sprintf("genes with max abs rlog difference > 0.1: %d (rlog_geo_differing_genes.csv)", length(diff_genes)),
  sprintf("PCA without them: PC1 %.2f%% vs authors %.2f%%, PC2 %.2f%% vs authors %.2f%%",
          p_new[1], p_ref[1], p_new[2], p_ref[2]))
write(extra, file.path(out, "rlog_geo_comparison.txt"), append = TRUE)
cat(extra, sep = "\n")