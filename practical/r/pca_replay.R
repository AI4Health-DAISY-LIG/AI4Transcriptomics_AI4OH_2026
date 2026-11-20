# Replay of the paper's Figure 3 (PCA) from the authors' saved rlog matrix.
#
# Inputs
#   data/raw/authors/sanitized.rlg.rda    authors' rlog matrix, TapscottLab/FSHD-Biopsy-RNA-seq
#                                         commit 90cd343, sha256 70b1758029229dbb89c33f3f0abc2dea5165d43f61e0da102566994c2dc18239
#   data/raw/GSE115650_biopsy.dds.rda.gz  GEO's DESeq2 object, used only to link the authors'
#                                         library IDs to GEO titles
# Method (authors' makeManuscripts_supFigs.R L82): DESeq2 plotPCA, i.e. the 500 genes with the
# highest variance of the rlog values, centred, not scaled.
# Outputs (results/step2_qc/)
#   pca_replay_variance.csv   variance explained by PC1 to PC10 (%)
#   pca_replay_scores.csv     PC1 to PC10 per sample, with GEO title and library ID
#   pca_replay_genes.txt      the 500 genes used, for the same-gene comparison in Python
#   r_session_info.txt        R and package versions
# Run from practical/:  pixi run --manifest-path r/pixi.toml Rscript r/pca_replay.R

suppressMessages(library(DESeq2))
out <- "results/step2_qc"
stopifnot(dir.exists(out))

rlg <- get(load("data/raw/authors/sanitized.rlg.rda"))

# The GEO object is a gzipped .rda: decompress it to a temporary file, then load it.
tmp <- tempfile(fileext = ".rda")
status <- system2("gunzip", c("-c", "data/raw/GSE115650_biopsy.dds.rda.gz"), stdout = tmp)
stopifnot(status == 0)
geo <- get(load(tmp))

# Link library IDs to GEO titles through the TopHat input read counts, unique per library;
# the mapped reads and the size factors must then agree as well.
a <- as.data.frame(colData(rlg))
g <- as.data.frame(colData(geo))
stopifnot(!anyDuplicated(a$Input), !anyDuplicated(g$Input))
idx <- match(a$Input, g$Input)
stopifnot(!anyNA(idx),
          identical(as.character(a$Mapped), as.character(g$Mapped[idx])),
          isTRUE(all.equal(a$sizeFactor, g$sizeFactor[idx])))
geo_title <- rownames(g)[idx]

# PCA as plotPCA computes it: top 500 genes by variance of the rlog values, prcomp centred.
mat <- assay(rlg)
rv <- matrixStats::rowVars(mat)
top <- order(rv, decreasing = TRUE)[seq_len(min(500, length(rv)))]
pca <- prcomp(t(mat[top, ]))
pct <- 100 * pca$sdev^2 / sum(pca$sdev^2)

# Cross-check against DESeq2's own plotPCA.
pp <- plotPCA(rlg, intgroup = "pheno_type", ntop = 500, returnData = TRUE)
stopifnot(isTRUE(all.equal(unname(pp$PC1), unname(pca$x[, 1]))),
          isTRUE(all.equal(attr(pp, "percentVar")[1:2], pct[1:2] / 100)))

k <- seq_len(min(10, ncol(pca$x)))
write.csv(data.frame(component = paste0("PC", k), variance_explained_pct = round(pct[k], 4)),
          file.path(out, "pca_replay_variance.csv"), row.names = FALSE)
write.csv(data.frame(geo_title = geo_title, library = colnames(rlg), round(pca$x[, k], 6)),
          file.path(out, "pca_replay_scores.csv"), row.names = FALSE)
writeLines(rownames(mat)[top], file.path(out, "pca_replay_genes.txt"))
writeLines(capture.output(sessionInfo()), file.path(out, "r_session_info.txt"))

cat(sprintf("PC1 %.2f%%  PC2 %.2f%%\n", pct[1], pct[2]))
cat("PC2 extremes:", geo_title[order(abs(pca$x[, 2]), decreasing = TRUE)[1:2]], "\n")