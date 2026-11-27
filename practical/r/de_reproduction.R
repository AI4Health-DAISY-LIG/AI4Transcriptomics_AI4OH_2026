# Step 3, Track A: reproduction of the authors' differential expression (tables S4, S5, S6,
# S8, S10) with their own DESeq2 object from GEO, as close to their code as possible.
#
# Authors' code (TapscottLab/FSHD-Biopsy-RNA-seq @ 90cd343):
#   about_cluster_deseq.Rmd L55-L185   one fit per contrast, on (group + 9 controls) subsets
#   misc_tools.R L320-L359 (do_deseq)  DESeqDataSet(subset, ~ pheno_type); DESeq(); results()
#                                       with defaults; report padj < 0.05 and |log2FC| >= 2
#   makeSE.R L164-L166                  gene filter rowSums > 20, size factors on the full cohort
# A subset keeps colData$sizeFactor, and DESeq() then reuses it ("using pre-existing size
# factors"), so each contrast runs on the full-cohort size factors. This script checks the
# message and the values, before and after DESeq().
#
# Inputs
#   data/raw/GSE115650_biopsy.dds.rda.gz          authors' DESeqDataSet deposited on GEO
#   results/step1_audit/sample_table.csv          GEO groups, to check RNA.Cluster against
# Outputs (results/step3_de/track_a/)
#   de_<contrast>.csv.gz                          all genes, ordered by padj (gzipped: ~27,500 rows each)
#   size_factors.csv                              full-cohort size factors: deposited vs recomputed
#   summary.txt                                   per contrast: samples, size-factor checks,
#                                                 gene counts vs the published targets
#   r_session_info.txt
# Run from practical/:  pixi run --manifest-path r/pixi.toml Rscript r/de_reproduction.R

suppressMessages(library(DESeq2))
out <- "results/step3_de/track_a"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

tmp <- tempfile(fileext = ".rda")
system2("gunzip", c("-c", "data/raw/GSE115650_biopsy.dds.rda.gz"), stdout = tmp)
dds <- get(load(tmp)[1])
stopifnot(identical(dim(dds), c(27506L, 43L)), !is.null(sizeFactors(dds)))

# RNA.Cluster in the authors' object must match the groups deposited on GEO (FSHD biopsies).
geo <- read.csv("results/step1_audit/sample_table.csv", check.names = FALSE)
geo <- geo[match(colnames(dds), geo$title), ]
stopifnot(!anyNA(geo$title))
fshd <- dds$pheno_type == "FSHD"
stopifnot(all(as.numeric(as.character(dds$RNA.Cluster[fshd])) ==
              as.numeric(geo[["char:rnaseq.rank"]][fshd])))

# Full-cohort size factors: deposited vs recomputed from the deposited counts.
sf_full <- sizeFactors(dds)
fresh <- estimateSizeFactors(DESeqDataSetFromMatrix(counts(dds),
                             data.frame(row.names = colnames(dds)), design = ~ 1))
sf <- data.frame(sample = colnames(dds), size_factor_deposited = sf_full,
                 size_factor_recomputed = sizeFactors(fresh))
write.csv(sf, file.path(out, "size_factors.csv"), row.names = FALSE)

gene_name <- if ("gene_name" %in% names(mcols(dds))) mcols(dds)$gene_name else NA

# Contrasts as in about_cluster_deseq.Rmd; targets are the gene rows of the supplement tables
# (S4-S10) or, for the two notebook-only contrasts, the counts written in the notebook.
ctrl <- dds$pheno_type == "Control"
contrasts <- list(
  S4_G1_outliers  = list(idx = colnames(dds) %in% c("01-0037", "32-0016") | ctrl, target = 5680),
  S5_G1           = list(idx = dds$RNA.Cluster == "1" | ctrl, target = 393),
  G1_no_outliers  = list(idx = (dds$RNA.Cluster == "1" | ctrl) &
                           !colnames(dds) %in% c("01-0037", "32-0016"), target = 40),
  S6_G2           = list(idx = dds$RNA.Cluster == "2" | ctrl, target = 303),
  G2_active       = list(idx = colnames(dds) %in% c("01-0030", "32-0007", "32-0013") | ctrl,
                         target = 540),
  S8_G3           = list(idx = dds$RNA.Cluster == "3" | ctrl, target = 816),
  S10_G4          = list(idx = dds$RNA.Cluster == "4" | ctrl, target = 2826)
)

lines <- c(sprintf("R %s, DESeq2 %s (authors: R 3.4.3, DESeq2 1.16.1)",
                   getRversion(), packageVersion("DESeq2")),
           sprintf("full-cohort size factors, max abs difference deposited vs recomputed: %.3g",
                   max(abs(sf$size_factor_deposited - sf$size_factor_recomputed))),
           "")

for (name in names(contrasts)) {
  idx <- contrasts[[name]]$idx
  sub <- DESeqDataSet(dds[, idx], design = ~ pheno_type)   # as do_deseq
  sf_before <- isTRUE(all.equal(sizeFactors(sub), sf_full[colnames(sub)]))
  msgs <- character()
  fit <- withCallingHandlers(DESeq(sub),
                             message = function(m) msgs <<- c(msgs, conditionMessage(m)))
  sf_after <- isTRUE(all.equal(sizeFactors(fit), sf_full[colnames(fit)]))
  stopifnot(identical(resultsNames(fit)[2], "pheno_type_FSHD_vs_Control"))
  res <- results(fit)                                        # defaults, as do_deseq
  res$gene_name <- if (length(gene_name) > 1) gene_name else NA
  res <- res[order(res$padj), ]
  write.csv(as.data.frame(res), gzfile(file.path(out, paste0("de_", name, ".csv.gz"))))

  sig <- which(res$padj < 0.05)
  cut <- which(res$padj < 0.05 & abs(res$log2FoldChange) >= 2)
  n_cut <- length(cut)
  lines <- c(lines,
    sprintf("[%s] samples: FSHD %d, Control %d (%s)", name, sum(sub$pheno_type == "FSHD"),
            sum(sub$pheno_type == "Control"),
            paste(colnames(sub)[sub$pheno_type == "FSHD"], collapse = " ")),
    sprintf("  size factors: full-cohort values before DESeq %s, after %s; pre-existing message: %s",
            sf_before, sf_after, any(grepl("using pre-existing size factors", msgs))),
    sprintf("  genes: %d, with padj: %d, padj < 0.05: %d", nrow(res), sum(!is.na(res$padj)),
            length(sig)),
    sprintf("  padj < 0.05 and |log2FC| >= 2: %d (up %d, down %d); published %d, difference %+d",
            n_cut, sum(res$log2FoldChange[cut] > 0), sum(res$log2FoldChange[cut] < 0),
            contrasts[[name]]$target, n_cut - contrasts[[name]]$target))
}

writeLines(lines, file.path(out, "summary.txt"))
writeLines(capture.output(sessionInfo()), file.path(out, "r_session_info.txt"))
cat(lines, sep = "\n")
