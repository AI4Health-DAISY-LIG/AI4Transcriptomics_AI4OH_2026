# Step 4, Track B: blind variance-stabilising transformation of the full count matrix, exported
# for the genes that scripts/step4_track_b.py needs. Independent analysis, not the reproduction.
#
# Method: DESeq2 varianceStabilizingTransformation(blind = TRUE, fitType = "parametric") on the
# GEO object (27,506 genes, 43 biopsies), with its deposited full-cohort size factors (equal to
# the recomputed ones, results/step3_de/track_a/size_factors.csv). Only the selected rows are
# written; the transformation itself uses every gene.
#
# Inputs
#   data/raw/GSE115650_biopsy.dds.rda.gz                 authors' DESeqDataSet from GEO
#   data/raw/authors/mk.rda                               authors' mapping of Yao 2014 Table S3
#   results/bonus/translator/yao2014_matching.csv         Yao 2014 Table S1, mapped features
#   results/step4_comparison/track_b/panel_genes.csv      Translator panel (step4_track_b.py --panel)
# Outputs (results/step4_comparison/track_b/)
#   vst_selected.csv.gz    VST values, selected genes x 43 biopsies
#   vst_features.csv       per gene: symbol, membership of S1, S3 and the panel
#   vst_info.txt           dispersion trend and counts
# Run from practical/:  pixi run --manifest-path r/pixi.toml Rscript r/vst_track_b.R

suppressMessages(library(DESeq2))
out <- "results/step4_comparison/track_b"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

tmp <- tempfile(fileext = ".rda")
system2("gunzip", c("-c", "data/raw/GSE115650_biopsy.dds.rda.gz"), stdout = tmp)
dds <- get(load(tmp)[1])
stopifnot(identical(dim(dds), c(27506L, 43L)), !is.null(sizeFactors(dds)))
sf <- read.csv("results/step3_de/track_a/size_factors.csv")
stopifnot(isTRUE(all.equal(unname(sizeFactors(dds)),
                           sf$size_factor_deposited[match(colnames(dds), sf$sample)])))
design(dds) <- ~ 1

vsd <- varianceStabilizingTransformation(dds, blind = TRUE, fitType = "parametric")
stopifnot(isTRUE(all.equal(unname(vsd$sizeFactor), unname(sizeFactors(dds)))))
coefs <- attr(dispersionFunction(estimateDispersionsFit(estimateDispersionsGeneEst(dds), fitType = "parametric")), "coefficients")  # the trend a blind VST fits (design ~1)

yao <- read.csv("results/bonus/translator/yao2014_matching.csv")
s1 <- unique(yao$feature[!is.na(yao$feature) & yao$feature != ""])
mk <- get(load("data/raw/authors/mk.rda")[1])
s3 <- intersect(rownames(mk), rownames(dds))
panel <- read.csv(file.path(out, "panel_genes.csv"))
panel_ids <- panel$gene_id[panel$status == "panel"]
missing <- setdiff(c(s1, panel_ids), rownames(dds))
stopifnot(length(missing) == 0)

keep <- unique(c(s1, s3, panel_ids))
symbol <- if ("gene_name" %in% names(rowData(dds))) rowData(dds)[keep, "gene_name"] else NA
features <- data.frame(gene_id = keep, gene_name = symbol, in_S1 = keep %in% s1,
                       in_S3 = keep %in% s3, in_panel = keep %in% panel_ids)
write.csv(features, file.path(out, "vst_features.csv"), row.names = FALSE)
values <- assay(vsd)[keep, ]
write.csv(round(values, 6), gzfile(file.path(out, "vst_selected.csv.gz")))

lines <- c(sprintf("R %s, DESeq2 %s", getRversion(), packageVersion("DESeq2")),
           "varianceStabilizingTransformation(blind = TRUE, fitType = 'parametric'), design ~ 1,",
           "deposited full-cohort size factors (unchanged by the transformation: TRUE)",
           sprintf("dispersion trend: asymptDisp %.5f, extraPois %.5f",
                   coefs[["asymptDisp"]], coefs[["extraPois"]]),
           sprintf("S1 (Yao 2014 Table S1, mapped): %d features; S3 (authors' mk.rda): %d of %d in the counts; panel: %d genes; exported: %d genes x %d biopsies",
                   length(s1), length(s3), nrow(mk), length(panel_ids), length(keep), ncol(values)),
           sprintf("S1 and S3 overlap: %d features", length(intersect(s1, s3))))
writeLines(lines, file.path(out, "vst_info.txt"))
cat(lines, sep = "\n")
