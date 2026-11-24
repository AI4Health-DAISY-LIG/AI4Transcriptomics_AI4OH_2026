# Step 4, Track A: reproduction of the authors' biomarker score, RNA groups, Figure 1 calls,
# Figure 2 and the MRI/pathology correlations, with their own DESeq2 object from GEO.
#
# Authors' code (TapscottLab/FSHD-Biopsy-RNA-seq @ 90cd343):
#   misc_tools.R L84-L129 (getRankScore)   rlog(fitType = "local") on the four-marker subset;
#                                            sum of (rlog - control mean); hclust on the rlog
#                                            values; calls = markers with relative rlog > 2
#   misc_tools.R L185-L198 (getRankSummary) table sorted by the sum; GEO score = sum - minimum
#   misc_tools.R L350-L351 (do_deseq)       the four marker IDs
#   addColDataToDDS.R L14-L31               relabelling of the raw clusters into groups 1-4
#   misc_tools.R L335-L348 (do_deseq)       DE call for the DUX4 candidates: the 55 features of
#                                            data/mk.rda intersected with the DE results, then
#                                            padj < 0.05 and |log2FC| >= 2. TRIM43CP is not in the
#                                            deposited counts, which leaves 54 (its explicit removal
#                                            in rankScore.R L25 concerns the clustering only)
#   makeManuscripts_supFig_scatterPlot.R    Figure 1 scatter panel: coloured by the final RNA group
#
# Inputs
#   data/raw/GSE115650_biopsy.dds.rda.gz          authors' DESeqDataSet deposited on GEO
#   data/raw/GSE115650_biopsy_FPKM.csv.gz          GEO FPKM (used if the object has no FPKM assay)
#   data/raw/authors/mk.rda                        authors' DUX4 candidate features (55)
#   results/step1_audit/sample_table.csv           GEO groups, scores, MRI and pathology
#   results/step3_de/track_a/size_factors.csv      full-cohort size factors (Step 3)
#   results/step3_de/track_a/de_S8_G3.csv.gz, de_S10_G4.csv.gz   Track A DE (Step 3)
# Outputs (results/step4_comparison/track_a/)
#   score.csv              per biopsy: sum, GEO-convention score, GEO score, difference, calls,
#                          raw cluster, recomputed group, GEO group
#   groups_crosstab.txt    raw cut and recomputed groups against GEO
#   controls_fpkm.csv      FPKM of the four markers in the nine controls
#   figure2.csv            score by MRI status: n, mean, population and sample SD
#   correlations.csv       Spearman correlations with the paper's values, by subset
#   dux4_targets_de.csv    DE status of each DUX4 candidate in groups 3 and 4
#   figure1_style.png      Figure 1 scatter panel: sum by biopsy, colour = recomputed group,
#                          dot size = calls + 1 (the heat map panel is not redrawn)
#   summary.txt, r_session_info.txt
# Run from practical/:  pixi run --manifest-path r/pixi.toml Rscript r/score_reproduction.R

suppressMessages(library(DESeq2))
out <- "results/step4_comparison/track_a"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
lines <- sprintf("R %s, DESeq2 %s (authors: R 3.4.x, DESeq2 1.16.1)",
                 getRversion(), packageVersion("DESeq2"))
say <- function(...) lines <<- c(lines, sprintf(...))
num <- function(v) suppressWarnings(as.numeric(as.character(v)))

tmp <- tempfile(fileext = ".rda")
system2("gunzip", c("-c", "data/raw/GSE115650_biopsy.dds.rda.gz"), stdout = tmp)
dds <- get(load(tmp)[1])
stopifnot(identical(dim(dds), c(27506L, 43L)))
sf <- read.csv("results/step3_de/track_a/size_factors.csv")
stopifnot(isTRUE(all.equal(unname(sizeFactors(dds)),
                           sf$size_factor_deposited[match(colnames(dds), sf$sample)])))
geo <- read.csv("results/step1_audit/sample_table.csv", check.names = FALSE)
geo <- geo[match(colnames(dds), geo$title), ]
stopifnot(!anyNA(geo$title))
ctrl <- dds$pheno_type == "Control"
fshd <- !ctrl
say("biopsies: %d FSHD, %d controls; size factors equal to the Step 3 values: TRUE",
    sum(fshd), sum(ctrl))

# --- 1.1 Score -------------------------------------------------------------------------------
ids <- c(PRAMEF2 = "ENSG00000120952.4", TRIM43 = "ENSG00000144015.4",
         KHDC1L = "ENSG00000256980.4", LEUTX = "ENSG00000213921.6")
stopifnot(all(ids %in% rownames(dds)))
symbols <- if ("gene_name" %in% names(rowData(dds))) rowData(dds)[ids, "gene_name"] else NA
say("markers: %s; symbols in the object: %s; match: %s",
    paste(ids, collapse = " "), paste(symbols, collapse = " "),
    isTRUE(all(symbols == names(ids))))

x <- dds[rownames(dds) %in% ids, ]                      # as getRankScore: subset, dds row order
lg <- tryCatch(assay(rlog(x, fitType = "local")),
               error = function(e) stop("rlog on the four-marker subset failed: ",
                                        conditionMessage(e)))
relative <- lg - rowMeans(lg[, ctrl])                   # rlog minus control mean, per gene
rownames(relative) <- paste0("relative_", names(ids)[match(rownames(relative), ids)])
log_sum <- colSums(relative)
score <- log_sum - min(log_sum)                         # GEO convention: minimum over all 43
calls <- colSums(relative > 2)
geo_score <- num(geo[["char:rnaseq.score"]])
diff <- score - geo_score
same_ranks <- identical(unname(rank(score[fshd])), unname(rank(geo_score[fshd])))
max_diff <- max(abs(diff[fshd]))
say("score: minimum of the sum at %s (%.5f); FSHD vs GEO: max |difference| %.2e, identical ranks %s, within 1e-3: %s",
    colnames(dds)[which.min(log_sum)], min(log_sum), max_diff, same_ranks, max_diff < 1e-3)
say("score: controls, GEO-convention values %s (GEO gives no score for controls)",
    paste(sprintf("%.3f", score[ctrl]), collapse = " "))

# --- 1.2 Groups ------------------------------------------------------------------------------
raw <- cutree(hclust(dist(t(lg))), k = 4)               # Euclidean, complete, on the rlog values
level <- rank(tapply(log_sum, raw, median))             # raw clusters ordered by their median sum
# Authors' rule, in their order, on the table sorted by the sum (addColDataToDDS.R L20-L26):
# cluster 2 -> 1; rows 20-24 -> 2; cluster 1 -> 3; clusters 3 and 4 -> 4 (the last two steps
# would overwrite a row 20-24 that is not in cluster 2). Their raw table is not public, so each
# raw cluster is also identified by its level (lowest -> 1, second -> 3, top two -> 4); the
# literal numeric rule shows whether the cluster numbering matches theirs.
sorted <- order(log_sum)
relabel <- function(cluster_sorted, as_g1, as_g3, as_g4) {
  g <- rep(1, length(cluster_sorted))
  g[cluster_sorted %in% as_g1] <- 1
  g[20:24] <- 2
  g[cluster_sorted %in% as_g3] <- 3
  g[cluster_sorted %in% as_g4] <- 4
  out <- numeric(length(g))
  out[sorted] <- g                                     # back to the object's column order
  out
}
literal <- relabel(raw[sorted], 2, 1, c(3, 4))
by_level <- relabel(level[as.character(raw[sorted])], 1, 2, c(3, 4))
geo_group <- num(geo[["char:rnaseq.rank"]])
dep_group <- num(as.character(dds$RNA.Cluster))
say("groups: raw cluster sizes %s; recomputed groups equal GEO for the 34 FSHD: %s; controls all in group 1: %s; equal to RNA.Cluster for all 43: %s; literal numeric rule gives the same: %s",
    paste(table(raw), collapse = "/"), identical(by_level[fshd], geo_group[fshd]),
    all(by_level[ctrl] == 1), identical(by_level, dep_group), identical(literal, by_level))
say("groups: sizes among FSHD (1-4) %s; paper 10/5/13/6",
    paste(tabulate(by_level[fshd], 4), collapse = "/"))
writeLines(c("raw cut (rows) vs GEO group (columns; 0 = control):",
             capture.output(table(raw = raw, geo = ifelse(ctrl, 0, geo_group))), "",
             "recomputed group (rows) vs GEO group (columns; 0 = control):",
             capture.output(table(recomputed = by_level, geo = ifelse(ctrl, 0, geo_group)))),
           file.path(out, "groups_crosstab.txt"))

write.csv(data.frame(sample = colnames(dds), condition = dds$pheno_type, batch = dds$Batch,
                     log_sum = log_sum, score_recomputed = score, score_geo = geo_score,
                     difference = diff, calls = calls, raw_cluster = raw,
                     group_recomputed = by_level, group_geo = geo_group,
                     t(round(relative, 5)), check.names = FALSE)[sorted, ],
          file.path(out, "score.csv"), row.names = FALSE)

# --- 1.3 Figure 1-style plot -----------------------------------------------------------------
draw <- function() {
  par(mar = c(6, 4, 2, 1))
  plot(seq_along(sorted), log_sum[sorted], col = by_level[sorted] + 1, pch = 19,
       cex = 0.5 * (calls[sorted] + 1), xaxt = "n", xlab = "",
       ylab = "Sum of relative expression (rlog - control mean)",
       main = "Figure 1 scatter panel, recomputed (dot size = markers above 2, plus 1)")
  axis(1, at = seq_along(sorted), cex.axis = 0.5, las = 2,
       labels = paste0(ifelse(ctrl[sorted], "", "*"), colnames(dds)[sorted]))
  legend("topleft", legend = paste("group", 1:4), col = 2:5, pch = 19, cex = 0.7)
}
ok <- tryCatch({ png(file.path(out, "figure1_style.png"), width = 1600, height = 900, res = 150,
                     type = "cairo"); TRUE }, error = function(e) FALSE)
if (!ok) pdf(file.path(out, "figure1_style.pdf"), width = 11, height = 6)
draw()
invisible(dev.off())

# --- 1.4 Biomarkers in the controls ----------------------------------------------------------
fpkm <- if ("FPKM" %in% assayNames(dds)) {
  src <- "FPKM assay of the GEO object"; assay(dds, "FPKM")[ids, ]
} else {
  src <- "GEO FPKM file"
  f <- as.matrix(read.csv(gzfile("data/raw/GSE115650_biopsy_FPKM.csv.gz"), row.names = 1,
                          check.names = FALSE))
  f[ids, colnames(dds)]
}
cf <- t(fpkm[, ctrl])
colnames(cf) <- names(ids)
write.csv(data.frame(sample = rownames(cf), cf, check.names = FALSE),
          file.path(out, "controls_fpkm.csv"), row.names = FALSE)
say("controls (%s): max FPKM per marker %s; values >= 0.1: %d of %d",
    src, paste(sprintf("%s %.3f", names(ids), apply(cf, 2, max)), collapse = ", "),
    sum(cf >= 0.1), length(cf))

# --- 1.5 and 1.6 MRI status ------------------------------------------------------------------
stir <- num(geo[["char:stir"]])
t1 <- num(geo[["char:t1"]])
pos <- fshd & !is.na(stir) & stir > 0
normal <- fshd & !is.na(stir) & stir == 0 & !is.na(t1) & t1 == 0
in_g24 <- by_level %in% 2:4
no16 <- colnames(dds) != "32-0016"
say("STIR+ biopsies in groups 2-4 (the paper's definition of expressing DUX4 targets): %d/%d; without 32-0016: %d/%d (paper 20/23 and 91%%)",
    sum(pos & in_g24), sum(pos), sum(pos & in_g24 & no16), sum(pos & no16))
say("STIR+ biopsies with at least one marker above relative rlog 2 (Figure 1 calls >= 1), another definition: %d/%d",
    sum(pos & calls >= 1), sum(pos))
say("Figure 2 groups: MRI-normal (STIR 0 and T1 0) %d, STIR+ %d, in neither: %s",
    sum(normal), sum(pos), paste(colnames(dds)[fshd & !pos & !normal], collapse = " "))
sdp <- function(v) sqrt(mean((v - mean(v))^2))
fig2 <- do.call(rbind, lapply(c("score_geo", "score_recomputed"), function(s) {
  v <- if (s == "score_geo") geo_score else score
  do.call(rbind, lapply(list(MRI_normal = normal, STIR_positive = pos), function(i)
    data.frame(score = s, n = sum(i), mean = mean(v[i]), sd_population = sdp(v[i]),
               sd_sample = sd(v[i]))))
}))
fig2 <- cbind(group = rep(c("MRI_normal", "STIR_positive"), 2), fig2,
              paper_mean = c(1.4, 7.7), paper_sd = c(2.4, 5.7))
write.csv(fig2, file.path(out, "figure2.csv"), row.names = FALSE)
welch <- t.test(geo_score[pos], geo_score[normal])$p.value
stopifnot(!anyDuplicated(geo_score[pos | normal]))       # no ties: the exact test applies
wilc_exact <- wilcox.test(geo_score[pos], geo_score[normal], exact = TRUE)$p.value
wilc <- suppressWarnings(wilcox.test(geo_score[pos], geo_score[normal], exact = FALSE,
                                     correct = TRUE)$p.value)
say("Figure 2 (GEO score): means %.3f / %.3f, population SD %.4f / %.4f, sample SD %.4f / %.4f (paper 1.4 / 7.7, SD 2.4 / 5.7)",
    fig2$mean[1], fig2$mean[2], fig2$sd_population[1], fig2$sd_population[2],
    fig2$sd_sample[1], fig2$sd_sample[2])
say("Figure 2 test, unnamed in the paper (P < 0.001); candidates, two-sided, none confirmed: Welch t-test p = %.2g; exact Wilcoxon rank-sum p = %.2g; Wilcoxon normal approximation with continuity correction p = %.2g",
    welch, wilc_exact, wilc)

# --- 1.7 Spearman correlations ---------------------------------------------------------------
endpoints <- c(stir = "char:stir", t1 = "char:t1", fat_fraction = "char:fat.fraction",
               pathology = "char:path.score")
paper <- list(stir = c(score = 0.61, group = 0.63), t1 = c(score = 0.49, group = 0.47),
              fat_fraction = c(score = 0.34, group = 0.32), pathology = c(score = 0.63))
path_ok <- !is.na(num(geo[["char:path.score"]]))
subsets <- list(all = fshd, without_32_0016 = fshd & no16, pathology_complete = fshd & path_ok)
rna <- list(score_recomputed = score, score_geo = geo_score, group = geo_group)
corr <- list()
for (e in names(endpoints)) for (r in names(rna)) for (s in names(subsets)) {
  y <- num(geo[[endpoints[[e]]]])
  i <- subsets[[s]] & !is.na(y) & !is.na(rna[[r]])
  ct <- suppressWarnings(cor.test(rna[[r]][i], y[i], method = "spearman", exact = FALSE))
  key <- if (r == "group") "group" else "score"
  corr[[length(corr) + 1]] <- data.frame(
    endpoint = e, rna = r, subset = s, n = sum(i), rho = unname(ct$estimate), p = ct$p.value,
    paper_rho = if (key %in% names(paper[[e]])) paper[[e]][[key]] else NA)
}
corr <- do.call(rbind, corr)
write.csv(corr, file.path(out, "correlations.csv"), row.names = FALSE)
say("Spearman (average ranks for ties; p from the asymptotic t approximation), all endpoint-complete FSHD:")
main <- corr[corr$subset == "all" & corr$rna != "score_recomputed", ]
lines <- c(lines, sprintf("  %-12s %-9s n %2d  rho %.3f  p %.2g  paper %s", main$endpoint, main$rna,
                          main$n, main$rho, main$p, ifelse(is.na(main$paper_rho), "-",
                                                           sprintf("%.2f", main$paper_rho))))

# --- 1.8 DUX4 candidates in groups 3 and 4 ---------------------------------------------------
mk <- get(load("data/raw/authors/mk.rda")[1])
cand <- rownames(mk)
cand_names <- if ("gene_name" %in% names(rowData(mk))) rowData(mk)$gene_name else NA
trim43cp <- "ENSG00000144188.9"
de_tab <- data.frame(gene_id = cand, gene_name = cand_names, in_counts = cand %in% rownames(dds),
                     is_trim43cp = cand == trim43cp)
for (contrast in c("S8_G3", "S10_G4")) {
  res <- read.csv(gzfile(sprintf("results/step3_de/track_a/de_%s.csv.gz", contrast)),
                  row.names = 1)
  hit <- res[cand, ]
  de_tab[[paste0("log2FC_", contrast)]] <- hit$log2FoldChange
  de_tab[[paste0("padj_", contrast)]] <- hit$padj
  de_tab[[paste0("de_", contrast)]] <- ifelse(de_tab$in_counts,             # NA = not assessed
    !is.na(hit$padj) & hit$padj < 0.05 & abs(hit$log2FoldChange) >= 2, NA)
}
write.csv(de_tab, file.path(out, "dux4_targets_de.csv"), row.names = FALSE)
keep54 <- de_tab$in_counts
say("DUX4 candidates (mk.rda): %d features, %d in the counts, TRIM43CP present: %s; not in the counts: %s",
    nrow(de_tab), sum(de_tab$in_counts), any(de_tab$is_trim43cp),
    paste(de_tab$gene_id[!de_tab$in_counts], collapse = " "))
say("DE (padj < 0.05, |log2FC| >= 2) among the %d candidates in the counts: group 3 %d, group 4 %d (paper 47/54, 54/54); not assessed: %s",
    sum(keep54), sum(de_tab$de_S8_G3[keep54]), sum(de_tab$de_S10_G4[keep54]),
    paste(de_tab$gene_name[!keep54], collapse = " "))

writeLines(lines, file.path(out, "summary.txt"))
writeLines(capture.output(sessionInfo()), file.path(out, "r_session_info.txt"))
cat(lines, sep = "\n")
