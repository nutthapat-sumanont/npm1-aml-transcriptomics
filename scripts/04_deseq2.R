# 04_deseq2.R
# Differential expression, NPM1-mutated vs NPM1-wild-type AML, main cohort
# (276 patients after dropping BA2294R in 02_cohort.R).
#
#  1. Build the DESeq2 dataset with the approved design
#     ~ wave + specimen_type + sex + NPM1_status (wt = reference), fit it, and
#     shrink the NPM1 log2 fold changes with apeglm.
#  2. Save the full results table to results/deseq2_npm1.csv.
#  3. Count genes up and down at padj < 0.05 and |log2FC| > 1.
#  4. Volcano plot, MA plot, heatmap of the top 50 genes.
#  5. Positive-control check (HOX genes, MEIS1, PBX3, CD34) and a check that
#     sex genes (XIST, RPS4Y1, ...) are NOT among the top hits.
#  6. Sensitivity check: add FLT3-ITD to the design and see how many of the
#     top 50 genes stay.
#
# Run from the repo root:  Rscript scripts/04_deseq2.R   (about 10 minutes)
# Needs data/ from scripts/01_download.R and results/cohort.csv from 02_cohort.R.

suppressPackageStartupMessages({
  library(DESeq2)
  library(BiocParallel)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
})

set.seed(2026)
register(MulticoreParam(max(1, min(4, parallel::detectCores()))))
dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

counts_file <- "data/beataml_waves1to4_counts_dbgap.txt"
npm1_cols   <- c(wt = "grey55", mut = "#D55E00")
padj_cut    <- 0.05
lfc_cut     <- 1
coef_name   <- "NPM1_status_mut_vs_wt"

save_fig <- function(p, file, w = 7, h = 5) {
  ggsave(file.path("figures", file), p, width = w, height = h, dpi = 300)
  cat("Saved figures/", file, "\n", sep = "")
}

# ---- 1. DESeq2 dataset and fit -----------------------------------------------
# Same factor coding as 03_qc.R. Leukapheresis counts as Blood; wave keeps its
# three levels (Waves1+2, Waves3+4, Both).

cohort <- read.csv("results/cohort.csv") %>%
  filter(main_cohort) %>%
  mutate(
    NPM1_status   = factor(NPM1, levels = c("wt", "mut")),   # wt = reference
    specimen_type = factor(ifelse(specimenType == "Bone Marrow Aspirate",
                                  "BM", "Blood"), levels = c("BM", "Blood")),
    wave          = factor(cohort, levels = c("Waves1+2", "Waves3+4", "Both")),
    sex           = factor(consensus_sex, levels = c("Female", "Male")),
    FLT3_ITD      = factor(FLT3_ITD, levels = c("negative", "positive"))
  )
rownames(cohort) <- cohort$dbgap_rnaseq_sample
if (nrow(cohort) != 276 || "BA2294R" %in% rownames(cohort)) {
  stop("Expected 276 main-cohort samples without BA2294R; re-run 02_cohort.R")
}
cat("Main cohort:", nrow(cohort), "samples\n")
print(table(NPM1 = cohort$NPM1_status, wave = cohort$wave))

raw <- fread(counts_file, check.names = FALSE)
gene_info <- as.data.frame(raw[, .(stable_id, display_label, biotype)])
counts <- as.matrix(raw[, cohort$dbgap_rnaseq_sample, with = FALSE])
rownames(counts) <- gene_info$stable_id
if (!identical(colnames(counts), rownames(cohort))) {
  stop("Sample order in the count matrix does not match the cohort table")
}
if (any(counts < 0) || any(counts != round(counts))) {
  stop("Counts must be non-negative whole numbers")
}
storage.mode(counts) <- "integer"

design_main <- ~ wave + specimen_type + sex + NPM1_status
dds <- DESeqDataSetFromMatrix(counts, colData = cohort, design = design_main)
rowData(dds)$symbol  <- gene_info$display_label
rowData(dds)$biotype <- gene_info$biotype

# Low-count filter, same rule as 03_qc.R: >= 10 reads in >= 100 samples
# (100 = size of the smaller group)
keep <- rowSums(counts(dds) >= 10) >= min(table(cohort$NPM1_status))
dds  <- dds[keep, ]
cat("Genes kept after low-count filter:", nrow(dds), "\n")

# DESeq() does three things: size factors (sequencing depth), a dispersion
# (gene-level noise) per gene shrunk towards a trend, then a negative binomial
# GLM per gene with a Wald test for each coefficient.
cat("\nFitting", deparse(design_main), "...\n")
dds <- DESeq(dds, parallel = TRUE, quiet = TRUE)
cat("Coefficients:", paste(resultsNames(dds), collapse = ", "), "\n")

# Unshrunken results: Wald statistic, p-value, BH-adjusted p-value.
# alpha = 0.05 tells independent filtering which FDR we plan to use.
res <- results(dds, name = coef_name, alpha = padj_cut)
# apeglm shrinks noisy fold changes (low counts, high dispersion) towards 0 but
# leaves well-measured large effects almost untouched. Better for ranking and
# plotting. It does not change the p-values.
shr <- lfcShrink(dds, coef = coef_name, type = "apeglm", res = res,
                 parallel = TRUE, quiet = TRUE)

# ---- 2. Results table --------------------------------------------------------
# Sorted by |Wald stat|, which orders genes like the p-value but does not run
# into ties when very small p-values round to 0.

tab <- data.frame(
  gene_id         = rownames(res),
  symbol          = rowData(dds)$symbol,
  biotype         = rowData(dds)$biotype,
  baseMean        = res$baseMean,
  log2FC          = shr$log2FoldChange,        # apeglm-shrunken
  lfcSE           = shr$lfcSE,
  log2FC_unshrunk = res$log2FoldChange,
  stat            = res$stat,                  # Wald statistic, unshrunken
  pvalue          = res$pvalue,
  padj            = res$padj
) %>% arrange(desc(abs(stat)))
tab$rank <- seq_len(nrow(tab))
tab$sig  <- !is.na(tab$padj) & tab$padj < padj_cut & abs(tab$log2FC) > lfc_cut
tab$direction <- ifelse(!tab$sig, "ns", ifelse(tab$log2FC > 0, "up", "down"))
write.csv(tab, "results/deseq2_npm1.csv", row.names = FALSE)
cat("Saved results/deseq2_npm1.csv (", nrow(tab), " genes)\n", sep = "")

# ---- 3. How many genes change ------------------------------------------------

n_up   <- sum(tab$direction == "up")
n_down <- sum(tab$direction == "down")
n_fdr  <- sum(tab$padj < padj_cut, na.rm = TRUE)
cat("\nGenes with padj <", padj_cut, "(any fold change):", n_fdr, "\n")
cat("Genes with padj <", padj_cut, "and |shrunken log2FC| >", lfc_cut, ":",
    n_up + n_down, "| up in NPM1-mut:", n_up, "| down:", n_down, "\n")
cat("Genes with padj = NA (independent filtering or Cook's outliers):",
    sum(is.na(tab$padj)), "\n")
summ <- data.frame(
  genes_tested = nrow(tab), padj_lt_0.05 = n_fdr,
  up = n_up, down = n_down,
  up_protein_coding = sum(tab$direction == "up" & tab$biotype == "protein_coding"),
  down_protein_coding = sum(tab$direction == "down" & tab$biotype == "protein_coding"))
write.csv(summ, "results/deseq2_npm1_summary.csv", row.names = FALSE)

# ---- 4. Figures --------------------------------------------------------------

controls <- c("HOXA9", "HOXA10", "MEIS1", "PBX3", "CD34")
top10    <- head(tab$symbol, 10)
# -log10(padj), with padj floored at 1e-300 so genes whose p rounds to 0 still plot
tab$neglog10_padj <- -log10(pmax(tab$padj, 1e-300))
tab$label <- ifelse(tab$symbol %in% union(controls, top10), tab$symbol, NA)
dir_cols <- c(up = "#D55E00", down = "#0072B2", ns = "grey75")

p <- ggplot(filter(tab, !is.na(padj)), aes(log2FC, neglog10_padj, colour = direction)) +
  geom_point(size = 0.6, alpha = 0.6) +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = -log10(padj_cut), linetype = "dashed", colour = "grey40") +
  geom_label_repel(aes(label = label), size = 2.8, colour = "black",
                   fill = alpha("white", 0.8), max.overlaps = Inf,
                   min.segment.length = 0, box.padding = 0.4, na.rm = TRUE) +
  scale_colour_manual(values = dir_cols,
                      labels = c(up = paste0("Up (", n_up, ")"),
                                 down = paste0("Down (", n_down, ")"), ns = "Not sig.")) +
  labs(x = "log2 fold change, NPM1-mut vs wt (apeglm)", y = "-log10 adjusted p-value",
       colour = NULL, title = "NPM1-mutated vs wild-type AML (n = 276)",
       subtitle = "~ wave + specimen_type + sex + NPM1_status; padj < 0.05 and |log2FC| > 1") +
  theme_bw()
save_fig(p, "04_volcano.png", w = 8, h = 6.5)

p <- ggplot(tab, aes(log10(baseMean), log2FC, colour = direction)) +
  geom_point(size = 0.6, alpha = 0.6) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", colour = "grey40") +
  geom_text_repel(data = filter(tab, symbol %in% controls), aes(label = symbol),
                  size = 3, colour = "black", min.segment.length = 0) +
  scale_colour_manual(values = dir_cols, guide = "none") +
  labs(x = "log10 mean normalized count", y = "log2 fold change (apeglm)",
       title = "MA plot, NPM1-mutated vs wild-type",
       subtitle = "Coloured: padj < 0.05 and |log2FC| > 1; labelled: positive controls") +
  theme_bw()
save_fig(p, "04_ma_plot.png", w = 7, h = 5)

# Heatmap: VST with blind = FALSE (uses the design when estimating the
# variance trend, the usual choice for plots after the fit), then each gene
# centred and scaled across samples (z-score) so colours compare samples.
vsd <- vst(dds, blind = FALSE)
top50 <- head(tab, 50)
mat <- assay(vsd)[top50$gene_id, ]
mat <- t(scale(t(mat)))
mat[mat > 3] <- 3; mat[mat < -3] <- -3
rownames(mat) <- make.unique(top50$symbol)
ann <- data.frame(NPM1 = cohort$NPM1_status, FLT3_ITD = cohort$FLT3_ITD,
                  sex = cohort$sex, row.names = rownames(cohort))
ann_cols <- list(NPM1 = npm1_cols,
                 FLT3_ITD = c(negative = "grey90", positive = "#009E73"),
                 sex = c(Female = "#CC79A7", Male = "#56B4E9"))
ph <- pheatmap(mat, annotation_col = ann, annotation_colors = ann_cols,
               show_colnames = FALSE, fontsize_row = 7, treeheight_row = 20,
               color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(101),
               breaks = seq(-3, 3, length.out = 102),
               main = "Top 50 genes, NPM1-mut vs wt (VST z-score, capped at +/-3)",
               silent = TRUE)
ggsave("figures/04_heatmap_top50.png", ph$gtable, width = 10.8, height = 9, dpi = 300,
       bg = "white")
cat("Saved figures/04_heatmap_top50.png\n")

# ---- 5. Positive controls and sex-gene check ---------------------------------

hoxb <- sort(grep("^HOXB[0-9]+$", tab$symbol, value = TRUE))
pc_genes <- c("HOXA9", "HOXA10", hoxb, "MEIS1", "PBX3", "CD34")
expected <- setNames(ifelse(pc_genes == "CD34", "down", "up"), pc_genes)
pc <- tab %>% filter(symbol %in% pc_genes) %>%
  mutate(expected = expected[symbol],
         as_expected = ifelse(expected == "up", log2FC > 0, log2FC < 0) &
                       !is.na(padj) & padj < padj_cut) %>%
  select(symbol, gene_id, baseMean, log2FC, padj, rank, direction, expected, as_expected) %>%
  arrange(match(symbol, pc_genes))
cat("\nPositive controls (HOXB genes that passed the low-count filter only):\n")
print(pc %>% mutate(baseMean = round(baseMean), log2FC = round(log2FC, 2),
                    padj = signif(padj, 3)), row.names = FALSE)
cat("HOXB genes removed by the low-count filter:",
    paste(setdiff(paste0("HOXB", c(1:9, 13)), hoxb), collapse = ", "), "\n")
write.csv(pc, "results/deseq2_positive_controls.csv", row.names = FALSE)
if (!all(pc$as_expected[pc$symbol %in% controls])) {
  warning("A core positive control does not go the expected way; look for a bug")
}

# Sex genes. With sex in the design these should NOT be top hits. XIST and
# RPS4Y1 were asked for; the other Y genes are a wider net.
sex_genes <- c("XIST", "RPS4Y1", "KDM5D", "DDX3Y", "UTY", "EIF1AY", "ZFY", "USP9Y")
sx <- tab %>% filter(symbol %in% sex_genes) %>%
  select(symbol, log2FC, padj, rank, direction)
cat("\nSex genes (rank out of", nrow(tab), "; 1 = most significant):\n")
print(sx %>% mutate(log2FC = round(log2FC, 2), padj = signif(padj, 3)), row.names = FALSE)
cat("Any sex gene in the top 50:", any(sx$rank <= 50), "\n")
write.csv(sx, "results/deseq2_sex_genes.csv", row.names = FALSE)

# ---- 6. Sensitivity: add FLT3-ITD --------------------------------------------
# FLT3-ITD is more common in NPM1-mutated AML (41% vs 19% here). If part of
# the NPM1 signal is really FLT3-ITD, adjusting for it would change the top list.

design_flt3 <- ~ wave + specimen_type + sex + FLT3_ITD + NPM1_status
cat("\nFLT3-ITD by NPM1 status:\n"); print(table(cohort$FLT3_ITD, cohort$NPM1_status))
dds_f <- dds
design(dds_f) <- design_flt3
cat("Fitting", deparse(design_flt3), "...\n")
dds_f <- DESeq(dds_f, parallel = TRUE, quiet = TRUE)
res_f <- results(dds_f, name = coef_name, alpha = padj_cut)
shr_f <- lfcShrink(dds_f, coef = coef_name, type = "apeglm", res = res_f,
                   parallel = TRUE, quiet = TRUE)
tab_f <- data.frame(gene_id = rownames(res_f), symbol = rowData(dds_f)$symbol,
                    log2FC = shr_f$log2FoldChange, stat = res_f$stat,
                    padj = res_f$padj) %>%
  arrange(desc(abs(stat))) %>% mutate(rank = row_number())
tab_f$sig <- !is.na(tab_f$padj) & tab_f$padj < padj_cut & abs(tab_f$log2FC) > lfc_cut

top50_f   <- head(tab_f$gene_id, 50)
kept      <- intersect(top50$gene_id, top50_f)
sig_main  <- tab$gene_id[tab$sig]
sig_f     <- tab_f$gene_id[tab_f$sig]
both      <- inner_join(tab[, c("gene_id", "log2FC")], tab_f[, c("gene_id", "log2FC")],
                        by = "gene_id", suffix = c("", "_flt3"))
cat("\nTop 50 genes (main model) still in the top 50 with FLT3-ITD added:",
    length(kept), "of 50\n")
cat("Top 50 genes (main model) still significant with FLT3-ITD added:",
    sum(top50$gene_id %in% sig_f), "of 50\n")
cat("Significant genes: main", length(sig_main), "| with FLT3-ITD", length(sig_f),
    "| shared", length(intersect(sig_main, sig_f)), "\n")
cat("Correlation of shrunken log2FC, all genes:",
    round(cor(both$log2FC, both$log2FC_flt3), 3), "\n")
cat("Top 50 genes that drop out of the top 50:",
    paste(top50$symbol[!top50$gene_id %in% kept], collapse = ", "), "\n")

sens <- top50 %>% select(gene_id, symbol, rank, log2FC, padj) %>%
  left_join(tab_f %>% select(gene_id, rank_flt3 = rank, log2FC_flt3 = log2FC,
                             padj_flt3 = padj, sig_flt3 = sig), by = "gene_id") %>%
  mutate(in_top50_flt3 = rank_flt3 <= 50)
write.csv(sens, "results/deseq2_sensitivity_flt3_top50.csv", row.names = FALSE)
write.csv(data.frame(
  top50_kept_in_top50 = length(kept), top50_still_sig = sum(top50$gene_id %in% sig_f),
  sig_main = length(sig_main), sig_flt3 = length(sig_f),
  sig_shared = length(intersect(sig_main, sig_f)),
  lfc_correlation = cor(both$log2FC, both$log2FC_flt3)),
  "results/deseq2_sensitivity_flt3_summary.csv", row.names = FALSE)

# ---- Save --------------------------------------------------------------------
# Fitted objects go to data/ (gitignored); re-create them by re-running this script
saveRDS(dds, "data/dds_deseq2_main.rds")
saveRDS(vsd, "data/vsd_deseq2_main.rds")
writeLines(capture.output(sessionInfo()), "results/sessionInfo_04_deseq2.txt")
cat("\nDone.\n")
