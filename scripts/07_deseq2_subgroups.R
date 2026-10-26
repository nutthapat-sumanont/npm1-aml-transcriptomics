# 07_deseq2_subgroups.R
# Sensitivity DESeq2 (extra model only; the main model in 04_deseq2.R is
# unchanged). Splits NPM1-mutated AML into the two Session 4 groups:
#
#   group = wt (reference) | NPM1_high (primitive-like) | NPM1_low (committed-like)
#   design: ~ wave + specimen_type + sex + group
#
#  (a) Are HOXA9, HOXA10, MEIS1 and PBX3 up, and CD34 down, in BOTH subgroups
#      vs wild-type?
#  (b) Scatter of log2FC high-vs-wt against low-vs-wt, HOX genes labelled.
#  (c) Of the 1,019 main-model genes (padj < 0.05, |log2FC| > 1), how many are
#      shared by both subgroups and how many are driven mostly by one?
#
# Run from the repo root:  Rscript scripts/07_deseq2_subgroups.R  (~10 minutes)
# Needs data/ (01), results/cohort.csv (02), results/deseq2_npm1.csv (04) and
# results/block_nek10_fto_scores.csv (05).

suppressPackageStartupMessages({
  library(DESeq2)
  library(BiocParallel)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
})

set.seed(2026)
register(MulticoreParam(max(1, min(4, parallel::detectCores()))))
dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

counts_file <- "data/beataml_waves1to4_counts_dbgap.txt"
padj_cut    <- 0.05
lfc_cut     <- 1
ratio_cut   <- 2      # (c): one subgroup's |log2FC| at least 2x the other's
controls    <- c("HOXA9", "HOXA10", "MEIS1", "PBX3", "CD34")

save_fig <- function(p, file, w = 7, h = 5) {
  ggsave(file.path("figures", file), p, width = w, height = h, dpi = 300)
  cat("Saved figures/", file, "\n", sep = "")
}

# ---- Cohort and design -------------------------------------------------------
# Same factor coding as 04_deseq2.R. The high/low split is the Session 4
# k-means split of the NEK10/FTO/USP44 block score (39 high, 61 low).

block <- read.csv("results/block_nek10_fto_scores.csv")
cohort <- read.csv("results/cohort.csv") %>%
  filter(main_cohort) %>%
  left_join(block[, c("dbgap_rnaseq_sample", "block_group")],
            by = "dbgap_rnaseq_sample") %>%
  mutate(
    group = factor(case_when(NPM1 == "wt" ~ "wt",
                             block_group == "high" ~ "NPM1_high",
                             block_group == "low"  ~ "NPM1_low"),
                   levels = c("wt", "NPM1_high", "NPM1_low")),
    specimen_type = factor(ifelse(specimenType == "Bone Marrow Aspirate",
                                  "BM", "Blood"), levels = c("BM", "Blood")),
    wave = factor(cohort, levels = c("Waves1+2", "Waves3+4", "Both")),
    sex  = factor(consensus_sex, levels = c("Female", "Male")))
rownames(cohort) <- cohort$dbgap_rnaseq_sample
stopifnot(nrow(cohort) == 276, !anyNA(cohort$group))
cat("Groups:\n"); print(table(cohort$group))

raw <- fread(counts_file, check.names = FALSE)
counts <- as.matrix(raw[, cohort$dbgap_rnaseq_sample, with = FALSE])
rownames(counts) <- raw$stable_id
storage.mode(counts) <- "integer"
sym <- setNames(raw$display_label, raw$stable_id)

design_sub <- ~ wave + specimen_type + sex + group
dds <- DESeqDataSetFromMatrix(counts, colData = cohort, design = design_sub)
# Same low-count filter as the main model (>= 10 reads in >= 100 samples), so
# the gene set is identical and the 1,019 genes can be compared one to one.
dds <- dds[rowSums(counts(dds) >= 10) >= 100, ]
cat("Genes kept:", nrow(dds), "\n")

cat("Fitting", deparse(design_sub), "...\n")
dds <- DESeq(dds, parallel = TRUE, quiet = TRUE)
cat("Coefficients:", paste(resultsNames(dds), collapse = ", "), "\n")

one_contrast <- function(coef) {
  res <- results(dds, name = coef, alpha = padj_cut)
  shr <- lfcShrink(dds, coef = coef, type = "apeglm", res = res,
                   parallel = TRUE, quiet = TRUE)
  data.frame(log2FC = shr$log2FoldChange, lfcSE = shr$lfcSE,
             stat = res$stat, padj = res$padj)
}
hi <- one_contrast("group_NPM1_high_vs_wt")
lo <- one_contrast("group_NPM1_low_vs_wt")

main <- read.csv("results/deseq2_npm1.csv")
tab <- data.frame(gene_id = rownames(dds), symbol = sym[rownames(dds)],
                  log2FC_high_vs_wt = hi$log2FC, padj_high_vs_wt = hi$padj,
                  stat_high_vs_wt = hi$stat,
                  log2FC_low_vs_wt = lo$log2FC, padj_low_vs_wt = lo$padj,
                  stat_low_vs_wt = lo$stat)
tab <- tab %>%
  left_join(main %>% select(gene_id, log2FC_main = log2FC, padj_main = padj,
                            main_direction = direction), by = "gene_id")
stopifnot(nrow(tab) == nrow(main))

# ---- (c) Shared vs subgroup-driven, among the 1,019 main-model genes ---------
# Two views, because they answer slightly different questions:
#  1. Significance: is the gene padj < 0.05 in the main-model direction in each
#     contrast? Note the low group (n = 61) has more power than the high
#     group (n = 39), so this view favours "low only".
#  2. Fold change (doesn't depend on group size): shared if neither
#     subgroup's |log2FC| is more than 2x the other's (and both go the
#     main-model way); otherwise driven mostly by the larger one.

sig_dir <- function(padj, lfc, dir) {
  !is.na(padj) & padj < padj_cut & sign(lfc) == ifelse(dir == "up", 1, -1)
}
tab$main_sig <- tab$main_direction %in% c("up", "down")
tab$sig_high <- sig_dir(tab$padj_high_vs_wt, tab$log2FC_high_vs_wt, tab$main_direction)
tab$sig_low  <- sig_dir(tab$padj_low_vs_wt,  tab$log2FC_low_vs_wt,  tab$main_direction)
tab$class_by_padj <- ifelse(!tab$main_sig, NA,
  ifelse(tab$sig_high & tab$sig_low, "both",
  ifelse(tab$sig_high, "high only", ifelse(tab$sig_low, "low only", "neither"))))

s <- ifelse(tab$main_direction == "up", 1, -1)
h <- s * tab$log2FC_high_vs_wt     # effect in the main-model direction
l <- s * tab$log2FC_low_vs_wt
tab$class_by_lfc <- ifelse(!tab$main_sig, NA,
  ifelse(h <= 0 & l <= 0, "neither",
  ifelse(h >= ratio_cut * pmax(l, 0), "mostly high",
  ifelse(l >= ratio_cut * pmax(h, 0), "mostly low", "shared"))))

write.csv(tab %>% arrange(desc(abs(stat_high_vs_wt) + abs(stat_low_vs_wt))),
          "results/deseq2_subgroups.csv", row.names = FALSE)
cat("Saved results/deseq2_subgroups.csv (", nrow(tab), " genes)\n", sep = "")

m <- tab[tab$main_sig, ]
cat("\nMain-model genes:", nrow(m), "\n")
count_tab <- rbind(
  data.frame(view = "significance", direction = m$main_direction, class = m$class_by_padj),
  data.frame(view = "fold change",  direction = m$main_direction, class = m$class_by_lfc)) %>%
  count(view, class, direction) %>%
  tidyr::pivot_wider(names_from = direction, values_from = n, values_fill = 0) %>%
  mutate(total = up + down)
print(as.data.frame(count_tab))
write.csv(count_tab, "results/deseq2_subgroups_shared_counts.csv", row.names = FALSE)
cat("\nSignificance view vs fold-change view:\n")
print(table(significance = m$class_by_padj, fold_change = m$class_by_lfc))

# ---- (a) Positive controls in each subgroup ----------------------------------
ctrl <- tab %>% filter(symbol %in% controls) %>%
  mutate(expected = ifelse(symbol == "CD34", "down", "up"),
         ok_high = sig_dir(padj_high_vs_wt, log2FC_high_vs_wt, expected),
         ok_low  = sig_dir(padj_low_vs_wt,  log2FC_low_vs_wt,  expected)) %>%
  select(symbol, expected, log2FC_high_vs_wt, padj_high_vs_wt, ok_high,
         log2FC_low_vs_wt, padj_low_vs_wt, ok_low, log2FC_main) %>%
  arrange(match(symbol, controls))
write.csv(ctrl, "results/deseq2_subgroups_controls.csv", row.names = FALSE)
cat("\nPositive controls, each subgroup vs wild-type:\n")
print(ctrl, row.names = FALSE, digits = 3)

# ---- (b) Scatter -------------------------------------------------------------
hox <- grepl("^HOX[AB][0-9]+$", tab$symbol) | tab$symbol %in% c("MEIS1", "PBX3")
tab$label <- ifelse((hox & tab$main_sig) | tab$symbol %in% c(controls, "NEK10", "FTO", "USP44"),
                    tab$symbol, NA)
class_cols <- c(shared = "#0072B2", "mostly high" = "#D55E00",
                "mostly low" = "#009E73", neither = "grey30")
plot_tab <- tab %>% mutate(class = factor(class_by_lfc, levels = names(class_cols)))
lim <- range(c(tab$log2FC_high_vs_wt, tab$log2FC_low_vs_wt), na.rm = TRUE)
n_lab <- table(factor(m$class_by_lfc, levels = names(class_cols)))
p <- ggplot(plot_tab, aes(log2FC_low_vs_wt, log2FC_high_vs_wt)) +
  geom_hline(yintercept = 0, colour = "grey70") +
  geom_vline(xintercept = 0, colour = "grey70") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(data = filter(plot_tab, !main_sig), colour = "grey85", size = 0.4) +
  geom_point(data = filter(plot_tab, main_sig), aes(colour = class),
             size = 0.8, alpha = 0.8) +
  geom_text_repel(aes(label = label), size = 2.8, max.overlaps = Inf,
                  min.segment.length = 0, segment.colour = "grey40",
                  na.rm = TRUE) +
  scale_colour_manual(values = class_cols, drop = FALSE,
                      labels = sprintf("%s (%d)", names(class_cols), n_lab),
                      name = sprintf("%d main-model genes", nrow(m))) +
  coord_equal(xlim = lim, ylim = lim) +
  labs(title = "Which NPM1 genes are shared by both subgroups?",
       subtitle = paste0("Each subgroup vs NPM1 wild-type (apeglm log2FC); grey = not in the main-model list\n",
                         "'Mostly' = one subgroup's fold change is at least 2x the other's"),
       x = "log2FC, NPM1-low (committed-like, n = 61) vs wild-type",
       y = "log2FC, NPM1-high (primitive-like, n = 39) vs wild-type") +
  theme_bw(base_size = 10)
save_fig(p, "06_subgroup_shared_vs_specific.png", w = 8, h = 7.5)

writeLines(c(paste("Design:", deparse(design_sub)), capture.output(sessionInfo())),
           "results/sessionInfo_07_deseq2_subgroups.txt")
cat("\nDone.\n")
