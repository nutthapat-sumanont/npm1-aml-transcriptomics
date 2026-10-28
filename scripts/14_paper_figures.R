# 14_paper_figures.R
# Session 14 (optional): four multi-panel figures and one supplementary workbook
# for a short paper, built only from tables already in results/.
#
# No new analysis. Nothing in results/ or figures/ is rewritten; outputs go to
# figures/paper/ and results/paper/ only. Where a panel redraws a curve from
# saved values (ROC curves from saved predictions, Kaplan-Meier curves from the
# saved survival table) or a number printed on a panel comes from a table, the
# script checks that number against the value reported in LEARNING.md and stops
# if it doesn't match.
#
#  1. Figure 1, Beat AML signature: PCA, volcano, gene-set bars, HOXA9 / MEIS1.
#  2. Figure 2, robustness and replication: cell-maturity sensitivity, Beat AML
#     vs TCGA fold changes, TCGA signature score, classifier ROC.
#  3. Figure 3, NEK10/FTO/USP44-high subgroup: heat strip with mutation tracks,
#     primitive score, shared vs subgroup-specific genes, FTO locus / IFEX9.
#  4. Figure 4, drug response and survival (exploratory).
#  5. results/paper/supplementary_tables.xlsx, with a README sheet that lists
#     every table and what each column means.
#
# Style: one font (Liberation Sans, Arial metrics, if installed; otherwise the
# system sans), 7 to 8 pt text, bold panel letters, 180 mm wide (two-column
# page). Each figure is saved as PDF (cairo) and 300 dpi PNG.
#
# Run from the repo root:  Rscript scripts/14_paper_figures.R   (under a minute)
# Needs only results/ (no data/ download).

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(survival)
  library(pROC)
  library(openxlsx)
})

if (!file.exists("CLAUDE.md")) stop("Run this from the repo root")
dir.create("figures/paper", showWarnings = FALSE, recursive = TRUE)
dir.create("results/paper", showWarnings = FALSE, recursive = TRUE)
set.seed(2026)   # jitter positions only

rd <- function(f) as.data.frame(fread(file.path("results", f)))

# Stop if a number shown on a figure doesn't match the one in LEARNING.md
# (compared after rounding to the digits LEARNING.md reports).
check <- function(what, value, expected, digits = 2) {
  if (!isTRUE(all.equal(round(value, digits), expected, tolerance = 1e-9))) {
    stop(sprintf("%s: got %s, LEARNING.md says %s", what, round(value, digits), expected))
  }
  invisible(value)
}
# Stop if a value redrawn from saved data differs from the saved table value
same <- function(what, value, saved) {
  if (abs(value - saved) > 1e-8) stop(sprintf("%s: got %.10g, table says %.10g", what, value, saved))
  invisible(value)
}

# ---- Style -------------------------------------------------------------------
fams <- tryCatch(system("fc-list : family", intern = TRUE), error = function(e) character())
font <- if (any(grepl("Liberation Sans", fams, fixed = TRUE))) "Liberation Sans" else "sans"

col_npm1  <- "#D55E00"   # NPM1-mutated (npm1_cols / group_cols in 03 to 10)
col_kmt2a <- "#0072B2"   # KMT2A-rearranged (group_cols in 05, 08)
col_other <- "grey55"    # other AML / NPM1 wild-type
col_high  <- "#A63D00"   # NPM1-mutated, NEK10/FTO/USP44-high (primitive-like)
col_low   <- "#F4A582"   # NPM1-mutated, NEK10/FTO/USP44-low (committed-like)
edge_low  <- "#B5622F"   # thin darker outline so the light orange shows on white
col_up    <- "#D55E00"   # up / down genes (dir_cols in 04, 07)
col_down  <- "#0072B2"
col_ns    <- "grey80"

group3_cols <- c("NPM1-mutated" = col_npm1, "KMT2A-rearranged" = col_kmt2a, "Other AML" = col_other)
sub_fill <- c(high = col_high, low = col_low)
sub_edge <- c(high = col_high, low = edge_low)

theme_paper <- function() {
  theme_classic(base_size = 8, base_family = font) +
    theme(axis.line = element_line(linewidth = 0.3),
          axis.ticks = element_line(linewidth = 0.3),
          axis.text = element_text(size = 7, colour = "black"),
          axis.title = element_text(size = 7.5),
          legend.text = element_text(size = 6.5),
          legend.title = element_text(size = 7),
          legend.key.size = unit(3, "mm"),
          legend.background = element_blank(),
          plot.title = element_text(size = 8, face = "bold"),
          plot.subtitle = element_text(size = 6.5, colour = "grey30"),
          strip.background = element_blank(),
          strip.text = element_text(size = 7.5, face = "bold"),
          plot.tag = element_text(size = 10, face = "bold", family = font))
}
theme_set(theme_paper())
txt <- 7 / .pt     # geom_text size equal to 7 pt
txt_small <- 6 / .pt

save_paper <- function(p, name, h_mm, w_mm = 180) {
  p <- p + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 10, face = "bold"))
  ggsave(file.path("figures/paper", paste0(name, ".pdf")), p, width = w_mm, height = h_mm,
         units = "mm", device = cairo_pdf)
  ggsave(file.path("figures/paper", paste0(name, ".png")), p, width = w_mm, height = h_mm,
         units = "mm", dpi = 300, type = "cairo")
  cat("Saved figures/paper/", name, ".pdf and .png\n", sep = "")
}
fmt_p <- function(p) ifelse(p < 0.001, sprintf("%.0e", p), sprintf("%.2g", p))

# ==== Figure 1. The NPM1 signature in Beat AML ================================

# A. PCA (Session 2, 277 samples before QC; BA2294R marked)
pca <- rd("qc_pca_scores.csv")
pca_r2 <- rd("qc_pca_r2.csv")
check("PC2 R2 with NPM1", pca_r2$NPM1_status[pca_r2$V1 == "PC2"], 0.43)
stopifnot(nrow(pca) == 277, "BA2294R" %in% pca$sample)
pca$grp <- ifelse(pca$sample == "BA2294R", "excluded",
                  ifelse(pca$NPM1_status == "mut", "mut", "wt"))
n_pca <- table(pca$grp)
check("PCA NPM1-mutated", n_pca[["mut"]], 100, 0); check("PCA wild-type", n_pca[["wt"]], 176, 0)
pca_lab <- c(mut = sprintf("NPM1-mutated (%d)", n_pca[["mut"]]),
             wt = sprintf("NPM1 wild-type (%d)", n_pca[["wt"]]))
# % variance isn't in results/; these are the Session 2 values (LEARNING.md)
p1a <- ggplot(filter(pca, grp != "excluded"), aes(PC1, PC2)) +
  geom_point(aes(fill = grp), shape = 21, size = 1.3, stroke = 0.15, colour = "white") +
  geom_point(data = filter(pca, grp == "excluded"), shape = 21, size = 2, stroke = 0.6,
             fill = NA, colour = "black") +
  geom_text_repel(data = filter(pca, grp == "excluded"),
                  label = "BA2294R\nexcluded after QC:\nsex mismatch", size = txt_small,
                  family = font, lineheight = 0.9, nudge_x = 18, nudge_y = -4,
                  segment.size = 0.25, min.segment.length = 0) +
  scale_fill_manual(values = c(mut = col_npm1, wt = col_other), labels = pca_lab, name = NULL) +
  labs(x = "PC1 (19.7% of variance)", y = "PC2 (13.0% of variance)",
       title = "PCA, top 500 variable genes (n = 277)",
       subtitle = "NPM1 status explains 43% of PC2") +
  theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, 0))

# B. Volcano (Session 3)
de <- rd("deseq2_npm1.csv")
n_up <- sum(de$sig & de$direction == "up"); n_down <- sum(de$sig & de$direction == "down")
check("Significant up", n_up, 307, 0); check("Significant down", n_down, 712, 0)
de$y <- -log10(pmax(de$padj, 1e-300))
de$dir <- factor(ifelse(de$sig, de$direction, "ns"), levels = c("up", "down", "ns"))
vol_genes <- c("HOXA9", "HOXA10", "MEIS1", "PBX3", "CD34", "NEK10", "FTO")
de$label <- ifelse(de$symbol %in% vol_genes,
                   ifelse(de$symbol == "FTO", "FTO locus\n(IFEX9)", de$symbol), NA)
p1b <- ggplot(filter(de, !is.na(padj)), aes(log2FC, y)) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "grey50", linewidth = 0.25) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey50", linewidth = 0.25) +
  geom_point(data = filter(de, !is.na(padj), dir == "ns"), colour = col_ns, size = 0.35) +
  geom_point(data = filter(de, dir != "ns"), aes(colour = dir), size = 0.45) +
  geom_point(data = filter(de, !is.na(label)), shape = 21, size = 1.4, stroke = 0.4,
             colour = "black", fill = NA) +
  geom_text_repel(aes(label = label), size = txt_small, family = font, lineheight = 0.85,
                  min.segment.length = 0, segment.size = 0.25, box.padding = 0.35,
                  max.overlaps = Inf, na.rm = TRUE) +
  scale_colour_manual(values = c(up = col_up, down = col_down), guide = "none") +
  labs(x = "log2 fold change, NPM1-mutated vs wild-type",
       y = "-log10 adjusted p",
       title = "Differential expression (n = 276)",
       subtitle = sprintf("%d up, %d down at padj < 0.05 and 2-fold (dashed lines)", n_up, n_down))

# C. Gene-set enrichment (Session 4): Verhaak UP / DN, then the top 5 up and
# top 5 down main collapsed GO BP terms by NES (rule fixed in the figure plan)
ver <- rd("gsea_verhaak_controls.csv")
check("Verhaak UP NES", ver$NES[ver$pathway == "VERHAAK_AML_WITH_NPM1_MUTATED_UP"], 3.31)
check("Verhaak DN NES", ver$NES[ver$pathway == "VERHAAK_AML_WITH_NPM1_MUTATED_DN"], -3.76)
go <- rd("gsea_gobp_collapsed.csv") %>% filter(main_pathway, padj < 0.05)
go_sel <- bind_rows(go %>% arrange(desc(NES)) %>% head(5), go %>% arrange(NES) %>% head(5))
clean_go <- function(x) {
  x <- tolower(gsub("_", " ", sub("^GOBP_", "", x)))
  x <- paste0(toupper(substring(x, 1, 1)), substring(x, 2))
  x <- gsub("\\bmhc\\b", "MHC", x); gsub("\\bii\\b", "II", x)
}
gs <- bind_rows(
  ver %>% transmute(label = ifelse(NES > 0, "Verhaak NPM1-mutated UP", "Verhaak NPM1-mutated DN"),
                    NES, padj, set = "Verhaak"),
  go_sel %>% transmute(label = clean_go(pathway), NES, padj, set = "GO BP")) %>%
  mutate(label = vapply(label, function(s) paste(strwrap(s, 38), collapse = "\n"), ""),
         set = factor(set, levels = c("Verhaak", "GO BP")),
         dir = ifelse(NES > 0, "up", "down"),
         fdr = paste("FDR", fmt_p(padj)))
gs$label <- factor(gs$label, levels = rev(unique(gs$label[order(gs$set, -gs$NES)])))
p1c <- ggplot(gs, aes(NES, label, fill = dir)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  geom_text(aes(x = 6.4, label = fmt_p(padj)), hjust = 1, size = txt_small, family = font) +
  geom_text(data = gs[gs$label == levels(gs$label)[nlevels(gs$label)], ], aes(x = 6.4, y = label), label = "FDR", hjust = 1, vjust = -1.4,
            size = txt_small, family = font, fontface = "bold", inherit.aes = FALSE) +
  facet_grid(rows = vars(set), scales = "free_y", space = "free_y") +
  scale_fill_manual(values = c(up = col_up, down = col_down), guide = "none") +
  scale_x_continuous(limits = c(-4.2, 6.5), breaks = seq(-4, 4, 2)) +
  coord_cartesian(clip = "off") +
  labs(x = "Normalized enrichment score (NES)", y = NULL,
       title = "Gene-set enrichment") +
  theme(axis.text.y = element_text(size = 6, lineheight = 0.85),
        strip.text.y = element_text(size = 6.5, angle = 0, hjust = 0))

# D. HOXA9 and MEIS1 in three groups (Session 4)
tg <- rd("three_group_hoxa9_meis1.csv")
tg$group3 <- factor(tg$group3, levels = names(group3_cols))
tg_med <- tg %>% group_by(gene, group3) %>% summarise(n = n(), med = median(vst), .groups = "drop")
check("HOXA9 median NPM1-mut", tg_med$med[tg_med$gene == "HOXA9" & tg_med$group3 == "NPM1-mutated"], 12.2, 1)
check("HOXA9 median KMT2A-r", tg_med$med[tg_med$gene == "HOXA9" & tg_med$group3 == "KMT2A-rearranged"], 12.3, 1)
check("HOXA9 median other", tg_med$med[tg_med$gene == "HOXA9" & tg_med$group3 == "Other AML"], 9.5, 1)
tg_n <- tg_med %>% filter(gene == "HOXA9")
tg_tests <- rd("three_group_tests.csv") %>% filter(comparison == "NPM1-mutated vs KMT2A-rearranged")
check("HOXA9 NPM1 vs KMT2A-r p", tg_tests$wilcox_p[tg_tests$gene == "HOXA9"], 0.42)
tg_ann <- tg_tests %>% transmute(gene, lab = sprintf("p = %.2f", wilcox_p), x = 1, xend = 2, y = 14.9)
tg$double <- tg$group3 == "KMT2A-rearranged" & tg$NPM1 == "mut"
p1d <- ggplot(tg, aes(group3, vst)) +
  geom_jitter(aes(fill = group3, shape = double), width = 0.22, height = 0, size = 0.8,
              stroke = 0.1, colour = "white") +
  geom_boxplot(outlier.shape = NA, fill = NA, width = 0.55, linewidth = 0.3) +
  geom_segment(data = tg_ann, aes(x = x, xend = xend, y = y, yend = y), inherit.aes = FALSE,
               linewidth = 0.25) +
  geom_text(data = tg_ann, aes(x = 1.5, y = 15.4, label = lab), inherit.aes = FALSE,
            size = txt_small, family = font) +
  facet_wrap(~ gene) +
  scale_fill_manual(values = group3_cols, guide = "none") +
  scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 24), labels = "KMT2A-r with NPM1 mutation",
                     breaks = "TRUE", name = NULL) +
  scale_x_discrete(labels = sprintf("%s\n(%d)", c("NPM1-\nmut", "KMT2A-\nr", "Other\nAML"), tg_n$n)) +
  scale_y_continuous(limits = c(NA, 15.8)) +
  guides(shape = guide_legend(override.aes = list(fill = col_kmt2a, size = 1.5))) +
  labs(x = NULL, y = "Expression (VST)", title = "HOXA9 and MEIS1 (menin targets)") +
  theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, 0),
        axis.text.x = element_text(size = 6, lineheight = 0.85))

fig1 <- (p1a | p1b) / (free(p1c) + p1d + plot_layout(widths = c(1.25, 1))) + plot_layout(heights = c(1, 1.15))
save_paper(fig1, "figure1_npm1_signature", h_mm = 165)

# ==== Figure 2. Robust to cell maturity, replicates in TCGA-LAML =============

# A. Cell-maturity sensitivity (Session 6): main design on the same 270 vs the
# two-score model; the all-six upper bound is not shown
ca <- rd("deseq2_celltype_adjusted.csv")
cs <- rd("deseq2_celltype_summary.csv")
r_ca <- cor(ca$log2FC_main_270, ca$log2FC_two_score, use = "complete.obs")
same("Two-score vs main_270 Pearson", r_ca, cs$lfc_pearson_vs_main_270[cs$model == "two_score"])
check("Two-score vs main_270 Pearson", r_ca, 0.89)
slope_ca <- cs$lfc_slope_vs_main_270[cs$model == "two_score"]
check("Two-score slope", slope_ca, 0.85)
check("Two-score top 50 still significant", cs$top50_still_sig[cs$model == "two_score"], 49, 0)
ca$dir <- ifelse(ca$sig_main, ifelse(ca$log2FC_main > 0, "up", "down"), "ns")
ca_lab <- c(HOXA9 = "HOXA9", PBX3 = "PBX3", CD34 = "CD34", MEIS1 = "MEIS1 (partly\ncircular)")
ca$label <- unname(ca_lab[ca$symbol])
p2a <- ggplot(ca, aes(log2FC_main_270, log2FC_two_score)) +
  geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.25) +
  geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.25) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.3) +
  geom_point(data = filter(ca, dir == "ns"), colour = col_ns, size = 0.3) +
  geom_point(data = filter(ca, dir != "ns"), aes(colour = dir), size = 0.4) +
  geom_text_repel(aes(label = label), size = txt_small, family = font, lineheight = 0.85,
                  min.segment.length = 0, segment.size = 0.25, box.padding = 0.4, na.rm = TRUE) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.3, size = txt, family = font,
           label = sprintf("Pearson r = %.2f\nslope (1,019 genes) = %.2f", r_ca, slope_ca)) +
  scale_colour_manual(values = c(up = col_up, down = col_down), guide = "none") +
  coord_equal() +
  labs(x = "log2FC, main design (n = 270)", y = "log2FC, + HSC-like + monocyte-like",
       title = "Adjusting for cell maturity",
       subtitle = "Coloured: the 1,019 main-model genes")

# B. Beat AML vs TCGA fold changes (Session 7), genes tested in both
tt <- rd("tcga_deseq2_npm1.csv")
cmp <- inner_join(de %>% select(gene_id, symbol, lfc_beat = log2FC, sig_beat = sig),
                  tt %>% select(gene_id, lfc_tcga = log2FC), by = "gene_id")
tvb <- rd("tcga_vs_beataml_summary.csv")
check("Genes tested in both", nrow(cmp), 17498, 0)
r_all <- cor(cmp$lfc_beat, cmp$lfc_tcga)
r_sig <- cor(cmp$lfc_beat[cmp$sig_beat], cmp$lfc_tcga[cmp$sig_beat])
check("Beat vs TCGA r, all", r_all, 0.74); check("Beat vs TCGA r, significant", r_sig, 0.84)
check("Same direction %", tvb$value[tvb$metric == "Beat AML significant: same direction in TCGA (%)"], 95.4, 1)
cmp$dir <- ifelse(cmp$sig_beat, ifelse(cmp$lfc_beat > 0, "up", "down"), "ns")
cmp$label <- ifelse(cmp$symbol %in% c("HOXA9", "MEIS1", "PBX3", "CD34", "NEK10"), cmp$symbol, NA)
p2b <- ggplot(cmp, aes(lfc_beat, lfc_tcga)) +
  geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.25) +
  geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.25) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.3) +
  geom_point(data = filter(cmp, dir == "ns"), colour = col_ns, size = 0.3) +
  geom_point(data = filter(cmp, dir != "ns"), aes(colour = dir), size = 0.4) +
  geom_text_repel(aes(label = label), size = txt_small, family = font, min.segment.length = 0,
                  segment.size = 0.25, box.padding = 0.4, na.rm = TRUE) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.3, size = txt, family = font,
           label = sprintf("r = %.2f (all %s genes)\nr = %.2f (Beat AML significant)",
                           r_all, format(nrow(cmp), big.mark = ","), r_sig)) +
  scale_colour_manual(values = c(up = col_up, down = col_down), guide = "none") +
  coord_equal() +
  labs(x = "log2FC, Beat AML (n = 276)", y = "log2FC, TCGA-LAML (n = 111)",
       title = "Replication in TCGA-LAML",
       subtitle = "Coloured: significant in Beat AML")

# C. TCGA signature score by group (Session 7)
ss <- rd("tcga_signature_scores.csv")
ss$group3 <- factor(ss$group3, levels = names(group3_cols))
main_ss <- filter(ss, group3 != "KMT2A-rearranged")
roc_ss <- roc(main_ss$NPM1, main_ss$score, levels = c("wt", "mut"), direction = "<", quiet = TRUE)
ci_ss <- ci.auc(roc_ss, method = "delong")
check("TCGA signature AUC", as.numeric(auc(roc_ss)), 0.992, 3)
check("TCGA signature AUC CI low", ci_ss[1], 0.981, 3)
check("TCGA patients in ROC", nrow(main_ss), 111, 0)
ss_n <- table(ss$group3)
p2c <- ggplot(ss, aes(group3, score)) +
  geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.25) +
  geom_jitter(aes(fill = group3), shape = 21, width = 0.22, height = 0, size = 1, stroke = 0.1,
              colour = "white") +
  geom_boxplot(outlier.shape = NA, fill = NA, width = 0.55, linewidth = 0.3) +
  annotate("text", x = 3.45, y = Inf, hjust = 1, vjust = 1.3, size = txt, family = font,
           label = sprintf("NPM1-mutated vs other AML\nROC AUC %.3f (95%% CI %.3f to %.3f)",
                           auc(roc_ss), ci_ss[1], min(ci_ss[3], 1))) +
  scale_fill_manual(values = group3_cols, guide = "none") +
  scale_x_discrete(labels = sprintf("%s\n(%d)", sub("-rearranged", "-r", names(ss_n)), as.integer(ss_n))) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
  labs(x = NULL, y = "Signature score (mean z, top 50 up - top 50 down)",
       title = "Beat AML signature scored in TCGA-LAML") +
  theme(axis.text.x = element_text(size = 6.5))

# D. Classifier ROC (Session 13): solid thin = Beat AML nested CV (one line per
# repeat), dashed thick = TCGA locked test (rescaled); black = LASSO, grey = HOXA9 + CD34
cvp <- rd("classifier_cv_predictions.csv")
tcp <- rd("classifier_tcga_predictions.csv") %>% filter(scaling == "rescaled")
cvs <- rd("classifier_cv_summary.csv")
tperf <- rd("classifier_tcga_performance.csv") %>% filter(scaling == "rescaled")
roc_df <- function(y, p) {
  r <- roc(y, p, levels = c(0, 1), direction = "<", quiet = TRUE)
  list(auc = as.numeric(auc(r)),
       df = data.frame(fpr = 1 - r$specificities, tpr = r$sensitivities) %>% arrange(fpr, tpr))
}
cv_curves <- list(); cv_auc <- list(lasso = c(), simple = c())
for (rr in sort(unique(cvp$repeat_id))) {
  d <- cvp[cvp$repeat_id == rr, ]
  for (m in c("lasso", "simple")) {
    r <- roc_df(d$y, d[[m]])
    cv_auc[[m]] <- c(cv_auc[[m]], r$auc)
    cv_curves[[paste(m, rr)]] <- cbind(r$df, model = m, cohort = "beat", grp = paste(m, rr))
  }
}
same("CV LASSO AUC", mean(cv_auc$lasso), cvs$auc_mean[cvs$model == "LASSO (lambda.1se)"])
same("CV HOXA9 + CD34 AUC", mean(cv_auc$simple), cvs$auc_mean[cvs$model == "HOXA9 + CD34"])
check("CV LASSO AUC", mean(cv_auc$lasso), 0.985, 3)
check("CV HOXA9 + CD34 AUC", mean(cv_auc$simple), 0.950, 3)
tc <- list()
tc_auc <- c()
for (m in c("lasso", "hoxa9_cd34")) {
  d <- tcp[tcp$model == m, ]
  r <- roc_df(as.integer(d$NPM1 == "mut"), d$p)
  tc_auc[m] <- r$auc
  tc[[m]] <- cbind(r$df, model = ifelse(m == "lasso", "lasso", "simple"), cohort = "tcga", grp = m)
}
same("TCGA LASSO AUC", tc_auc[["lasso"]], tperf$auc[tperf$model == "lasso"])
same("TCGA HOXA9 + CD34 AUC", tc_auc[["hoxa9_cd34"]], tperf$auc[tperf$model == "hoxa9_cd34"])
check("TCGA LASSO AUC", tc_auc[["lasso"]], 0.996, 3)
check("TCGA HOXA9 + CD34 AUC", tc_auc[["hoxa9_cd34"]], 0.966, 3)
curves <- bind_rows(c(cv_curves, tc)) %>%
  mutate(key = factor(paste(cohort, model),
                      levels = c("beat lasso", "beat simple", "tcga lasso", "tcga simple")))
f3 <- function(x) sprintf("%.3f", x)
lp <- tperf[tperf$model == "lasso", ]; sp <- tperf[tperf$model == "hoxa9_cd34", ]
roc_labs <- c(
  "beat lasso"  = sprintf("Beat AML CV, LASSO (33 genes): AUC %s (repeats %s to %s)",
                          f3(mean(cv_auc$lasso)), f3(min(cv_auc$lasso)), f3(max(cv_auc$lasso))),
  "beat simple" = sprintf("Beat AML CV, HOXA9 + CD34: AUC %s (repeats %s to %s)",
                          f3(mean(cv_auc$simple)), f3(min(cv_auc$simple)), f3(max(cv_auc$simple))),
  "tcga lasso"  = sprintf("TCGA locked, LASSO: AUC %s (95%% CI %s to %s)",
                          f3(lp$auc), f3(lp$auc_lo), f3(lp$auc_hi)),
  "tcga simple" = sprintf("TCGA locked, HOXA9 + CD34: AUC %s (95%% CI %s to %s)",
                          f3(sp$auc), f3(sp$auc_lo), f3(sp$auc_hi)))
roc_cols <- c("beat lasso" = "black", "beat simple" = "grey60", "tcga lasso" = "black", "tcga simple" = "grey60")
roc_lty  <- c("beat lasso" = "solid", "beat simple" = "solid", "tcga lasso" = "22", "tcga simple" = "22")
roc_lwd  <- c("beat lasso" = 0.2, "beat simple" = 0.2, "tcga lasso" = 0.75, "tcga simple" = 0.75)
roc_plot <- function(d, xmax = 1, ymin = 0, legend = TRUE) {
  ggplot(d, aes(fpr, tpr, group = grp, colour = key, linetype = key, linewidth = key)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey85", linewidth = 0.25) +
    geom_path() +
    scale_colour_manual(values = roc_cols, labels = roc_labs, name = NULL) +
    scale_linetype_manual(values = roc_lty, labels = roc_labs, name = NULL) +
    scale_linewidth_manual(values = roc_lwd, labels = roc_labs, name = NULL) +
    coord_cartesian(xlim = c(0, xmax), ylim = c(ymin, 1), expand = FALSE) +
    theme(legend.position = if (legend) "bottom" else "none")
}
p2d_main <- roc_plot(curves) +
  guides(colour = guide_legend(ncol = 1), linetype = guide_legend(ncol = 1),
         linewidth = guide_legend(ncol = 1, override.aes = list(linewidth = c(0.4, 0.4, 0.75, 0.75)))) +
  scale_x_continuous(breaks = seq(0, 1, 0.25)) +
  labs(x = "False positive rate (1 - specificity)", y = "True positive rate (sensitivity)",
       title = "Predicting NPM1 status from expression") +
  theme(aspect.ratio = 1, legend.key.width = unit(6, "mm"), legend.text = element_text(size = 6),
        legend.margin = margin(0, 0, 0, 0))
p2d_zoom <- roc_plot(curves, xmax = 0.25, ymin = 0.75, legend = FALSE) +
  scale_x_continuous(breaks = c(0, 0.1, 0.2)) + scale_y_continuous(breaks = c(0.8, 0.9, 1)) +
  labs(x = NULL, y = NULL) +
  theme(axis.text = element_text(size = 5), plot.background = element_rect(colour = "grey50", linewidth = 0.3),
        panel.border = element_rect(fill = NA, colour = NA))
# wrapped so the panel and its inset take one tag (D)
p2d <- wrap_elements(full = p2d_main + inset_element(p2d_zoom, left = 0.42, bottom = 0.04, right = 0.98,
                                                     top = 0.6, ignore_tag = TRUE))

fig2 <- (p2a | p2b) / (p2c | p2d) + plot_layout(heights = c(1, 1.25))
save_paper(fig2, "figure2_robustness_replication", h_mm = 190)

# ==== Figure 3. The NEK10/FTO/USP44-high subgroup ============================

blk <- rd("block_nek10_fto_scores.csv")
sps <- rd("subgroup_per_sample.csv")
com <- rd("subgroup_comutations.csv")
stopifnot(nrow(blk) == 100, setequal(blk$dbgap_rnaseq_sample, sps$dbgap_rnaseq_sample))
n_sub <- table(blk$block_group)
check("High group", n_sub[["high"]], 37, 0); check("Low group", n_sub[["low"]], 63, 0)
check("High primitive-like", sum(sps$block_group == "high" & sps$pc_call == "primitive-like"), 32, 0)
check("Low committed-like", sum(sps$block_group == "low" & sps$pc_call == "committed-like"), 54, 0)
ord <- blk$dbgap_rnaseq_sample[order(blk$core_score)]
sub_lab <- c(high = sprintf("NEK10/FTO/USP44-high, primitive-like (%d)", n_sub[["high"]]),
             low = sprintf("NEK10/FTO/USP44-low, committed-like (%d)", n_sub[["low"]]))

# A. Heat strip (z-score of VST across the 100) plus annotation tracks
heat <- blk %>% select(sample = dbgap_rnaseq_sample, NEK10, FTO, USP44) %>%
  pivot_longer(-sample, names_to = "gene", values_to = "vst") %>%
  group_by(gene) %>% mutate(z = (vst - mean(vst)) / sd(vst)) %>% ungroup() %>%
  mutate(sample = factor(sample, levels = ord),
         gene = factor(ifelse(gene == "FTO", "FTO locus (IFEX9)", gene),
                       levels = rev(c("NEK10", "FTO locus (IFEX9)", "USP44"))))
cut_x <- n_sub[["low"]] + 0.5
# counts, not %, so no rounding (7 of 61 = 11.5%) can disagree with the legend
pct <- function(f) { r <- com[com$feature == f, ]; sprintf("%d/%d vs %d/%d", r$high_n, r$high_of, r$low_n, r$low_of) }
check("TET2 high %", com$high_pct[com$feature == "TET2"], 33, 0)
check("IDH2 high %", com$high_pct[com$feature == "IDH2"], 39, 0)
check("PTPN11 low %", com$low_pct[com$feature == "PTPN11"], 25, 0)
tracks_def <- c(FLT3_ITD = "FLT3_ITD_pos", mut_TET2 = "TET2", mut_IDH2 = "IDH2",
                mut_DNMT3A = "DNMT3A", mut_NRAS = "NRAS", mut_PTPN11 = "PTPN11")
track_names <- c(FLT3_ITD = "FLT3-ITD", mut_TET2 = "TET2", mut_IDH2 = "IDH2",
                 mut_DNMT3A = "DNMT3A", mut_NRAS = "NRAS", mut_PTPN11 = "PTPN11")
trk <- sps %>% mutate(FLT3_ITD = FLT3_ITD == "positive") %>%
  select(sample = dbgap_rnaseq_sample, pc_call, block_group, all_of(names(tracks_def))) %>%
  mutate(across(all_of(names(tracks_def)), ~ ifelse(is.na(.x), "no exome data",
                                                     ifelse(.x, "mutated", "not mutated")))) %>%
  mutate(group = block_group, mer = pc_call) %>%
  select(-block_group, -pc_call) %>%
  pivot_longer(-sample, names_to = "track", values_to = "value")
trk_lab <- c(group = "Subgroup (k-means)", mer = "Mer et al. call",
             setNames(sprintf("%s  (%s)", track_names, vapply(tracks_def, pct, "")), names(tracks_def)))
trk <- trk %>% mutate(sample = factor(sample, levels = ord),
                      track = factor(track, levels = rev(c("group", "mer", names(tracks_def))),
                                     labels = rev(trk_lab[c("group", "mer", names(tracks_def))])))
trk_fill <- c(high = col_high, low = col_low, "primitive-like" = "grey25", "committed-like" = "grey80",
              mutated = "black", "not mutated" = "grey93", "no exome data" = "white")
p3a_heat <- ggplot(heat, aes(sample, gene, fill = z)) +
  geom_tile() +
  geom_vline(xintercept = cut_x, linewidth = 0.4, linetype = "dashed") +
  annotate("text", x = cut_x / 2, y = 3.9, label = sprintf("Low (%d)", n_sub[["low"]]), size = txt,
           family = font) +
  annotate("text", x = cut_x + n_sub[["high"]] / 2, y = 3.9, label = sprintf("High (%d)", n_sub[["high"]]),
           size = txt, family = font) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", limits = c(-2.5, 2.5),
                       oob = scales::squish, name = "z-score\n(VST)") +
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = expansion(add = c(0.5, 1))) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = NULL, title = "100 NPM1-mutated samples, ordered by the NEK10/FTO/USP44 score") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line = element_blank(),
        legend.key.height = unit(2.5, "mm"))
p3a_trk <- ggplot(trk, aes(sample, track, fill = value)) +
  geom_tile(colour = "white", linewidth = 0.1) +
  geom_vline(xintercept = cut_x, linewidth = 0.4, linetype = "dashed") +
  scale_fill_manual(values = trk_fill, breaks = names(trk_fill), name = NULL,
                    labels = c(high = "High (primitive-like)", low = "Low (committed-like)",
                               "primitive-like" = "Primitive-like", "committed-like" = "Committed-like",
                               mutated = "Mutated / positive", "not mutated" = "Not mutated",
                               "no exome data" = "No exome data")) +
  scale_x_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL, caption = "Brackets: mutated / with data, high vs low group") +
  guides(fill = guide_legend(ncol = 2, override.aes = list(colour = "grey60", linewidth = 0.2))) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line = element_blank(),
        axis.text.y = element_text(size = 6.5), legend.key.size = unit(2.5, "mm"),
        plot.caption = element_text(size = 6, colour = "grey30", hjust = 0))
p3a <- wrap_elements(full = p3a_heat / p3a_trk + plot_layout(heights = c(3, 8)))

# B. Primitive score vs the NEK10/FTO/USP44 core score (Session 5, rerun in 8)
prim <- rd("subgroup_primitive_committed_scores.csv")
rho <- prim$spearman_vs_core[prim$score == "primitive_score"]
# LEARNING.md (Sessions 5, 8) says 0.91; the committed table says 0.9048, so the
# figure shows 0.90. A rounding slip (0.9048 -> 0.905 -> 0.91), noted in Session 14.
check("Primitive vs core Spearman", rho, 0.90)
p3b <- ggplot(sps, aes(core_score, primitive_score)) +
  geom_point(aes(fill = block_group, colour = block_group, shape = pc_call), size = 1.3, stroke = 0.3) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.3, size = txt, family = font,
           label = sprintf("Spearman %.2f\n32 of 37 high primitive-like\n54 of 63 low committed-like", rho)) +
  scale_fill_manual(values = sub_fill, labels = sub_lab, name = NULL) +
  scale_colour_manual(values = sub_edge, labels = sub_lab, name = NULL) +
  scale_shape_manual(values = c("primitive-like" = 21, "committed-like" = 24), name = "Mer et al. call") +
  guides(fill = guide_legend(ncol = 1, override.aes = list(shape = 21)), colour = "none",
         shape = guide_legend(ncol = 1, override.aes = list(fill = "grey70", colour = "grey30"))) +
  labs(x = "NEK10/FTO/USP44 score (mean z)", y = "Primitive score (Mer et al. markers, mean z)",
       title = "The high group is the primitive subtype") +
  theme(legend.position = "bottom", legend.box = "vertical", legend.spacing.y = unit(0, "mm"),
        legend.margin = margin(0, 0, 0, 0))

# C. Shared vs subgroup-specific (Session 5, rerun in 8)
ds <- rd("deseq2_subgroups.csv")
cnt <- rd("deseq2_subgroups_shared_counts.csv") %>% filter(view == "fold change")
n_cls <- setNames(cnt$total, cnt$class)
check("Shared", n_cls[["shared"]], 521, 0); check("Mostly high", n_cls[["mostly high"]], 381, 0)
check("Mostly low", n_cls[["mostly low"]], 117, 0)
cls_fill <- c(shared = "grey35", "mostly high" = col_high, "mostly low" = col_low)
cls_edge <- c(shared = "grey35", "mostly high" = col_high, "mostly low" = edge_low)
cls_lab <- setNames(sprintf("%s (%d)", names(cls_fill), n_cls[names(cls_fill)]), names(cls_fill))
ds$label <- ifelse(ds$symbol %in% c("HOXA9", "HOXA10", "MEIS1", "PBX3", "CD34", "NEK10", "USP44"),
                   ds$symbol, ifelse(ds$symbol == "FTO", "FTO locus", NA))
lim <- range(c(ds$log2FC_high_vs_wt, ds$log2FC_low_vs_wt), na.rm = TRUE)
p3c <- ggplot(ds, aes(log2FC_low_vs_wt, log2FC_high_vs_wt)) +
  geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.25) +
  geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.25) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.3) +
  geom_point(data = filter(ds, !main_sig), colour = "grey88", size = 0.3) +
  geom_point(data = filter(ds, main_sig), aes(fill = class_by_lfc, colour = class_by_lfc),
             shape = 21, size = 0.8, stroke = 0.15) +
  geom_text_repel(aes(label = label), size = txt_small, family = font, min.segment.length = 0,
                  segment.size = 0.25, box.padding = 0.4, max.overlaps = Inf, na.rm = TRUE) +
  scale_fill_manual(values = cls_fill, labels = cls_lab, name = "1,019 main-model genes") +
  scale_colour_manual(values = cls_edge, labels = cls_lab, name = "1,019 main-model genes") +
  coord_equal(xlim = lim, ylim = lim) +
  guides(fill = guide_legend(override.aes = list(size = 1.6))) +
  labs(x = "log2FC, low (committed-like) vs wild-type", y = "log2FC, high (primitive-like) vs wild-type",
       title = "HOX shared, NEK10/FTO/USP44 high only",
       subtitle = "'Mostly' = one subgroup's fold change at least 2x the other's") +
  theme(legend.position = "bottom", legend.direction = "vertical", legend.margin = margin(0, 0, 0, 0))

# D. FTO locus in TCGA-LAML, per patient (Session 10)
fto <- rd("fto_locus_per_sample.csv")
fg_n <- table(fto$fto_group)
check("TCGA FTO-high", fg_n[["FTO_high"]], 6, 0); check("TCGA FTO-low", fg_n[["FTO_low"]], 25, 0)
med_of <- function(v, g) median(fto[[v]][fto$fto_group == g])
check("IFEX9 exon 1 median, high", med_of("cov_ifex9_exon1", "FTO_high"), 256, 0)
check("IFEX9 exon 1 median, low", med_of("cov_ifex9_exon1", "FTO_low"), 1.1, 1)
# Exons 1-8 and exon 9 medians aren't printed on the panel. (The table has 12.82
# and 87.6 for FTO-high, LEARNING.md Session 10 says 12.9 and 87; noted in Session 14.)
check("Exons 1-8 median, low", med_of("cov_exons_1_8", "FTO_low"), 18.2, 1)
check("IFEX9 junction min, high", min(fto$jx_ifex9_to_exon9[fto$fto_group == "FTO_high"]), 61, 0)
check("IFEX9 junction max, high", max(fto$jx_ifex9_to_exon9[fto$fto_group == "FTO_high"]), 434, 0)
fto_long <- fto %>%
  transmute(patient, fto_group,
            `FTO exons 1-8` = cov_exons_1_8, `FTO exon 9` = cov_exon_9,
            `IFEX9 exon 1` = cov_ifex9_exon1, `IFEX9 exon 1 to\nexon 9 junction reads` = jx_ifex9_to_exon9) %>%
  pivot_longer(-c(patient, fto_group)) %>%
  mutate(name = factor(name, levels = unique(name)),
         fto_group = factor(fto_group, levels = c("FTO_high", "FTO_low", "NPM1_wt")))
fto_fill <- c(FTO_high = col_high, FTO_low = col_low, NPM1_wt = col_other)
fto_edge <- c(FTO_high = col_high, FTO_low = edge_low, NPM1_wt = col_other)
fto_xlab <- c(FTO_high = sprintf("High\n(%d)", fg_n[["FTO_high"]]), FTO_low = sprintf("Low\n(%d)", fg_n[["FTO_low"]]),
              NPM1_wt = sprintf("NPM1 wt\n(%d)", fg_n[["NPM1_wt"]]))
p3d <- ggplot(fto_long, aes(fto_group, value)) +
  geom_jitter(aes(fill = fto_group, colour = fto_group), shape = 21, width = 0.22, height = 0,
              size = 0.9, stroke = 0.2) +
  geom_boxplot(outlier.shape = NA, fill = NA, width = 0.55, linewidth = 0.3) +
  facet_wrap(~ name, nrow = 1) +
  scale_fill_manual(values = fto_fill, guide = "none") +
  scale_colour_manual(values = fto_edge, guide = "none") +
  scale_x_discrete(labels = fto_xlab) +
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1, base = 10),
                     breaks = c(0, 1, 10, 100, 1000), labels = c("0", "1", "10", "100", "1,000")) +
  labs(x = NULL, y = "Depth-adjusted coverage / reads",
       title = "TCGA-LAML NPM1-mutated: the FTO rise is IFEX9, not FTO itself") +
  theme(axis.text.x = element_text(size = 6), strip.text = element_text(size = 6.5, lineheight = 0.9))

fig3 <- p3a / (p3b | p3c) / p3d + plot_layout(heights = c(1.15, 1.45, 0.85))
save_paper(fig3, "figure3_primitive_subgroup", h_mm = 225)

# ==== Figure 4. Drug response and survival (exploratory) =====================

# A. Drug volcano, NPM1-mutated vs wild-type (Session 8)
dn <- rd("drug_npm1_lm.csv")
check("Drugs tested", nrow(dn), 123, 0); check("Drugs at FDR < 0.05", sum(dn$fdr < 0.05), 0, 0)
dn$short <- sub(" \\(.*\\)$", "", dn$inhibitor)
dn$label <- ifelse(dn$short == "Venetoclax" | rank(dn$p) <= 6, dn$short, NA)
dn$ven <- dn$short == "Venetoclax"
p4a <- ggplot(dn, aes(auc_diff, -log10(fdr))) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey50", linewidth = 0.25) +
  geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.25) +
  geom_point(shape = 21, size = 1.3, fill = "grey75", colour = "grey30", stroke = 0.2) +
  geom_point(data = filter(dn, ven), shape = 21, size = 1.6, fill = "black", colour = "black") +
  geom_text_repel(aes(label = label), size = txt_small, family = font, min.segment.length = 0,
                  segment.size = 0.25, box.padding = 0.35, max.overlaps = Inf, na.rm = TRUE) +
  annotate("text", x = Inf, y = -log10(0.05), hjust = 1, vjust = -0.5, size = txt_small,
           family = font, colour = "grey40", label = "FDR 0.05") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.3, size = txt, family = font,
           label = sprintf("%d of %d drugs at FDR < 0.05", sum(dn$fdr < 0.05), nrow(dn))) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.1))) +
  labs(x = "AUC difference, NPM1-mutated - wild-type\n(adjusted for monocyte-like score)",
       y = "-log10 FDR",
       title = "Ex vivo drug response by NPM1 status",
       subtitle = "Lower AUC = more sensitive")

# B. Forest plot, high minus low within NPM1-mutated (Session 8)
dsu <- rd("drug_subgroup_npm1.csv")
check("Subgroup drugs at FDR < 0.05", sum(dsu$fdr < 0.05), 15, 0)
check("Venetoclax high - low", dsu$auc_diff[dsu$inhibitor == "Venetoclax"], -100, 0)
mer_drugs <- c("Sunitinib", "Sorafenib", "Quizartinib (AC220)")
stopifnot(all(mer_drugs %in% dsu$inhibitor))
check("Sunitinib FDR", dsu$fdr[dsu$inhibitor == "Sunitinib"], 0.32)
check("Sorafenib FDR", dsu$fdr[dsu$inhibitor == "Sorafenib"], 0.41)
check("Quizartinib FDR", dsu$fdr[dsu$inhibitor == "Quizartinib (AC220)"], 0.99)
das <- rd("drug_analysis_samples.csv")
n_drug_sub <- table(das$block_group[das$NPM1 == "mut"])
check("Drug set high", n_drug_sub[["high"]], 32, 0); check("Drug set low", n_drug_sub[["low"]], 48, 0)
fp <- dsu %>% filter(fdr < 0.05 | inhibitor %in% mer_drugs) %>%
  mutate(short = sub(" \\(.*\\)$", "", inhibitor),
         set = factor(ifelse(fdr < 0.05, "FDR < 0.05", "Mer et al. 2021 drugs"),
                      levels = c("FDR < 0.05", "Mer et al. 2021 drugs")),
         side = ifelse(auc_diff < 0, "high", "low"),
         fdr_lab = fmt_p(fdr)) %>%
  arrange(set, auc_diff)
fp$short <- factor(fp$short, levels = rev(unique(fp$short)))
fdr_x <- max(fp$ci_high) + 45
fdr_head <- data.frame(set = factor("FDR < 0.05", levels = levels(fp$set)), x = fdr_x,
                       short = levels(fp$short)[nlevels(fp$short)])
p4b <- ggplot(fp, aes(auc_diff, short)) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.25) +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), width = 0, orientation = "y", linewidth = 0.3) +
  geom_point(aes(fill = side, colour = side), shape = 21, size = 1.5, stroke = 0.3) +
  geom_text(aes(x = fdr_x, label = fdr_lab), size = txt_small, family = font, hjust = 1) +
  geom_text(data = fdr_head, aes(x = x, y = short), label = "FDR", size = txt_small, family = font,
            hjust = 1, vjust = -1.2, fontface = "bold") +
  facet_grid(rows = vars(set), scales = "free_y", space = "free_y") +
  scale_fill_manual(values = sub_fill, guide = "none") +
  scale_colour_manual(values = sub_edge, guide = "none") +
  coord_cartesian(xlim = c(min(fp$ci_low) - 5, fdr_x + 2), clip = "off") +
  labs(x = "AUC difference, high - low (95% CI)\n< 0: more sensitive in primitive-like", y = NULL,
       title = "Within NPM1-mutated: primitive-like vs committed-like",
       subtitle = sprintf("Analysis set %d high, %d low; n per drug varies", n_drug_sub[["high"]],
                          n_drug_sub[["low"]])) +
  theme(axis.text.y = element_text(size = 6), strip.text.y = element_text(size = 6.5, angle = 0, hjust = 0))

# Kaplan-Meier with numbers at risk underneath (same layout as 10_survival.R)
km_panel <- function(fit, levels, cols, ltys, title, ann, labels = levels, outline = NULL,
                     text_cols = cols, xmax = 60, by = 12) {
  if (is.null(names(labels))) names(labels) <- levels   # match labels by name, never by position
  strata <- rep(sub("^[^=]*=", "", names(fit$strata)), fit$strata)
  d <- data.frame(time = fit$time, surv = fit$surv, n.censor = fit$n.censor, group = strata)
  d <- rbind(data.frame(time = 0, surv = 1, n.censor = 0, group = unique(strata)), d)
  d$group <- factor(d$group, levels = levels)
  p <- ggplot(d, aes(time, surv, colour = group, linetype = group))
  # a thin darker line under a light curve, like the outline on light points
  for (g in names(outline)) {
    p <- p + geom_step(data = d[d$group == g, ], colour = outline[[g]], linewidth = 0.8,
                       linetype = "solid", show.legend = FALSE)
  }
  p <- p +
    geom_step(linewidth = 0.5) +
    geom_point(data = d[d$n.censor > 0, ], shape = 3, size = 0.7, show.legend = FALSE) +
    annotate("text", x = xmax, y = 1, hjust = 1, vjust = 1, size = txt_small, family = font,
             label = ann, lineheight = 0.95) +
    scale_colour_manual(values = cols, labels = labels, name = NULL) +
    scale_linetype_manual(values = ltys, labels = labels, name = NULL) +
    scale_x_continuous(breaks = seq(0, xmax, by)) +
    scale_y_continuous(labels = scales::percent, breaks = seq(0, 1, 0.25)) +
    coord_cartesian(xlim = c(0, xmax), ylim = c(0, 1)) +
    guides(colour = guide_legend(ncol = 2), linetype = guide_legend(ncol = 2)) +
    labs(x = "Months from diagnosis", y = "Overall survival", title = title) +
    theme(legend.position = "bottom", legend.key.width = unit(5, "mm"), legend.margin = margin(0, 0, 0, 0))
  r <- summary(fit, times = seq(0, xmax, by), extend = TRUE)
  risk <- data.frame(time = r$time, n.risk = r$n.risk, group = sub("^[^=]*=", "", as.character(r$strata)))
  risk$group <- factor(risk$group, levels = rev(levels))
  tab <- ggplot(risk, aes(time, group, label = n.risk, colour = group)) +
    geom_text(size = txt_small, family = font, show.legend = FALSE) +
    scale_colour_manual(values = text_cols) +
    scale_y_discrete(labels = labels) +
    scale_x_continuous(breaks = seq(0, xmax, by)) +
    coord_cartesian(xlim = c(0, xmax)) +
    labs(x = NULL, y = NULL, subtitle = "Number at risk") +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(),
          axis.text.y = element_text(size = 5.5), plot.subtitle = element_text(size = 6))
  wrap_elements(full = p / tab + plot_layout(heights = c(4, 1.3)))
}

sv <- rd("survival_beataml.csv")
cox <- rd("survival_cox.csv")
check("Survival set", nrow(sv), 275, 0); check("Deaths", sum(sv$death), 133, 0)
sv$age10 <- sv$ageAtDiagnosis / 10
# Refit only to check the saved HRs (survival_cox.csv) match; nothing new is reported
fav <- coxph(Surv(os_months, death) ~ I(grp4 == "NPM1-mut / ITD-neg") + age10, data = sv)
hr_fav <- exp(coef(fav))[1]
row_fav <- cox[cox$model == "NPM1-mut / ITD-neg vs everyone else + age" & cox$term == "npm1_favNPM1-mut / ITD-neg", ]
same("HR NPM1-mut / ITD-neg vs everyone else", hr_fav, row_fav$HR)
check("HR NPM1-mut / ITD-neg vs everyone else", row_fav$HR, 0.41)
check("HR CI low", row_fav$CI_low, 0.25); check("HR CI high", row_fav$CI_high, 0.67)
row_npm1 <- cox[cox$model == "NPM1 + age" & cox$term == "NPM1_statusmut", ]
check("HR NPM1 + age", row_npm1$HR, 0.70); check("p NPM1 + age", row_npm1$p, 0.06)

# C. Four NPM1 x FLT3-ITD groups: colour = NPM1, line type = FLT3-ITD
g4_lev <- c("NPM1-mut / ITD-neg", "NPM1-mut / ITD-pos", "NPM1-wt / ITD-neg", "NPM1-wt / ITD-pos")
g4_cols <- setNames(c(col_npm1, col_npm1, col_other, col_other), g4_lev)
g4_lty  <- setNames(c("solid", "22", "solid", "22"), g4_lev)
km_g4 <- survfit(Surv(os_months, death) ~ grp4, data = sv)
p4c <- km_panel(km_g4, g4_lev, g4_cols, g4_lty, "NPM1 x FLT3-ITD (n = 275)",
                sprintf("NPM1-mut / ITD-neg vs everyone else:\nHR %.2f (%.2f to %.2f), age-adjusted\nNPM1 alone: HR %.2f, p = %.2f",
                        row_fav$HR, row_fav$CI_low, row_fav$CI_high, row_npm1$HR, row_npm1$p))

# D. High vs low within NPM1-mutated
mut <- sv %>% filter(NPM1 == "mut") %>%
  inner_join(blk %>% select(dbgap_rnaseq_sample, block_group), by = "dbgap_rnaseq_sample") %>%
  mutate(block_group = factor(block_group, levels = c("low", "high")))
check("NPM1-mutated in survival", nrow(mut), 100, 0); check("NPM1-mutated deaths", sum(mut$death), 44, 0)
blk_fit <- coxph(Surv(os_months, death) ~ block_group + age10 + FLT3_ITD, data = mut)
row_blk <- cox[cox$model == "NPM1-mut: block high vs low + age + FLT3-ITD" & cox$term == "block_grouphigh", ]
same("HR high vs low", exp(coef(blk_fit))[["block_grouphigh"]], row_blk$HR)
check("HR high vs low", row_blk$HR, 1.19)
pw <- rd("survival_power.csv")
check("Power for HR 1.5", pw$power_HR_1.5[pw$cohort == "Beat AML"], 0.25)
km_blk <- survfit(Surv(os_months, death) ~ block_group, data = mut)
blk_lev <- c("high", "low")
p4d <- km_panel(km_blk, blk_lev, setNames(c(col_high, col_low), blk_lev), setNames(c("solid", "solid"), blk_lev),
                "Within NPM1-mutated (n = 100)",
                labels = c(high = sprintf("High, primitive-like (%d)", n_sub[["high"]]),
                           low = sprintf("Low, committed-like (%d)", n_sub[["low"]])),
                outline = c(low = edge_low), text_cols = setNames(c(col_high, edge_low), blk_lev),
                sprintf("High vs low: HR %.2f (%.2f to %.2f),\nadjusted for age and FLT3-ITD\n%d deaths: %.0f%% power for HR 1.5",
                        row_blk$HR, row_blk$CI_low, row_blk$CI_high, sum(mut$death),
                        100 * pw$power_HR_1.5[pw$cohort == "Beat AML"]))

fig4 <- (p4a | p4b) / (p4c | p4d) + plot_layout(heights = c(1, 1.2))
save_paper(fig4, "figure4_drugs_survival", h_mm = 190)

# ==== 5. Supplementary tables ================================================

sheets <- list(
  list("S1_cohort_beataml", "cohort_flow.csv", "Cohort flow, Beat AML 2.0 (steps 0 to 9)", "1, 3"),
  list("S2_cohort_tcga", "tcga_cohort_flow.csv", "Cohort flow, TCGA-LAML (steps 0 to 9)", "7"),
  list("S3_deseq2_npm1", "deseq2_npm1.csv", "Main DESeq2 model, NPM1-mutated vs wild-type, Beat AML (276 samples, 19,398 genes)", "3"),
  list("S4a_controls_main", "deseq2_positive_controls.csv", "Positive controls, main model (Beat AML)", "3"),
  list("S4b_controls_tcga", "tcga_positive_controls.csv", "Positive controls, TCGA-LAML model, with the Beat AML values", "7"),
  list("S4c_controls_celltype", "deseq2_celltype_controls.csv", "Positive controls in the cell-type sensitivity models (270 samples)", "6"),
  list("S4d_controls_subgroups", "deseq2_subgroups_controls.csv", "Positive controls in the subgroup model (high and low vs wild-type)", "5, 8"),
  list("S5a_gsea_verhaak", "gsea_verhaak_controls.csv", "GSEA, Verhaak et al. 2005 NPM1 gene sets (Beat AML)", "4"),
  list("S5b_gsea_hallmark", "gsea_hallmark.csv", "GSEA, MSigDB Hallmark (Beat AML)", "4"),
  list("S5c_gsea_gobp_collapsed", "gsea_gobp_collapsed.csv", "GSEA, GO Biological Process terms at FDR < 0.05, with collapsePathways result (Beat AML)", "4, 5"),
  list("S6a_sensitivity_flt3", "deseq2_sensitivity_flt3_summary.csv", "FLT3-ITD sensitivity model vs main model", "3"),
  list("S6b_sensitivity_celltype", "deseq2_celltype_summary.csv", "Cell-type sensitivity models vs main model (all_six is an upper bound only)", "6"),
  list("S7a_tcga_vs_beataml", "tcga_vs_beataml_summary.csv", "Beat AML vs TCGA-LAML agreement", "7"),
  list("S7b_tcga_three_group", "tcga_three_group_summary.csv", "HOXA9 and MEIS1 by group, TCGA-LAML", "7"),
  list("S7c_tcga_three_group_tests", "tcga_three_group_tests.csv", "HOXA9 and MEIS1 group comparisons, TCGA-LAML", "7"),
  list("S7d_tcga_deseq2_npm1", "tcga_deseq2_npm1.csv", "TCGA-LAML DESeq2 model, NPM1-mutated vs wild-type (111 patients, 26,069 genes)", "7"),
  list("S8a_subgroup_comutations", "subgroup_comutations.csv", "Co-mutations, NEK10/FTO/USP44-high vs -low (Beat AML NPM1-mutated)", "5, 8"),
  list("S8b_celltype_score_tests", "celltype_score_tests.csv", "Cell-type scores: NPM1-mutated vs wild-type, and high vs low", "6, 8"),
  list("S8c_subgroup_shared_counts", "deseq2_subgroups_shared_counts.csv", "Main-model genes shared by both subgroups or driven by one", "5, 8"),
  list("S9_fto_ifex9", "fto_key_feature_tests.csv", "FTO locus key features, TCGA-LAML (FTO-high 6, FTO-low 25, wild-type 80)", "10"),
  list("S10a_drugs_npm1", "drug_npm1_lm.csv", "Ex vivo drug response, NPM1-mutated vs wild-type, per drug", "8"),
  list("S10b_drugs_subgroup", "drug_subgroup_npm1.csv", "Ex vivo drug response, high vs low within NPM1-mutated, per drug (report only)", "8"),
  list("S11a_survival_cox", "survival_cox.csv", "Cox models, Beat AML and TCGA-LAML (exploratory)", "9"),
  list("S11b_survival_km", "survival_km_summary.csv", "Kaplan-Meier summaries (exploratory)", "9"),
  list("S11c_survival_power", "survival_power.csv", "Power of the high vs low survival comparison", "9"),
  list("S12a_classifier_model", "classifier_final_model.csv", "Final LASSO model, 33 genes (locked)", "13"),
  list("S12b_classifier_cv", "classifier_cv_summary.csv", "Nested cross-validation in Beat AML", "13"),
  list("S12c_classifier_cv_by_size", "classifier_cv_by_size.csv", "Nested CV with models capped at k genes", "13"),
  list("S12d_classifier_tcga", "classifier_tcga_performance.csv", "Locked models tested once in TCGA-LAML", "13"))

# Column descriptions. Shared names first; then "file:column" entries for
# names that mean something different in one table.
de_cols <- c(
  gene_id = "Ensembl gene ID without version (Beat AML: Ensembl 75, GRCh37; TCGA matched to it by ID)",
  symbol = "Gene symbol in the Beat AML annotation (Ensembl 75, 2014 names)",
  biotype = "Ensembl gene biotype",
  baseMean = "Mean DESeq2-normalized count across the samples in the model",
  log2FC = "log2 fold change, NPM1-mutated vs wild-type, apeglm-shrunken",
  lfcSE = "Standard error of the shrunken log2FC",
  log2FC_unshrunk = "log2 fold change before shrinkage (maximum likelihood)",
  stat = "Wald statistic (unshrunken log2FC / its SE); used to rank genes for GSEA",
  pvalue = "Wald test p-value",
  padj = "Benjamini-Hochberg adjusted p-value",
  rank = "Rank by absolute Wald statistic (1 = most significant)",
  sig = "TRUE if padj < 0.05 and abs(log2FC) > 1",
  direction = "up / down in NPM1-mutated if sig, otherwise ns",
  expected = "Expected direction in NPM1-mutated AML (up or down)",
  as_expected = "TRUE if padj < 0.05 and log2FC goes the expected way",
  pathway = "Gene set name (MSigDB 2026.1.Hs)",
  pval = "fgsea p-value",
  log2err = "fgsea's estimate of the error of log2(p-value)",
  ES = "Enrichment score (furthest point of the running sum from zero)",
  NES = "Normalized enrichment score; > 0 = towards higher in NPM1-mutated, < 0 = lower",
  size = "Number of genes of the set in the ranked list",
  leadingEdge = "Leading-edge genes (symbols, semicolon-separated)",
  inhibitor = "Drug name as in the Beat AML drug AUC file",
  drug_class = "Broad drug class, one per drug, by a fixed priority order ('Not in family table' = drug not in beataml_drug_families.xlsx)",
  ci_low = "Lower 95% confidence limit of auc_diff",
  ci_high = "Upper 95% confidence limit of auc_diff",
  p = "p-value",
  fdr = "Benjamini-Hochberg FDR across the drugs tested",
  cohort = "Beat AML or TCGA (TCGA-LAML)",
  model = "Model name",
  HR = "Hazard ratio",
  CI_low = "Lower 95% confidence limit of HR",
  CI_high = "Upper 95% confidence limit of HR",
  n = "Patients in the model",
  deaths = "Deaths in the model",
  comparison = "Which groups are compared",
  gene = "Gene symbol",
  group3 = "NPM1-mutated, KMT2A-rearranged or Other AML",
  wilcox_p = "Wilcoxon rank-sum p-value",
  median_diff = "Difference in medians, first group minus second group named in comparison",
  auc_mean = "Mean ROC AUC over 10 CV repeats (each from that repeat's out-of-fold predictions)",
  auc_min = "Lowest ROC AUC over the 10 repeats",
  auc_max = "Highest ROC AUC over the 10 repeats",
  genes_median = "Median number of genes per outer fit (50 fits)"
)
specific <- c(
  "cohort_flow.csv:step" = "Filter step (0 = all samples; 6 to 8 added after review; 9 = QC exclusion of BA2294R)",
  "cohort_flow.csv:description" = "What the step keeps",
  "cohort_flow.csv:samples" = "Samples left after the step",
  "cohort_flow.csv:patients" = "Patients left after the step",
  "cohort_flow.csv:npm1_positive" = "Samples NPM1-mutated (clinical consensus call)",
  "cohort_flow.csv:npm1_negative" = "Samples NPM1 wild-type",
  "cohort_flow.csv:npm1_mut_pct" = "% NPM1-mutated",
  "tcga_cohort_flow.csv:step" = "Filter step, in the same order as Beat AML",
  "tcga_cohort_flow.csv:filter" = "What the step keeps",
  "tcga_cohort_flow.csv:patients" = "Patients left after the step",
  "tcga_cohort_flow.csv:NPM1_mut" = "Patients NPM1-mutated (TCGA paper calls, see Session 7 rule)",
  "tcga_cohort_flow.csv:NPM1_wt" = "Patients NPM1 wild-type",
  "tcga_cohort_flow.csv:pct_mut" = "% NPM1-mutated",
  "tcga_cohort_flow.csv:note" = "Note on the step",
  "tcga_deseq2_npm1.csv:gene_id" = "Ensembl gene ID without version (GENCODE v36, GRCh38)",
  "tcga_deseq2_npm1.csv:symbol" = "Gene symbol in GENCODE v36 (may differ from the 2014 Beat AML symbol)",
  "tcga_deseq2_npm1.csv:biotype" = "GENCODE v36 gene type",
  "tcga_deseq2_npm1.csv:baseMean" = "Mean DESeq2-normalized count across the 111 TCGA main-cohort patients",
  "tcga_deseq2_npm1.csv:log2FC" = "log2 fold change, NPM1-mutated vs wild-type (~ sex + NPM1_status), apeglm-shrunken",
  "tcga_deseq2_npm1.csv:tested_in_beataml" = "TRUE if the same Ensembl ID was tested in the Beat AML model",
  "tcga_positive_controls.csv:log2FC_beataml" = "Shrunken log2FC in the Beat AML main model",
  "tcga_positive_controls.csv:padj_beataml" = "padj in the Beat AML main model",
  "tcga_positive_controls.csv:baseMean" = "Mean DESeq2-normalized count, TCGA",
  "tcga_positive_controls.csv:log2FC" = "Shrunken log2FC, TCGA",
  "tcga_positive_controls.csv:padj" = "padj, TCGA",
  "tcga_positive_controls.csv:rank" = "Rank by absolute Wald statistic, TCGA",
  "deseq2_celltype_controls.csv:model" = "main (276), main_270 (main design on the 270 with scores), two_score (+ HSC-like + Monocyte-like; the sensitivity model), all_six (upper bound only, not a result)",
  "deseq2_celltype_controls.csv:symbol" = "Gene symbol",
  "deseq2_celltype_controls.csv:pct_of_main_270_lfc" = "log2FC as a % of the main design's log2FC on the same 270 samples",
  "deseq2_celltype_controls.csv:survives" = "TRUE if padj < 0.05 in the expected direction",
  "deseq2_celltype_controls.csv:partly_circular" = "TRUE if the gene is one of the van Galen marker genes behind a score in this model (MEIS1 in HSC-like; TRH, POU4F1 in GMP-like)",
  "deseq2_celltype_controls.csv:core_control" = "TRUE for the five core positive controls (HOXA9, HOXA10, MEIS1, PBX3, CD34)",
  "deseq2_subgroups_controls.csv:symbol" = "Gene symbol",
  "deseq2_subgroups_controls.csv:log2FC_high_vs_wt" = "Shrunken log2FC, NEK10/FTO/USP44-high (primitive-like, n = 37) vs wild-type",
  "deseq2_subgroups_controls.csv:padj_high_vs_wt" = "padj, high vs wild-type",
  "deseq2_subgroups_controls.csv:ok_high" = "TRUE if padj < 0.05 in the expected direction, high vs wild-type",
  "deseq2_subgroups_controls.csv:log2FC_low_vs_wt" = "Shrunken log2FC, NEK10/FTO/USP44-low (committed-like, n = 63) vs wild-type",
  "deseq2_subgroups_controls.csv:padj_low_vs_wt" = "padj, low vs wild-type",
  "deseq2_subgroups_controls.csv:ok_low" = "TRUE if padj < 0.05 in the expected direction, low vs wild-type",
  "deseq2_subgroups_controls.csv:log2FC_main" = "Shrunken log2FC in the main model",
  "gsea_verhaak_controls.csv:padj" = "BH-adjusted p-value within the Verhaak collection",
  "gsea_hallmark.csv:padj" = "BH-adjusted p-value within Hallmark",
  "gsea_gobp_collapsed.csv:padj" = "BH-adjusted p-value within GO BP (all 3,540 sets tested)",
  "gsea_gobp_collapsed.csv:main_pathway" = "TRUE if kept as a main term by fgsea::collapsePathways",
  "gsea_gobp_collapsed.csv:parent" = "For a folded term, the main term it was folded into",
  "gsea_gobp_collapsed.csv:n_children" = "For a main term, how many terms were folded into it",
  "deseq2_sensitivity_flt3_summary.csv:top50_kept_in_top50" = "Main-model top 50 genes still in the top 50 with FLT3-ITD in the model",
  "deseq2_sensitivity_flt3_summary.csv:top50_still_sig" = "Main-model top 50 genes still significant (padj < 0.05, abs(log2FC) > 1)",
  "deseq2_sensitivity_flt3_summary.csv:sig_main" = "Significant genes, main model",
  "deseq2_sensitivity_flt3_summary.csv:sig_flt3" = "Significant genes, model with FLT3-ITD",
  "deseq2_sensitivity_flt3_summary.csv:sig_shared" = "Significant in both",
  "deseq2_sensitivity_flt3_summary.csv:lfc_correlation" = "Pearson correlation of shrunken log2FC, all genes",
  "deseq2_celltype_summary.csv:model" = "main, main_270, two_score (sensitivity model) or all_six (upper bound only)",
  "deseq2_celltype_summary.csv:role" = "What the model is for",
  "deseq2_celltype_summary.csv:top50_in_top50" = "Main-model top 50 genes still in this model's top 50",
  "deseq2_celltype_summary.csv:top50_still_sig" = "Main-model top 50 genes still significant (padj < 0.05, abs(log2FC) > 1)",
  "deseq2_celltype_summary.csv:top50_same_direction" = "Main-model top 50 genes with the same sign",
  "deseq2_celltype_summary.csv:n_sig" = "Significant genes in this model",
  "deseq2_celltype_summary.csv:n_up" = "Of which up",
  "deseq2_celltype_summary.csv:n_down" = "Of which down",
  "deseq2_celltype_summary.csv:sig_shared_with_main" = "Significant here and in the main model",
  "deseq2_celltype_summary.csv:lfc_pearson_vs_main" = "Pearson r of shrunken log2FC vs the main model (276)",
  "deseq2_celltype_summary.csv:lfc_spearman_vs_main" = "Spearman rho of shrunken log2FC vs the main model (276)",
  "deseq2_celltype_summary.csv:lfc_pearson_vs_main_270" = "Pearson r of shrunken log2FC vs the main design on the same 270",
  "deseq2_celltype_summary.csv:lfc_slope_vs_main_270" = "Slope through zero of log2FC vs main_270, across the 1,019 main-model genes",
  "deseq2_celltype_summary.csv:top50_dropped" = "Main-model top 50 genes no longer significant (semicolon-separated)",
  "tcga_vs_beataml_summary.csv:metric" = "What is measured (genes tested in both cohorts unless stated)",
  "tcga_vs_beataml_summary.csv:value" = "Value",
  "tcga_three_group_summary.csv:n" = "Patients in the group",
  "tcga_three_group_summary.csv:median_vst" = "Median VST expression (blind VST, TCGA main cohort + KMT2A-r)",
  "tcga_three_group_summary.csv:q25" = "25th percentile of VST",
  "tcga_three_group_summary.csv:q75" = "75th percentile of VST",
  "tcga_three_group_tests.csv:median_diff" = "Difference in median VST, first group minus second",
  "subgroup_comutations.csv:feature" = "Mutation (exome calls; FLT3_ITD_pos from the clinical assay; cohesin = SMC1A, SMC3, RAD21 or STAG2)",
  "subgroup_comutations.csv:n_tested" = "Samples with data (97 with exome data, 100 for FLT3-ITD)",
  "subgroup_comutations.csv:low_n" = "Mutated samples in the low (committed-like) group",
  "subgroup_comutations.csv:low_of" = "Low-group samples with data",
  "subgroup_comutations.csv:high_n" = "Mutated samples in the high (primitive-like) group",
  "subgroup_comutations.csv:high_of" = "High-group samples with data",
  "subgroup_comutations.csv:low_pct" = "% mutated, low group",
  "subgroup_comutations.csv:high_pct" = "% mutated, high group",
  "subgroup_comutations.csv:odds_ratio" = "Fisher's exact test odds ratio, high vs low",
  "subgroup_comutations.csv:fisher_p" = "Fisher's exact test p-value",
  "subgroup_comutations.csv:fisher_padj" = "BH-adjusted across the 10 features",
  "celltype_score_tests.csv:comparison" = "NPM1-mut vs wt (reference = wild-type) or NPM1-high vs NPM1-low (reference = low)",
  "celltype_score_tests.csv:cell_type" = "Bottomly et al. 2022 malignant cell-type score",
  "celltype_score_tests.csv:n_ref" = "Samples with scores, reference group",
  "celltype_score_tests.csv:n_alt" = "Samples with scores, other group",
  "celltype_score_tests.csv:median_ref" = "Median score, reference group",
  "celltype_score_tests.csv:median_alt" = "Median score, other group",
  "celltype_score_tests.csv:median_diff" = "Median other minus median reference",
  "celltype_score_tests.csv:glass_delta" = "Glass's delta: mean difference in reference-group SDs",
  "celltype_score_tests.csv:padj" = "BH-adjusted across the six scores, within the comparison",
  "deseq2_subgroups_shared_counts.csv:view" = "fold change (2x rule) or significance (padj < 0.05 in the main-model direction)",
  "deseq2_subgroups_shared_counts.csv:class" = "shared, mostly high, mostly low (fold change view); both, high only, low only, neither (significance view)",
  "deseq2_subgroups_shared_counts.csv:down" = "Genes down in the main model",
  "deseq2_subgroups_shared_counts.csv:up" = "Genes up in the main model",
  "deseq2_subgroups_shared_counts.csv:total" = "All genes in the class (of 1,019)",
  "fto_key_feature_tests.csv:measure" = "cov_* = mean depth-adjusted coverage of the feature; log2_* = log2 ratio to FTO exons 1-8; jx_ifex9_to_exon9 = IFEX9 exon 1 to FTO exon 9 junction reads",
  "fto_key_feature_tests.csv:med_high" = "Median in FTO-high (n = 6), on log2(value + 0.05) for cov_* and jx_*, as is for log2_*",
  "fto_key_feature_tests.csv:med_low" = "Median in FTO-low (n = 25), same scale",
  "fto_key_feature_tests.csv:med_wt" = "Median in NPM1 wild-type (n = 80), same scale",
  "fto_key_feature_tests.csv:p_high_vs_low" = "Wilcoxon p-value, high vs low (smallest possible with 6 vs 25 is 2.7e-6)",
  "fto_key_feature_tests.csv:p_low_vs_wt" = "Wilcoxon p-value, low vs wild-type",
  "drug_npm1_lm.csv:n_wt" = "Wild-type samples tested with the drug",
  "drug_npm1_lm.csv:n_mut" = "NPM1-mutated samples tested with the drug",
  "drug_npm1_lm.csv:median_auc_wt" = "Median AUC, wild-type (lower = more sensitive)",
  "drug_npm1_lm.csv:median_auc_mut" = "Median AUC, NPM1-mutated",
  "drug_npm1_lm.csv:auc_diff" = "NPM1-mutated minus wild-type mean AUC from auc ~ NPM1_status + Monocyte_like",
  "drug_npm1_lm.csv:se" = "Standard error of auc_diff",
  "drug_npm1_lm.csv:p" = "p-value of auc_diff",
  "drug_npm1_lm.csv:auc_diff_unadjusted" = "Same difference from auc ~ NPM1_status",
  "drug_npm1_lm.csv:p_unadjusted" = "p-value of the unadjusted difference",
  "drug_npm1_lm.csv:monocyte_coef" = "AUC change per unit of monocyte-like score (same model)",
  "drug_npm1_lm.csv:monocyte_p" = "p-value of monocyte_coef",
  "drug_npm1_lm.csv:fdr" = "BH FDR across the 123 drugs tested",
  "drug_npm1_lm.csv:families" = "Target families from beataml_drug_families.xlsx (semicolon-separated; empty if the drug isn't in it)",
  "drug_subgroup_npm1.csv:n_low" = "Low (committed-like) samples tested with the drug",
  "drug_subgroup_npm1.csv:n_high" = "High (primitive-like) samples tested with the drug",
  "drug_subgroup_npm1.csv:median_auc_low" = "Median AUC, low group",
  "drug_subgroup_npm1.csv:median_auc_high" = "Median AUC, high group",
  "drug_subgroup_npm1.csv:auc_diff" = "High minus low mean AUC from auc ~ block_group (the primary model)",
  "drug_subgroup_npm1.csv:p" = "p-value of auc_diff",
  "drug_subgroup_npm1.csv:auc_diff_adj_flt3" = "Same difference adjusted for FLT3-ITD",
  "drug_subgroup_npm1.csv:p_adj_flt3" = "p-value, FLT3-ITD-adjusted",
  "drug_subgroup_npm1.csv:auc_diff_adj_monocyte" = "Same difference adjusted for the monocyte-like score",
  "drug_subgroup_npm1.csv:p_adj_monocyte" = "p-value, monocyte-adjusted",
  "drug_subgroup_npm1.csv:fdr" = "BH FDR across the 113 drugs tested, primary model only",
  "survival_cox.csv:model" = "Cox model: groups compared and covariates (age per 10 years, FLT3-ITD)",
  "survival_cox.csv:term" = "Model term (e.g. NPM1_statusmut = NPM1-mutated vs wild-type; age10 = per 10 years of age)",
  "survival_km_summary.csv:group" = "Group within the comparison",
  "survival_km_summary.csv:n" = "Patients in the group",
  "survival_km_summary.csv:deaths" = "Deaths in the group",
  "survival_km_summary.csv:median_OS_months" = "Kaplan-Meier median overall survival, months (NA = not reached)",
  "survival_km_summary.csv:median_CI_low" = "Lower 95% confidence limit of the median",
  "survival_km_summary.csv:median_CI_high" = "Upper 95% confidence limit of the median (NA = not reached)",
  "survival_km_summary.csv:OS_1yr" = "Kaplan-Meier overall survival at 1 year (proportion)",
  "survival_km_summary.csv:OS_2yr" = "Kaplan-Meier overall survival at 2 years (proportion)",
  "survival_km_summary.csv:OS_3yr" = "Kaplan-Meier overall survival at 3 years (proportion)",
  "survival_power.csv:deaths" = "Deaths in the comparison",
  "survival_power.csv:fraction_in_group" = "Share of patients in the high group",
  "survival_power.csv:power_HR_1.5" = "Power to detect HR 1.5 at alpha 0.05 (Schoenfeld 1983)",
  "survival_power.csv:power_HR_2" = "Power to detect HR 2",
  "survival_power.csv:power_HR_3" = "Power to detect HR 3",
  "survival_power.csv:min_HR_80pct_power" = "Smallest HR detectable with 80% power",
  "survival_power.csv:deaths_needed_HR_1.5" = "Deaths needed for 80% power at HR 1.5",
  "survival_power.csv:deaths_needed_HR_2" = "Deaths needed for 80% power at HR 2",
  "classifier_final_model.csv:term" = "(Intercept) or the gene's Ensembl ID",
  "classifier_final_model.csv:symbol" = "Gene symbol (Beat AML annotation)",
  "classifier_final_model.csv:coef" = "LASSO coefficient, log-odds of NPM1-mutated per VST unit",
  "classifier_final_model.csv:beat_mean" = "Mean VST in Beat AML (276); TCGA was rescaled to this mean",
  "classifier_final_model.csv:beat_sd" = "SD of VST in Beat AML; TCGA was rescaled to this SD",
  "classifier_final_model.csv:per_sd" = "coef x beat_sd: log-odds per Beat AML SD, comparable across genes",
  "classifier_final_model.csv:direction" = "higher or lower in NPM1-mutated",
  "classifier_final_model.csv:deseq2_rank" = "Rank in the main DESeq2 model",
  "classifier_final_model.csv:cv_pct_of_fits" = "% of the 50 outer CV fits that picked the gene",
  "classifier_final_model.csv:entry_order" = "Order in which the gene enters the LASSO path as the penalty is relaxed (1 = first)",
  "classifier_cv_summary.csv:model" = "LASSO with lambda.1se (main), LASSO with lambda.min, or HOXA9 + CD34",
  "classifier_cv_summary.csv:genes_min" = "Fewest genes in an outer fit",
  "classifier_cv_summary.csv:genes_max" = "Most genes in an outer fit",
  "classifier_cv_summary.csv:sens" = "Sensitivity at probability 0.5, averaged over repeats",
  "classifier_cv_summary.csv:spec" = "Specificity at probability 0.5, averaged over repeats",
  "classifier_cv_summary.csv:acc" = "Accuracy at probability 0.5, averaged over repeats",
  "classifier_cv_by_size.csv:max_genes" = "Largest number of genes allowed (LASSO path capped at k genes)",
  "classifier_tcga_performance.csv:model" = "lasso (33 genes) or hoxa9_cd34",
  "classifier_tcga_performance.csv:scaling" = "rescaled (main: per-gene mean and SD matched to Beat AML) or VST only",
  "classifier_tcga_performance.csv:auc" = "ROC AUC in TCGA (111 patients)",
  "classifier_tcga_performance.csv:auc_lo" = "Lower DeLong 95% confidence limit",
  "classifier_tcga_performance.csv:auc_hi" = "Upper DeLong 95% confidence limit",
  "classifier_tcga_performance.csv:sensitivity" = "Sensitivity at probability 0.5",
  "classifier_tcga_performance.csv:specificity" = "Specificity at probability 0.5",
  "classifier_tcga_performance.csv:accuracy" = "Accuracy at probability 0.5",
  "classifier_tcga_performance.csv:mean_predicted" = "Mean predicted probability of NPM1-mutated",
  "classifier_tcga_performance.csv:observed" = "Observed share NPM1-mutated",
  "classifier_tcga_performance.csv:wt_called_mut" = "Wild-type patients called NPM1-mutated",
  "classifier_tcga_performance.csv:mut_called_wt" = "NPM1-mutated patients called wild-type"
)
describe <- function(file, col) {
  key <- paste0(file, ":", col)
  if (key %in% names(specific)) return(specific[[key]])
  if (col %in% names(de_cols)) return(de_cols[[col]])
  NA_character_
}

wb <- createWorkbook()
hdr <- createStyle(textDecoration = "bold", border = "bottom", fgFill = "#F2F2F2")
fmt_dec <- createStyle(numFmt = "0.000")
fmt_sci <- createStyle(numFmt = "0.00E+00")
addWorksheet(wb, "README")
tables_index <- list(); columns_index <- list()
for (s in sheets) {
  name <- s[[1]]; file <- s[[2]]
  stopifnot(nchar(name) <= 31)
  d <- rd(file)
  desc <- vapply(names(d), function(cn) describe(file, cn), "")
  if (anyNA(desc)) stop("No description for ", file, ": ", paste(names(d)[is.na(desc)], collapse = ", "))
  addWorksheet(wb, name)
  writeData(wb, name, d, headerStyle = hdr)
  freezePane(wb, name, firstRow = TRUE)
  setColWidths(wb, name, cols = seq_along(d),
               widths = pmin(pmax(nchar(names(d)) + 2, 10), ifelse(names(d) %in% c("leadingEdge", "families"), 40, 45)))
  # Display formats only: the stored values are the CSV values, unrounded
  for (j in seq_along(d)) {
    x <- d[[j]]
    if (!is.numeric(x) || all(is.na(x)) || all(x == round(x), na.rm = TRUE)) next
    st <- if (grepl("^(p|pval|pvalue|padj|fdr|fisher_p|fisher_padj|wilcox_p|monocyte_p|spearman_p)$|^p_", names(d)[j]))
      fmt_sci else fmt_dec
    addStyle(wb, name, st, rows = seq_len(nrow(d)) + 1, cols = j, gridExpand = TRUE, stack = TRUE)
  }
  tables_index[[name]] <- data.frame(sheet = name, contents = s[[3]], source = paste0("results/", file),
                                     session = s[[4]], rows = nrow(d), columns = ncol(d))
  columns_index[[name]] <- data.frame(sheet = name, column = names(d), meaning = unname(desc))
}
ti <- bind_rows(tables_index); ci <- bind_rows(columns_index)
intro <- c("Supplementary tables: NPM1-mutated vs wild-type AML (Beat AML 2.0, replication in TCGA-LAML).",
           "Every sheet is a copy of a table in results/ of the project repository, values unchanged (cells are only displayed rounded).",
           "Lower part: what each column means. AUC in the drug sheets is ex vivo drug AUC (lower = more sensitive); AUC in the classifier sheets is ROC AUC.",
           "Data: Beat AML 2.0 (Bottomly et al., Cancer Cell 2022; Tyner et al., Nature 2018; CC-BY-4.0) and TCGA-LAML (NEJM 2013).")
writeData(wb, "README", data.frame(x = intro), colNames = FALSE, startRow = 1)
writeData(wb, "README", "Tables", startRow = length(intro) + 2)
addStyle(wb, "README", createStyle(textDecoration = "bold", fontSize = 12), rows = length(intro) + 2, cols = 1)
writeData(wb, "README", ti, startRow = length(intro) + 3, headerStyle = hdr)
col_start <- length(intro) + 3 + nrow(ti) + 2
writeData(wb, "README", "Columns", startRow = col_start)
addStyle(wb, "README", createStyle(textDecoration = "bold", fontSize = 12), rows = col_start, cols = 1)
writeData(wb, "README", ci, startRow = col_start + 1, headerStyle = hdr)
setColWidths(wb, "README", cols = 1:6, widths = c(28, 70, 45, 10, 8, 8))
saveWorkbook(wb, "results/paper/supplementary_tables.xlsx", overwrite = TRUE)
cat("Saved results/paper/supplementary_tables.xlsx:", length(sheets), "tables +", "README,",
    nrow(ci), "column descriptions\n")

writeLines(capture.output(sessionInfo()), "results/sessionInfo_14_paper_figures.txt")
cat("Done.\n")
