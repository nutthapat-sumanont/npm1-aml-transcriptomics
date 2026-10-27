# 07_cell_composition.R
# How much of the NPM1 signal is just cell maturity?
#
# Uses the six malignant cell-type scores from Bottomly et al., Cancer Cell
# 2022 (file beataml2_manuscript_vg_cts.xlsx). Each score summarises the top 30
# marker genes of one AML cell type from van Galen et al., Cell 2019
# (single-cell RNA-seq): HSC-like, Progenitor-like, GMP-like, Promono-like,
# Monocyte-like, cDC-like. A score is the first principal component of those
# 30 genes (like a WGCNA eigengene), computed separately for bone marrow and
# blood samples, so it says "how much does this sample look like cell type X,
# compared with other samples of the same specimen type".
#
#  1. Join the scores to the main cohort and report how many samples have one.
#  2. Get the van Galen top 30 genes per cell type and check which positive
#     controls and top 50 NPM1 genes sit inside a signature.
#  3. Compare each score between NPM1-mutated and wild-type (main cohort).
#  4. Compare each score between the NEK10/FTO/USP44-high and -low
#     NPM1-mutated samples (Session 4 split).
#  5. Figure: figures/07_celltype_scores.png.
#  6. How the scores relate to each other and to the positive-control genes.
#
# Run from the repo root:  Rscript scripts/07_cell_composition.R
# Needs data/ (01_download.R), results/cohort.csv (02), results/deseq2_npm1.csv
# (04) and results/block_nek10_fto_scores.csv (05).

suppressPackageStartupMessages({
  library(DESeq2)
  library(data.table)
  library(dplyr)
  library(readxl)
  library(ggplot2)
  library(patchwork)
})

set.seed(2026)
dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

counts_file <- "data/beataml_waves1to4_counts_dbgap.txt"
scores_file <- "data/beataml2_manuscript_vg_cts.xlsx"
vg_file     <- "data/van_galen_2019_table_s3.xlsx"
vg_url      <- "https://ars.els-cdn.com/content/image/1-s2.0-S0092867419300947-mmc3.xlsx"
padj_cut    <- 0.05
controls    <- c("HOXA9", "HOXA10", "MEIS1", "PBX3", "CD34")
# Cell types from most primitive to most mature
cell_types  <- c("HSC-like", "Progenitor-like", "GMP-like", "Promono-like",
                 "Monocyte-like", "cDC-like")
npm1_cols   <- c(wt = "grey55", mut = "#D55E00")
block_cols  <- c(low = "#009E73", high = "#D55E00")

save_fig <- function(p, file, w = 7, h = 5) {
  ggsave(file.path("figures", file), p, width = w, height = h, dpi = 300)
  cat("Saved figures/", file, "\n", sep = "")
}

# ---- 1. Scores and coverage --------------------------------------------------
# Bottomly et al. scored the 560 RNA-seq samples they used in the manuscript
# (clinical column manuscript_rnaseq == "yes"): one sample per patient, no
# remission samples. They scored Waves 1+2 and then projected Waves 3+4 onto
# the same PCs. No wave "Both" sample was scored.

scores <- as.data.frame(read_excel(scores_file))
stopifnot(identical(names(scores), c("dbgap_rnaseq_sample", "GMP-like", "HSC-like",
                                     "Monocyte-like", "Progenitor-like",
                                     "Promono-like", "cDC-like")))
cat("Cell-type score file:", nrow(scores), "samples x", length(cell_types), "scores\n")

block <- read.csv("results/block_nek10_fto_scores.csv")
cohort <- read.csv("results/cohort.csv") %>%
  filter(main_cohort) %>%
  left_join(block[, c("dbgap_rnaseq_sample", "block_group")], by = "dbgap_rnaseq_sample") %>%
  left_join(scores, by = "dbgap_rnaseq_sample") %>%
  mutate(NPM1_status = factor(NPM1, levels = c("wt", "mut")),
         block_group = factor(block_group, levels = c("low", "high")),
         has_scores  = !is.na(`HSC-like`))
stopifnot(nrow(cohort) == 276)

coverage <- cohort %>%
  group_by(NPM1_status) %>%
  summarise(n = n(), with_scores = sum(has_scores), .groups = "drop") %>%
  bind_rows(cohort %>% filter(NPM1 == "mut") %>% group_by(block_group) %>%
              summarise(n = n(), with_scores = sum(has_scores), .groups = "drop") %>%
              mutate(NPM1_status = paste0("mut_", block_group)) %>%
              select(-block_group)) %>%
  bind_rows(data.frame(NPM1_status = "main cohort", n = nrow(cohort),
                       with_scores = sum(cohort$has_scores)))
cat("\nMain-cohort samples with cell-type scores:\n"); print(as.data.frame(coverage))
cat("\nMissing scores (all have manuscript_rnaseq = 'no' in the clinical file):\n")
print(cohort %>% filter(!has_scores) %>%
        select(dbgap_subject_id, dbgap_rnaseq_sample, NPM1, cohort, specimenType),
      row.names = FALSE)
write.csv(coverage, "results/celltype_score_coverage.csv", row.names = FALSE)

sc <- cohort %>% filter(has_scores)

# ---- 2. Which genes are in each signature? -----------------------------------
# van Galen et al. 2019 Table S3, "Tumor-derived, per cell-type" columns. For
# GMP-like that column says "see combined signature", so (like Bottomly et
# al.) we take the tumor-derived combined GMP-like column. Downloaded once to
# data/; the 180 genes are saved to results/ so later sessions don't need it.

if (!file.exists(vg_file)) {
  cat("\nDownloading van Galen et al. 2019 Table S3 ...\n")
  download.file(vg_url, vg_file, mode = "wb", quiet = TRUE)
}
vg_raw <- as.data.frame(read_excel(vg_file, col_names = FALSE, .name_repair = "minimal"))
vg_cols <- c("HSC-like" = 8, "Progenitor-like" = 9, "GMP-like" = 6,
             "Promono-like" = 11, "Monocyte-like" = 12, "cDC-like" = 13)
stopifnot(vg_raw[2, 8] == "HSC-like", vg_raw[2, 6] == "GMP-like", vg_raw[2, 13] == "cDC-like")
vg_genes <- bind_rows(lapply(names(vg_cols), function(ct)
  data.frame(cell_type = ct, rank = 1:30, symbol = vg_raw[3:32, vg_cols[[ct]]])))
write.csv(vg_genes, "results/van_galen_top30_genes.csv", row.names = FALSE)

main <- read.csv("results/deseq2_npm1.csv")
top50 <- head(main, 50)
in_sig <- function(g) {
  hit <- vg_genes$cell_type[vg_genes$symbol == g]
  if (length(hit)) paste(hit, collapse = ";") else NA_character_
}
overlap <- data.frame(symbol = unique(c(controls, top50$symbol))) %>%
  mutate(positive_control = symbol %in% controls,
         main_rank = main$rank[match(symbol, main$symbol)],
         signature = vapply(symbol, in_sig, ""))
cat("\nPositive controls inside a van Galen signature:\n")
print(overlap %>% filter(positive_control) %>% select(symbol, main_rank, signature),
      row.names = FALSE)
cat("Top 50 NPM1 genes inside a signature:",
    sum(!is.na(overlap$signature) & overlap$main_rank <= 50), "of 50\n")
print(overlap %>% filter(!is.na(signature), main_rank <= 50) %>%
        select(symbol, main_rank, signature), row.names = FALSE)
write.csv(overlap, "results/celltype_signature_overlap.csv", row.names = FALSE)

# ---- 3 & 4. Score differences ------------------------------------------------
# Wilcoxon rank-sum test per score, BH-corrected across the six scores within
# each comparison. Effect sizes: difference in medians (in score units) and
# Glass's delta = (mean of group 2 - mean of group 1) / SD of group 1, the
# effect size Bottomly et al. report (group 1 = wild-type, or NPM1-low).

compare_scores <- function(d, grp, ref, alt, label) {
  bind_rows(lapply(cell_types, function(ct) {
    a <- d[[ct]][d[[grp]] == ref]
    b <- d[[ct]][d[[grp]] == alt]
    data.frame(comparison = label, cell_type = ct,
               n_ref = length(a), n_alt = length(b),
               median_ref = median(a), median_alt = median(b),
               median_diff = median(b) - median(a),
               glass_delta = (mean(b) - mean(a)) / sd(a),
               wilcox_p = wilcox.test(b, a)$p.value)
  })) %>% mutate(padj = p.adjust(wilcox_p, "BH"))
}

tests <- bind_rows(
  compare_scores(sc, "NPM1_status", "wt", "mut", "NPM1-mut vs wt"),
  compare_scores(filter(sc, NPM1 == "mut"), "block_group", "low", "high",
                 "NPM1-high vs NPM1-low"))
cat("\nCell-type scores, NPM1-mutated vs wild-type, and high vs low within NPM1-mutated:\n")
print(tests %>% mutate(across(c(median_ref, median_alt, median_diff, glass_delta), ~ round(.x, 2)),
                       across(c(wilcox_p, padj), ~ signif(.x, 2))), row.names = FALSE)
write.csv(tests, "results/celltype_score_tests.csv", row.names = FALSE)

# ---- 5. Figure ---------------------------------------------------------------

long <- sc %>%
  select(dbgap_rnaseq_sample, NPM1_status, block_group, all_of(cell_types)) %>%
  tidyr::pivot_longer(all_of(cell_types), names_to = "cell_type", values_to = "score") %>%
  mutate(cell_type = factor(cell_type, levels = cell_types))

p_lab <- function(cmp) {
  t <- tests %>% filter(comparison == cmp)
  setNames(sprintf("%s\nGlass d %+.2f, padj %s", t$cell_type, t$glass_delta,
                   ifelse(t$padj < 1e-3, formatC(t$padj, format = "e", digits = 0),
                          formatC(t$padj, format = "f", digits = 3))),
           t$cell_type)
}
box_panel <- function(d, x, cols, labs_x, cmp, title) {
  ggplot(d, aes(.data[[x]], score, fill = .data[[x]])) +
    geom_hline(yintercept = 0, colour = "grey80") +
    geom_boxplot(outlier.shape = NA, alpha = 0.5, width = 0.6) +
    geom_jitter(aes(colour = .data[[x]]), width = 0.18, size = 0.7, alpha = 0.7) +
    facet_wrap(~ cell_type, nrow = 1, labeller = as_labeller(p_lab(cmp))) +
    scale_fill_manual(values = cols, guide = "none") +
    scale_colour_manual(values = cols, guide = "none") +
    scale_x_discrete(labels = labs_x) +
    labs(x = NULL, y = "Cell-type score (PC1)", title = title) +
    theme_bw(base_size = 9) +
    theme(strip.text = element_text(size = 7.5))
}
n_g <- table(sc$NPM1_status)
n_b <- table(sc$block_group[sc$NPM1 == "mut"])
pA <- box_panel(long, "NPM1_status", npm1_cols,
                c(wt = sprintf("wt\n(%d)", n_g[["wt"]]), mut = sprintf("NPM1-mut\n(%d)", n_g[["mut"]])),
                "NPM1-mut vs wt", "A. NPM1-mutated vs wild-type, main cohort")
pB <- box_panel(filter(long, NPM1_status == "mut"), "block_group", block_cols,
                c(low = sprintf("low\n(%d)", n_b[["low"]]), high = sprintf("high\n(%d)", n_b[["high"]])),
                "NPM1-high vs NPM1-low",
                "B. Within NPM1-mutated: NEK10/FTO/USP44-high (primitive-like) vs -low (committed-like)")
p <- pA / pB + plot_annotation(
  title = "Malignant cell-type scores (Bottomly et al. 2022, van Galen et al. 2019 signatures)",
  subtitle = paste0("Cell types ordered from most primitive (left) to most mature (right). ",
                    "Wilcoxon test, BH across the six scores in each panel.\n",
                    "Glass d = mean difference in SDs of the reference group (wt, or NPM1-low)."),
  theme = theme(plot.title = element_text(size = 11), plot.subtitle = element_text(size = 8.5)))
save_fig(p, "07_celltype_scores.png", w = 11, h = 7.5)

# ---- 6. Scores vs each other and vs the positive-control genes ---------------
# The six scores aren't independent: they mostly line up along one
# primitive-to-mature axis. Spearman correlations between scores, and between
# each score and the positive-control genes (blind VST, same low-count rule as
# before), in all samples and within each NPM1 group.

cor_scores <- cor(sc[, cell_types], method = "spearman")
cat("\nSpearman correlation between scores (main cohort with scores):\n")
print(round(cor_scores, 2))
write.csv(round(cor_scores, 3), "results/celltype_score_correlations.csv")
pca <- prcomp(sc[, cell_types], scale. = TRUE)
cat("Share of variance in the six scores explained by their PC1:",
    round(summary(pca)$importance[2, 1], 2), "\n")

# How well do the scores predict NPM1 status on their own? (logistic
# regression, in-sample AUC: 0.5 = no better than a coin, 1 = perfect). If a
# set of scores nearly identifies NPM1 status, putting that set in the DESeq2
# design would mostly regress out NPM1 itself (over-adjustment).
auc <- function(p, y) {
  r <- rank(p)
  (sum(r[y]) - sum(y) * (sum(y) + 1) / 2) / (sum(y) * sum(!y))
}
score_sets <- c(as.list(setNames(cell_types, cell_types)), list(
  "HSC-like + Monocyte-like" = c("HSC-like", "Monocyte-like"),
  "HSC-like + Progenitor-like + Monocyte-like" = c("HSC-like", "Progenitor-like", "Monocyte-like"),
  "GMP-like + Promono-like" = c("GMP-like", "Promono-like"),
  "all six" = cell_types))
pred_npm1 <- bind_rows(lapply(names(score_sets), function(n) {
  f <- glm(reformulate(paste0("`", score_sets[[n]], "`"), "NPM1 == 'mut'"),
           data = sc, family = binomial)
  data.frame(scores = n, auc = auc(fitted(f), sc$NPM1 == "mut"),
             mcfadden_r2 = 1 - f$deviance / f$null.deviance)
}))
cat("\nHow well each set of scores predicts NPM1 status (logistic regression, in-sample):\n")
print(pred_npm1 %>% mutate(across(-scores, ~ round(.x, 3))), row.names = FALSE)
write.csv(pred_npm1, "results/celltype_score_predict_npm1.csv", row.names = FALSE)

raw <- fread(counts_file, check.names = FALSE)
counts <- as.matrix(raw[, sc$dbgap_rnaseq_sample, with = FALSE])
rownames(counts) <- raw$stable_id
storage.mode(counts) <- "integer"
counts <- counts[rowSums(counts >= 10) >= 100, ]
vsd <- vst(counts, blind = TRUE)
sym <- setNames(raw$display_label, raw$stable_id)
genes <- c(controls, "NEK10", "FTO", "USP44")
gene_ids <- names(sym)[match(genes, sym)]
stopifnot(!anyNA(gene_ids), all(gene_ids %in% rownames(vsd)))

gene_cor <- bind_rows(lapply(c("all", "wt", "mut"), function(g) {
  keep <- if (g == "all") rep(TRUE, nrow(sc)) else sc$NPM1 == g
  bind_rows(lapply(seq_along(genes), function(i)
    data.frame(samples = g, gene = genes[i],
               t(sapply(cell_types, function(ct)
                 cor(vsd[gene_ids[i], keep], sc[[ct]][keep], method = "spearman"))),
               check.names = FALSE)))
}))
cat("\nSpearman correlation, gene (VST) vs cell-type score:\n")
print(gene_cor %>% mutate(across(all_of(cell_types), ~ round(.x, 2))), row.names = FALSE)
write.csv(gene_cor, "results/celltype_score_gene_correlations.csv", row.names = FALSE)

# How much of the NPM1 difference in each control gene goes with the scores?
# R squared of NPM1 status alone, the six scores alone, and both together
# (linear model on VST). If adding NPM1 on top of the scores still explains a
# lot, the gene's NPM1 difference isn't just cell maturity.
r2 <- bind_rows(lapply(seq_along(genes), function(i) {
  y <- vsd[gene_ids[i], ]
  f_s <- as.formula(paste("y ~", paste0("`", cell_types, "`", collapse = " + ")))
  data.frame(gene = genes[i],
             r2_npm1   = summary(lm(y ~ NPM1_status, data = sc))$r.squared,
             r2_scores = summary(lm(f_s, data = sc))$r.squared,
             r2_both   = summary(lm(update(f_s, . ~ . + NPM1_status), data = sc))$r.squared)
})) %>% mutate(npm1_added_r2 = r2_both - r2_scores)
cat("\nVariance in each gene explained by NPM1 status, the six scores, or both (VST, linear model):\n")
print(r2 %>% mutate(across(-gene, ~ round(.x, 2))), row.names = FALSE)
write.csv(r2, "results/celltype_score_gene_r2.csv", row.names = FALSE)

# ---- 7. Sensitivity DESeq2 with cell-type scores -----------------------------
# Approved in Session 6. Extra models only; the main model (04_deseq2.R) is
# unchanged. All three fits use the same 270 samples (those with scores):
#
#   main_270:  ~ wave + specimen_type + sex + NPM1_status
#              (the main design refit on the 270, so we can tell what changes
#              because 6 samples are gone from what changes because of the
#              adjustment)
#   two_score: ~ wave + specimen_type + sex + HSC_like + Monocyte_like + NPM1_status
#              THE sensitivity model. HSC-like and Monocyte-like are the two
#              ends of the hierarchy; together they predict NPM1 with AUC 0.75.
#   all_six:   ~ wave + specimen_type + sex + all six scores + NPM1_status
#              UPPER BOUND ONLY, not a result. The six scores predict NPM1
#              with AUC 0.93, so this model removes NPM1 biology along with
#              cell maturity. It shows the most "cell type" could explain.
#
# Wave has only two levels here: none of the wave "Both" samples has scores.
# Scores are scaled to mean 0, SD 1 before fitting (a numerical convenience;
# it doesn't change the NPM1 coefficient). Same 19,398 genes as the main
# model, so every gene can be compared one to one.

suppressPackageStartupMessages(library(BiocParallel))
register(MulticoreParam(max(1, min(4, parallel::detectCores()))))
lfc_cut   <- 1
coef_name <- "NPM1_status_mut_vs_wt"

cd <- sc %>%
  mutate(specimen_type = factor(ifelse(specimenType == "Bone Marrow Aspirate",
                                       "BM", "Blood"), levels = c("BM", "Blood")),
         wave = droplevels(factor(cohort, levels = c("Waves1+2", "Waves3+4", "Both"))),
         sex  = factor(consensus_sex, levels = c("Female", "Male")))
score_vars <- gsub("-", "_", cell_types)          # HSC-like -> HSC_like
for (i in seq_along(cell_types)) cd[[score_vars[i]]] <- as.numeric(scale(cd[[cell_types[i]]]))
rownames(cd) <- cd$dbgap_rnaseq_sample
stopifnot(nrow(cd) == 270, nlevels(cd$wave) == 2)
cat("\nDESeq2 sensitivity: ", nrow(cd), " samples (", sum(cd$NPM1 == "mut"),
    " NPM1-mut, ", sum(cd$NPM1 == "wt"), " wt); wave levels: ",
    paste(levels(cd$wave), collapse = ", "), "\n", sep = "")

counts_dds <- as.matrix(raw[, cd$dbgap_rnaseq_sample, with = FALSE])
rownames(counts_dds) <- raw$stable_id
storage.mode(counts_dds) <- "integer"
counts_dds <- counts_dds[main$gene_id, ]          # same genes as the main model
stopifnot(identical(colnames(counts_dds), rownames(cd)))

designs <- list(
  main_270  = ~ wave + specimen_type + sex + NPM1_status,
  two_score = ~ wave + specimen_type + sex + HSC_like + Monocyte_like + NPM1_status,
  all_six   = as.formula(paste("~ wave + specimen_type + sex +",
                               paste(score_vars, collapse = " + "), "+ NPM1_status")))
model_signatures <- list(main_270 = character(0),
                         two_score = c("HSC-like", "Monocyte-like"),
                         all_six = cell_types)

fit_one <- function(design) {
  cat("Fitting", deparse(design, width.cutoff = 500), "...\n")
  dds <- DESeqDataSetFromMatrix(counts_dds, colData = cd, design = design)
  dds <- DESeq(dds, parallel = TRUE, quiet = TRUE)
  res <- results(dds, name = coef_name, alpha = padj_cut)
  shr <- lfcShrink(dds, coef = coef_name, type = "apeglm", res = res,
                   parallel = TRUE, quiet = TRUE)
  data.frame(gene_id = rownames(res), log2FC = shr$log2FoldChange,
             stat = res$stat, padj = res$padj) %>%
    arrange(desc(abs(stat))) %>%
    mutate(rank = row_number(),
           sig = !is.na(padj) & padj < padj_cut & abs(log2FC) > lfc_cut)
}
fits <- lapply(designs, fit_one)

# One wide table: main model (276, from 04_deseq2.R) plus the three new fits
wide <- main %>% select(gene_id, symbol, log2FC_main = log2FC, padj_main = padj,
                        rank_main = rank, sig_main = sig)
for (m in names(fits)) {
  f <- fits[[m]]
  names(f)[-1] <- paste0(names(f)[-1], "_", m)
  wide <- left_join(wide, f, by = "gene_id")
}
wide$in_signature <- vapply(wide$symbol, in_sig, "")
stopifnot(nrow(wide) == nrow(main), !anyNA(wide$log2FC_all_six))
write.csv(wide %>% arrange(rank_main), "results/deseq2_celltype_adjusted.csv", row.names = FALSE)
cat("Saved results/deseq2_celltype_adjusted.csv (", nrow(wide), " genes)\n", sep = "")

# ---- 8. What survives? -------------------------------------------------------
# "Survives" = padj < 0.05 in the expected direction (same rule as Session 3).
# A positive control is flagged "partly circular" in a model if it is one of
# the 30 genes behind a score in that model's design: the model then partly
# adjusts the gene for itself. Among all positive controls (HOXA9, HOXA10,
# HOXB cluster, MEIS1, PBX3, CD34) only MEIS1 is in a signature (HSC-like).

hoxb <- sort(grep("^HOXB[0-9]+$", main$symbol, value = TRUE))
pc_all <- c("HOXA9", "HOXA10", hoxb, "MEIS1", "PBX3", "CD34")
in_any <- intersect(pc_all, vg_genes$symbol)
cat("\nPositive controls (incl. HOXB cluster) inside a van Galen signature:",
    paste(in_any, collapse = ", "), "\n")

pc_long <- bind_rows(lapply(c("main", names(fits)), function(m) {
  w <- wide %>% filter(symbol %in% pc_all)
  lfc  <- w[[paste0("log2FC_", m)]]
  padj <- w[[paste0("padj_", m)]]
  expected <- ifelse(w$symbol == "CD34", "down", "up")
  circular <- if (m == "main") rep(FALSE, nrow(w)) else
    vapply(w$symbol, function(g) any(vg_genes$symbol[vg_genes$cell_type %in%
                                                       model_signatures[[m]]] == g), TRUE)
  data.frame(model = m, symbol = w$symbol, expected = expected,
             log2FC = lfc, padj = padj,
             pct_of_main_270_lfc = if (m == "main") NA else
               100 * lfc / w$log2FC_main_270,
             survives = !is.na(padj) & padj < padj_cut &
               sign(lfc) == ifelse(expected == "up", 1, -1),
             partly_circular = circular,
             core_control = w$symbol %in% controls)
})) %>% mutate(model = factor(model, levels = c("main", names(fits))),
               symbol = factor(symbol, levels = pc_all)) %>%
  arrange(symbol, model)
write.csv(pc_long, "results/deseq2_celltype_controls.csv", row.names = FALSE)
cat("\nCore positive controls by model (log2FC, padj; * = partly circular):\n")
print(pc_long %>% filter(core_control) %>%
        mutate(log2FC = round(log2FC, 2), padj = signif(padj, 2),
               pct_of_main_270_lfc = round(pct_of_main_270_lfc),
               flag = ifelse(partly_circular, "*", "")) %>%
        select(symbol, model, log2FC, padj, pct_of_main_270_lfc, survives, flag),
      row.names = FALSE)
cat("HOXB genes surviving: ", paste(sapply(names(fits), function(m)
  sprintf("%s %d/%d", m, sum(pc_long$survives[pc_long$model == m & grepl("^HOXB", pc_long$symbol)]),
          length(hoxb))), collapse = "; "), "\n", sep = "")

# Top 50 genes and fold-change agreement, each model vs the main model (276)
# and vs the main design on the same 270 samples.
top50_ids <- head(main$gene_id, 50)
summ7 <- bind_rows(lapply(names(fits), function(m) {
  f <- fits[[m]]
  lf <- wide[[paste0("log2FC_", m)]]
  data.frame(model = m,
             role = c(main_270 = "same samples, main design",
                      two_score = "sensitivity model",
                      all_six = "UPPER BOUND, not a result")[[m]],
             top50_in_top50 = sum(top50_ids %in% head(f$gene_id, 50)),
             top50_still_sig = sum(top50_ids %in% f$gene_id[f$sig]),
             top50_same_direction = sum(sign(wide$log2FC_main[1:50]) ==
                                        sign(lf[match(top50_ids, wide$gene_id)])),
             n_sig = sum(f$sig), n_up = sum(f$sig & f$log2FC > 0),
             n_down = sum(f$sig & f$log2FC < 0),
             sig_shared_with_main = length(intersect(f$gene_id[f$sig], main$gene_id[main$sig])),
             lfc_pearson_vs_main = cor(wide$log2FC_main, lf),
             lfc_spearman_vs_main = cor(wide$log2FC_main, lf, method = "spearman"),
             lfc_pearson_vs_main_270 = cor(wide$log2FC_main_270, lf),
             # slope through 0 over the 1,019 main-model genes, adjusted model
             # vs main design on the same 270: 1 = fold changes unchanged,
             # 0.6 = on average 40% smaller
             lfc_slope_vs_main_270 = unname(coef(lm(lf[wide$sig_main] ~
                                                      0 + wide$log2FC_main_270[wide$sig_main]))[1]),
             top50_dropped = paste(main$symbol[1:50][!top50_ids %in% f$gene_id[f$sig]],
                                   collapse = ";"))
}))
summ7 <- bind_rows(data.frame(model = "main", role = "main model, 276 samples (04_deseq2.R)",
                              top50_in_top50 = 50, top50_still_sig = sum(main$sig[1:50]),
                              top50_same_direction = 50,
                              n_sig = sum(main$sig), n_up = sum(main$direction == "up"),
                              n_down = sum(main$direction == "down"),
                              sig_shared_with_main = sum(main$sig),
                              lfc_pearson_vs_main = 1, lfc_spearman_vs_main = 1,
                              lfc_pearson_vs_main_270 = cor(wide$log2FC_main, wide$log2FC_main_270),
                              lfc_slope_vs_main_270 = NA, top50_dropped = ""),
                   summ7)
cat("\nSummary, each model vs the main model (top 50 = main-model top 50 by |Wald stat|):\n")
print(summ7 %>% mutate(across(where(is.double), ~ round(.x, 3))) %>% select(-top50_dropped),
      row.names = FALSE)
for (m in names(fits)) cat("Top 50 genes no longer significant in", m, ":",
                           summ7$top50_dropped[summ7$model == m], "\n")
write.csv(summ7, "results/deseq2_celltype_summary.csv", row.names = FALSE)

# Figure: shrunken log2FC, main model on 270 (x) vs each adjusted model (y).
# Colours = main-model up / down / ns; labels = positive controls (* = partly
# circular in that model) and top 10 genes.
lab_genes <- unique(c(controls, head(main$symbol, 10)))
plot_df <- bind_rows(lapply(c("two_score", "all_six"), function(m) {
  circ <- if (m == "two_score") intersect(lab_genes, vg_genes$symbol[vg_genes$cell_type %in% model_signatures$two_score])
          else intersect(lab_genes, vg_genes$symbol)
  data.frame(model = m, x = wide$log2FC_main_270, y = wide[[paste0("log2FC_", m)]],
             symbol = wide$symbol,
             direction = factor(ifelse(!wide$sig_main, "ns",
                                       ifelse(wide$log2FC_main > 0, "up", "down")),
                                levels = c("up", "down", "ns"))) %>%
    mutate(label = ifelse(symbol %in% lab_genes,
                          paste0(symbol, ifelse(symbol %in% circ, "*", "")), NA))
}))
s_two <- summ7[summ7$model == "two_score", ]
s_six <- summ7[summ7$model == "all_six", ]
plot_df$model <- factor(plot_df$model, levels = c("two_score", "all_six"), labels = c(
  sprintf("Sensitivity model: + HSC-like + Monocyte-like\nr = %.2f, slope %.2f; top 50 still significant: %d/50",
          cor(wide$log2FC_main_270, wide$log2FC_two_score), s_two$lfc_slope_vs_main_270,
          s_two$top50_still_sig),
  sprintf("UPPER BOUND (not a result): + all six scores\nr = %.2f, slope %.2f; top 50 still significant: %d/50",
          cor(wide$log2FC_main_270, wide$log2FC_all_six), s_six$lfc_slope_vs_main_270,
          s_six$top50_still_sig)))
dir_cols <- c(up = "#D55E00", down = "#0072B2", ns = "grey80")
lim <- range(c(plot_df$x, plot_df$y))
p <- ggplot(plot_df, aes(x, y)) +
  geom_hline(yintercept = 0, colour = "grey70") +
  geom_vline(xintercept = 0, colour = "grey70") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  geom_point(data = ~ filter(.x, direction == "ns"), colour = "grey80", size = 0.4) +
  geom_point(data = ~ filter(.x, direction != "ns"), aes(colour = direction),
             size = 0.6, alpha = 0.7) +
  ggrepel::geom_label_repel(aes(label = label), size = 2.5, max.overlaps = Inf,
                            min.segment.length = 0, fill = alpha("white", 0.8),
                            na.rm = TRUE) +
  facet_wrap(~ model) +
  scale_colour_manual(values = dir_cols, name = "Main model (276)",
                      breaks = c("up", "down")) +
  coord_equal(xlim = lim, ylim = lim) +
  labs(x = "log2FC NPM1-mut vs wt, main design on the same 270 samples (apeglm)",
       y = "log2FC NPM1-mut vs wt, adjusted for cell-type scores (apeglm)",
       title = "How much of the NPM1 signal survives adjustment for cell maturity?",
       subtitle = paste0("Dashed line = no change. Slope = average fold-change ratio across the 1,019 main-model genes (1 = unchanged).\n",
                         "* = partly circular: the gene is one of the 30 genes behind a score in that model.")) +
  theme_bw(base_size = 9) + theme(legend.position = "bottom")
save_fig(p, "07_celltype_adjusted_lfc.png", w = 11, h = 6.8)

writeLines(c(paste("Designs:", vapply(designs, function(d) deparse(d, width.cutoff = 500), "")),
             capture.output(sessionInfo())),
           "results/sessionInfo_07_cell_composition.txt")
cat("\nDone.\n")
