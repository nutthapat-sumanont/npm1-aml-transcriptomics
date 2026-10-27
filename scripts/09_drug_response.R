# 09_drug_response.R
# Session 8: the pharma angle. Beat AML tested each patient's leukaemia cells
# ex vivo against a panel of drugs and summarised each dose-response curve as
# an area under the curve (AUC). Lower AUC = fewer cells survive = more
# sensitive.
#
#  1. Match drug results to the specimen used for RNA-seq, and count matches.
#  2. Every drug: AUC in NPM1-mutated vs wild-type, linear model
#     auc ~ NPM1_status + Monocyte_like (Session 6 cell-type score). BH across
#     drugs. Drugs with fewer than 20 tested samples per group are skipped.
#  3. "Drug volcano": AUC difference vs -log10 FDR, coloured by drug family,
#     venetoclax labelled.
#  4. Venetoclax in detail: AUC by NPM1 status and against the monocyte-like
#     score, one figure.
#  5. Is any menin inhibitor in the panel?
#  6. Within NPM1-mutated: NEK10/FTO/USP44-high vs -low (three-gene split),
#     every drug, report only, with FDR. Close look at sorafenib, sunitinib and
#     quizartinib (Mer et al., Nat Commun 2021).
#  7. HOX program score (mean VST of HOXA9, HOXA10, MEIS1, PBX3): does it track
#     response to any drug beyond NPM1 status?
#
# Run from the repo root:  Rscript scripts/09_drug_response.R   (about 2 minutes)
# Needs data/ (01_download.R, incl. beataml_drug_families.xlsx),
# results/cohort.csv (02) and results/block_nek10_fto_scores.csv (05).

suppressPackageStartupMessages({
  library(DESeq2)
  library(data.table)
  library(dplyr)
  library(readxl)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

set.seed(2026)
dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

auc_file    <- "data/beataml_probit_curve_fits_v4_dbgap.txt"
family_file <- "data/beataml_drug_families.xlsx"
scores_file <- "data/beataml2_manuscript_vg_cts.xlsx"
counts_file <- "data/beataml_waves1to4_counts_dbgap.txt"
min_per_group <- 20     # skip a drug with fewer tested samples in either group
fdr_cut       <- 0.05
hox_genes     <- c("HOXA9", "HOXA10", "MEIS1", "PBX3")
npm1_cols     <- c(wt = "grey55", mut = "#D55E00")
block_cols    <- c(low = "#009E73", high = "#D55E00")

save_fig <- function(p, file, w = 7, h = 5) {
  ggsave(file.path("figures", file), p, width = w, height = h, dpi = 300)
  cat("Saved figures/", file, "\n", sep = "")
}

# ---- 1. Match drug results to the RNA-seq specimen ---------------------------
# Each row of the AUC file is one drug on one specimen. The file carries the
# RNA-seq sample ID of that same specimen (blank if the specimen had no
# RNA-seq), so joining on dbgap_rnaseq_sample means the drug and RNA results
# come from the same tube of cells. We use only those exact matches. Patients
# whose drug data are on a different specimen are counted but not used.

auc <- fread(auc_file)
stopifnot(all(auc$type == "single-agent"))
cohort <- read.csv("results/cohort.csv") %>% filter(main_cohort)
stopifnot(nrow(cohort) == 276)

matched <- cohort$dbgap_rnaseq_sample %in% auc$dbgap_rnaseq_sample
other_spec <- unique(auc[dbgap_subject_id %in% cohort$dbgap_subject_id[!matched],
                         .(dbgap_subject_id, dbgap_rnaseq_sample, dbgap_dnaseq_sample)])
match_tab <- data.frame(
  step = c("Main cohort", "Drug data on the same specimen as RNA-seq",
           "No drug data on that specimen",
           "  ...but drug data on another specimen of the same patient (not used)",
           "Same specimen and a monocyte-like score (analysis set)"),
  n = NA, NPM1_mut = NA, NPM1_wt = NA)
scores <- as.data.frame(read_excel(scores_file))
cohort$has_drug  <- matched
cohort$has_score <- cohort$dbgap_rnaseq_sample %in% scores$dbgap_rnaseq_sample
sets <- list(rep(TRUE, nrow(cohort)), matched, !matched,
             !matched & cohort$dbgap_subject_id %in% other_spec$dbgap_subject_id,
             matched & cohort$has_score)
for (i in seq_along(sets)) {
  s <- sets[[i]]
  match_tab[i, 2:4] <- c(sum(s), sum(s & cohort$NPM1 == "mut"), sum(s & cohort$NPM1 == "wt"))
}
write.csv(match_tab, "results/drug_specimen_matching.csv", row.names = FALSE)
cat("Specimen matching:\n"); print(match_tab, row.names = FALSE)
cat("\nPatients with drug data only on another specimen:\n"); print(other_spec)

# Analysis set: same-specimen drug data AND a monocyte-like score, so every
# model below (with or without the score) uses the same samples.
block <- read.csv("results/block_nek10_fto_scores.csv")
an <- cohort %>%
  filter(has_drug, has_score) %>%
  left_join(scores[, c("dbgap_rnaseq_sample", "Monocyte-like", "HSC-like")],
            by = "dbgap_rnaseq_sample") %>%
  rename(Monocyte_like = `Monocyte-like`, HSC_like = `HSC-like`) %>%
  left_join(block[, c("dbgap_rnaseq_sample", "block_group")], by = "dbgap_rnaseq_sample") %>%
  mutate(NPM1_status = factor(NPM1, levels = c("wt", "mut")),
         block_group = factor(block_group, levels = c("low", "high")),
         FLT3_ITD = factor(FLT3_ITD, levels = c("negative", "positive")))

d <- auc[dbgap_rnaseq_sample %in% an$dbgap_rnaseq_sample]
stopifnot(!anyDuplicated(d[, .(dbgap_rnaseq_sample, inhibitor)]))
cat("\nAUC rows in the analysis set:", nrow(d), "| non-converged fits dropped:",
    sum(!d$converged), "\n")
d <- d[converged == TRUE]
d <- merge(d[, .(dbgap_rnaseq_sample, inhibitor, auc)], as.data.table(an),
           by = "dbgap_rnaseq_sample")

# ---- Drug families -----------------------------------------------------------
# beataml_drug_families.xlsx, sheet drug_family: one row per drug-family pair
# (many drugs hit several families; 49 drugs aren't listed at all). For the
# colour in the volcano each drug gets ONE broad class, taken from its families
# in this priority order (most specific first). The full family list per drug
# is kept in the table.
fam <- as.data.table(read_excel(family_file, sheet = "drug_family"))
class_levels <- c("BCL-2 family", "Type III RTK (FLT3, KIT, PDGFR, CSF1R)",
                  "Other receptor tyrosine kinases", "Non-receptor tyrosine kinases",
                  "PI3K / AKT / mTOR", "Other kinases", "Non-kinase targets",
                  "Not in family table")
family_class <- function(f) {
  case_when(
    grepl("Bcl-2", f) ~ class_levels[1],
    grepl("^Type III RTKs", f) ~ class_levels[2],
    grepl("^Type [IVX]+ RTKs", f) ~ class_levels[3],
    grepl("^(Abl|Ack|Csk|Fer|Src|Syk|Tec) family|Janus kinase", f) ~ class_levels[4],
    grepl("Akt|PIKK|Phosphatidyl|Glycerophospholipid", f) ~ class_levels[5],
    grepl(paste0("HDAC|Histone|Proteasome|Heat shock|Chromatin|G protein|GPCR|",
                 "Carbonic|Lanosterol|Peptidase|PSMD14|ADP-ribose|Peptidyl-prolyl|",
                 "channel"), f) ~ class_levels[7],
    TRUE ~ class_levels[6])
}
fam[, class := factor(family_class(family), levels = class_levels)]
drug_class <- fam[, .(families = paste(sort(unique(family)), collapse = "; "),
                      drug_class = levels(class)[min(as.integer(class))]),
                  by = inhibitor]
drugs <- data.table(inhibitor = sort(unique(auc$inhibitor)))
drugs <- merge(drugs, drug_class, by = "inhibitor", all.x = TRUE)
drugs[is.na(drug_class), drug_class := "Not in family table"]
drugs[, drug_class := factor(drug_class, levels = class_levels)]
write.csv(drugs, "results/drug_family_classes.csv", row.names = FALSE)
cat("\nDrug classes (all", nrow(drugs), "drugs in the AUC file):\n"); print(table(drugs$drug_class))

# ---- 2. NPM1-mutated vs wild-type, every drug -------------------------------
# For each drug: auc ~ NPM1_status + Monocyte_like. The NPM1 coefficient is the
# difference in mean AUC (mutated minus wild-type) at the same monocyte-like
# score. Negative = NPM1-mutated cells are more sensitive. The monocyte-like
# score is in because maturity affects drug response (venetoclax most of all)
# and the two NPM1 subgroups differ in it. An unadjusted difference is kept
# next to it for comparison.

fit_drug <- function(dd, formula, term) {
  m <- lm(formula, data = dd)
  co <- summary(m)$coefficients
  ci <- confint(m)[term, ]
  data.frame(estimate = co[term, "Estimate"], se = co[term, "Std. Error"],
             ci_low = ci[1], ci_high = ci[2], p = co[term, "Pr(>|t|)"])
}

counts_by <- function(dd, var) {
  tab <- dcast(dd[, .N, by = c("inhibitor", var)], as.formula(paste("inhibitor ~", var)),
               value.var = "N", fill = 0)
  tab
}

n_npm1 <- counts_by(d, "NPM1_status")
setnames(n_npm1, c("wt", "mut"), c("n_wt", "n_mut"))
keep2 <- n_npm1[n_wt >= min_per_group & n_mut >= min_per_group, inhibitor]
cat("\nStep 2: drugs tested in the analysis set:", nrow(n_npm1),
    "| kept (>=", min_per_group, "per group):", length(keep2),
    "| skipped:", nrow(n_npm1) - length(keep2), "\n")

res2 <- rbindlist(lapply(keep2, function(dr) {
  dd <- d[inhibitor == dr]
  adj <- fit_drug(dd, auc ~ NPM1_status + Monocyte_like, "NPM1_statusmut")
  raw <- fit_drug(dd, auc ~ NPM1_status, "NPM1_statusmut")
  mono <- fit_drug(dd, auc ~ NPM1_status + Monocyte_like, "Monocyte_like")
  data.table(inhibitor = dr,
             median_auc_wt = median(dd$auc[dd$NPM1_status == "wt"]),
             median_auc_mut = median(dd$auc[dd$NPM1_status == "mut"]),
             auc_diff = adj$estimate, se = adj$se, ci_low = adj$ci_low,
             ci_high = adj$ci_high, p = adj$p,
             auc_diff_unadjusted = raw$estimate, p_unadjusted = raw$p,
             monocyte_coef = mono$estimate, monocyte_p = mono$p)
}))
res2[, fdr := p.adjust(p, method = "BH")]
res2 <- merge(n_npm1, res2, by = "inhibitor")
res2 <- merge(res2, drugs, by = "inhibitor")
setorder(res2, p)
fwrite(res2, "results/drug_npm1_lm.csv")
cat("Drugs at FDR <", fdr_cut, ":", sum(res2$fdr < fdr_cut),
    "(lower AUC in NPM1-mutated:", sum(res2$fdr < fdr_cut & res2$auc_diff < 0),
    "| higher:", sum(res2$fdr < fdr_cut & res2$auc_diff > 0), ")\n")
print(head(res2[, .(inhibitor, n_wt, n_mut, auc_diff, ci_low, ci_high, p, fdr,
                    auc_diff_unadjusted, drug_class)], 20), digits = 3)

# ---- 3. Drug volcano ---------------------------------------------------------
class_cols <- setNames(c("#CC79A7", "#D55E00", "#E69F00", "#56B4E9", "#009E73",
                         "#0072B2", "#7F3C8D", "grey70"), class_levels)
vol <- copy(res2)
vol[, label := ifelse(inhibitor == "Venetoclax" | (fdr < fdr_cut & rank(p) <= 10),
                      sub(" \\(.*\\)$", "", inhibitor), NA)]
vol[, is_ven := inhibitor == "Venetoclax"]
p_vol <- ggplot(vol, aes(auc_diff, -log10(fdr))) +
  geom_hline(yintercept = -log10(fdr_cut), linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = 0, colour = "grey80") +
  geom_point(aes(fill = drug_class), shape = 21, size = 2.4, colour = "grey20",
             stroke = 0.25, alpha = 0.9) +
  geom_point(data = vol[is_ven == TRUE], shape = 21, size = 4, colour = "black",
             stroke = 1, fill = NA) +
  geom_text_repel(aes(label = label, fontface = ifelse(is_ven, "bold", "plain")),
                  size = 3, max.overlaps = Inf, min.segment.length = 0,
                  segment.colour = "grey40", box.padding = 0.4, na.rm = TRUE) +
  scale_fill_manual(values = class_cols, drop = FALSE, name = "Drug class") +
  labs(title = "Ex vivo drug response, NPM1-mutated vs wild-type AML",
       subtitle = sprintf(paste0("%d drugs (>= %d tested samples per group); AUC ~ NPM1 status + ",
                                 "monocyte-like score\nDashed line: FDR = 0.05. Left = more ",
                                 "sensitive in NPM1-mutated"), nrow(vol), min_per_group),
       x = "AUC difference, NPM1-mutated minus wild-type (adjusted)",
       y = "-log10 FDR (BH across drugs)") +
  theme_bw(base_size = 10) +
  theme(legend.position = "right", legend.text = element_text(size = 8))
save_fig(p_vol, "09_drug_volcano.png", w = 9.5, h = 6)

# ---- 4. Venetoclax -----------------------------------------------------------
# Venetoclax blocks BCL-2. Immature (stem/progenitor-like) AML cells depend on
# BCL-2 to stay alive; monocytic cells lean on MCL-1 instead and lose BCL-2, so
# they are resistant (Pei et al. Cancer Discov 2020; Kuusanmaki et al.
# Haematologica 2020; Zhang et al. Nat Cancer 2020, in Beat AML). So the
# monocyte-like score should track venetoclax AUC, and it could hide or
# exaggerate an NPM1 effect. Plotted on the analysis set.
ven <- d[inhibitor == "Venetoclax"]
ven_fit <- lm(auc ~ NPM1_status + Monocyte_like, data = ven)
ven_int <- lm(auc ~ NPM1_status * Monocyte_like, data = ven)
ven_rho <- ven[, .(n = .N, spearman_rho = cor(auc, Monocyte_like, method = "spearman"),
                   spearman_p = cor.test(auc, Monocyte_like, method = "spearman",
                                         exact = FALSE)$p.value,
                   median_auc = median(auc), median_monocyte = median(Monocyte_like)),
               by = NPM1_status]
ven_sum <- rbind(
  data.frame(quantity = c("NPM1 mut - wt, adjusted for monocyte-like",
                          "NPM1 mut - wt, unadjusted",
                          "Monocyte-like (per score unit), adjusted for NPM1",
                          "Interaction NPM1 x monocyte-like"),
             estimate = c(coef(ven_fit)["NPM1_statusmut"], res2[inhibitor == "Venetoclax", auc_diff_unadjusted],
                          coef(ven_fit)["Monocyte_like"], coef(ven_int)["NPM1_statusmut:Monocyte_like"]),
             p = c(summary(ven_fit)$coefficients[c("NPM1_statusmut", "Monocyte_like"), 4][1],
                   res2[inhibitor == "Venetoclax", p_unadjusted],
                   summary(ven_fit)$coefficients["Monocyte_like", 4],
                   summary(ven_int)$coefficients["NPM1_statusmut:Monocyte_like", 4])))
ven_sum$fdr_step2 <- c(res2[inhibitor == "Venetoclax", fdr], NA, NA, NA)
write.csv(ven_sum, "results/drug_venetoclax_summary.csv", row.names = FALSE)
write.csv(ven_rho, "results/drug_venetoclax_by_npm1.csv", row.names = FALSE)
cat("\nVenetoclax:\n"); print(ven_sum, digits = 3); print(ven_rho, digits = 3)

ven[, NPM1_lab := factor(ifelse(NPM1_status == "mut", "NPM1-mutated", "NPM1 wild-type"),
                         levels = c("NPM1 wild-type", "NPM1-mutated"))]
lab_cols <- setNames(npm1_cols, c("NPM1 wild-type", "NPM1-mutated"))
pa <- ggplot(ven, aes(NPM1_lab, auc, fill = NPM1_lab)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, width = 0.6) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.7) +
  scale_fill_manual(values = lab_cols, guide = "none") +
  labs(title = "A. Venetoclax AUC by NPM1 status",
       subtitle = sprintf("n = %d wt, %d mut\nAdjusted difference %.0f (FDR %.2g)",
                          sum(ven$NPM1_status == "wt"), sum(ven$NPM1_status == "mut"),
                          coef(ven_fit)["NPM1_statusmut"], res2[inhibitor == "Venetoclax", fdr]),
       x = NULL, y = "Venetoclax AUC (lower = more sensitive)") +
  theme_bw(base_size = 10)
pb <- ggplot(ven, aes(Monocyte_like, auc, colour = NPM1_lab)) +
  geom_point(size = 1.3, alpha = 0.8) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.8) +
  scale_colour_manual(values = lab_cols, name = NULL) +
  labs(title = "B. Venetoclax AUC vs monocyte-like score",
       subtitle = sprintf("Spearman rho: wt %.2f, mut %.2f\nLines: linear fit per group",
                          ven_rho[NPM1_status == "wt", spearman_rho],
                          ven_rho[NPM1_status == "mut", spearman_rho]),
       x = "Monocyte-like score (Bottomly et al. 2022; higher = more monocytic)",
       y = "Venetoclax AUC") +
  theme_bw(base_size = 10) + theme(legend.position = "bottom")
save_fig(pa + pb + plot_layout(widths = c(1, 1.6)), "09_venetoclax.png", w = 10, h = 5)

# ---- 5. Menin inhibitors? ----------------------------------------------------
# Search every drug name in the AUC file (all samples, not just ours) for
# menin inhibitors by name or code, and for DOT1L inhibitors, the other drug
# class aimed at the same KMT2A/HOXA9 machinery.
menin_patterns <- c(
  menin = "menin|revumenib|SNDX|ziftomenib|KO-539|bleximenib|JNJ-7527|enzomenib|DSP-5336|MI-[0-9]|VTP|BMF-219|icovamenib",
  DOT1L = "DOT1L|pinometostat|EPZ")
menin_check <- data.frame(
  target = names(menin_patterns), pattern = unname(menin_patterns),
  drugs_searched = uniqueN(auc$inhibitor),
  matches = vapply(menin_patterns, function(p)
    paste(grep(p, unique(auc$inhibitor), ignore.case = TRUE, value = TRUE), collapse = "; "), ""))
menin_check$n_matches <- vapply(menin_patterns, function(p)
  length(grep(p, unique(auc$inhibitor), ignore.case = TRUE)), 0L)
write.csv(menin_check, "results/drug_menin_check.csv", row.names = FALSE)
cat("\nMenin / DOT1L inhibitor search:\n"); print(menin_check[, c("target", "drugs_searched", "n_matches")])

# ---- 6. NEK10/FTO/USP44-high vs -low within NPM1-mutated ---------------------
# Report only. Primary: auc ~ block_group (difference in mean AUC, high minus
# low), BH across drugs, same >= 20 per group rule. Two extra columns, not
# used for the FDR: adjusted for FLT3-ITD (more common in the high group, and
# sorafenib, sunitinib and quizartinib all hit FLT3) and adjusted for the
# monocyte-like score (the groups differ strongly in it, so this mostly asks
# whether the difference is "just" maturity). Wilcoxon p kept for comparison
# with Mer et al., who used Wilcoxon.
dm <- d[NPM1_status == "mut"]
n_blk <- counts_by(dm, "block_group")
setnames(n_blk, c("low", "high"), c("n_low", "n_high"))
keep6 <- n_blk[n_low >= min_per_group & n_high >= min_per_group, inhibitor]
cat("\nStep 6: NPM1-mutated samples in the analysis set:",
    paste(names(table(an$block_group)), table(an$block_group), collapse = ", "),
    "| drugs tested:", nrow(n_blk), "| kept:", length(keep6),
    "| skipped:", nrow(n_blk) - length(keep6), "\n")
res6 <- rbindlist(lapply(keep6, function(dr) {
  dd <- dm[inhibitor == dr]
  a <- fit_drug(dd, auc ~ block_group, "block_grouphigh")
  f <- fit_drug(dd, auc ~ block_group + FLT3_ITD, "block_grouphigh")
  m <- fit_drug(dd, auc ~ block_group + Monocyte_like, "block_grouphigh")
  data.table(inhibitor = dr,
             median_auc_low = median(dd$auc[dd$block_group == "low"]),
             median_auc_high = median(dd$auc[dd$block_group == "high"]),
             auc_diff = a$estimate, ci_low = a$ci_low, ci_high = a$ci_high, p = a$p,
             wilcox_p = wilcox.test(auc ~ block_group, data = dd, exact = FALSE)$p.value,
             auc_diff_adj_flt3 = f$estimate, p_adj_flt3 = f$p,
             auc_diff_adj_monocyte = m$estimate, p_adj_monocyte = m$p)
}))
res6[, fdr := p.adjust(p, method = "BH")]
res6 <- merge(n_blk, res6, by = "inhibitor")
res6 <- merge(res6, drugs[, .(inhibitor, drug_class)], by = "inhibitor")
setorder(res6, p)
fwrite(res6, "results/drug_subgroup_npm1.csv")
cat("Drugs at FDR <", fdr_cut, ":", sum(res6$fdr < fdr_cut), "\n")
print(head(res6[, .(inhibitor, n_low, n_high, auc_diff, p, fdr, wilcox_p,
                    auc_diff_adj_flt3, p_adj_flt3, auc_diff_adj_monocyte)], 15), digits = 3)
focus <- c("Sorafenib", "Sunitinib", "Quizartinib (AC220)", "Ruxolitinib (INCB018424)")
cat("\nMer et al. drugs (ruxolitinib was significant in their UHN cohort only):\n")
print(res6[inhibitor %in% focus, .(inhibitor, n_low, n_high, median_auc_low,
                                   median_auc_high, auc_diff, p, fdr, wilcox_p,
                                   auc_diff_adj_flt3, p_adj_flt3,
                                   auc_diff_adj_monocyte, p_adj_monocyte)], digits = 3)
# FLT3-ITD x group counts for the three drugs, to see where the samples sit
cat("\nFLT3-ITD by group, NPM1-mutated samples tested with sorafenib:\n")
print(table(dm[inhibitor == "Sorafenib", .(block_group, FLT3_ITD)]))

kin <- dm[inhibitor %in% focus[1:3]]
kin[, drug := factor(sub(" \\(.*\\)$", "", inhibitor), levels = c("Sorafenib", "Sunitinib", "Quizartinib"))]
lab6 <- res6[inhibitor %in% focus[1:3]]
lab6[, drug := factor(sub(" \\(.*\\)$", "", inhibitor), levels = levels(kin$drug))]
lab6[, txt := sprintf("diff %.0f, FDR %.2g\nFLT3-ITD-adj. p %.2g", auc_diff, fdr, p_adj_flt3)]
p6 <- ggplot(kin, aes(block_group, auc)) +
  geom_boxplot(aes(fill = block_group), outlier.shape = NA, alpha = 0.5, width = 0.6) +
  geom_jitter(aes(shape = FLT3_ITD), width = 0.15, size = 1.6, alpha = 0.8) +
  geom_text(data = lab6, aes(x = 1.5, y = Inf, label = txt), vjust = 1.2, size = 2.8) +
  facet_wrap(~ drug) +
  scale_fill_manual(values = block_cols, guide = "none") +
  scale_shape_manual(values = c(negative = 1, positive = 17), name = "FLT3-ITD") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.3))) +
  labs(title = "Kinase inhibitors in NEK10/FTO/USP44-low vs -high NPM1-mutated AML",
       subtitle = "High = primitive-like (Mer et al. 2021). Lower AUC = more sensitive",
       x = "NEK10 / FTO / USP44 group (three-gene split)", y = "AUC") +
  theme_bw(base_size = 10) + theme(legend.position = "bottom")
save_fig(p6, "09_subgroup_kinase_inhibitors.png", w = 8, h = 4.8)

# ---- 7. HOX program score ----------------------------------------------------
# Blind VST of the 276 main-cohort samples, same low-count rule as before
# (>= 10 reads in >= 100 samples). HOX score = mean VST of the four genes
# (roughly log2 scale). For each drug: auc ~ NPM1_status + HOX_score; the HOX
# coefficient is the AUC change per 1 VST unit of HOX score at the same NPM1
# status. BH across drugs. Extra column: also adjusted for the monocyte-like
# score. Same drug set rule as step 2.
raw <- fread(counts_file, check.names = FALSE)
cm <- as.matrix(raw[, cohort$dbgap_rnaseq_sample, with = FALSE])
rownames(cm) <- raw$stable_id
storage.mode(cm) <- "integer"
keep <- rowSums(cm >= 10) >= 100
vst_mat <- assay(vst(DESeqDataSetFromMatrix(cm[keep, ], colData = cohort, design = ~ 1),
                     blind = TRUE))
sym <- setNames(raw$display_label, raw$stable_id)
hox_ids <- names(sym)[sym %in% hox_genes & names(sym) %in% rownames(vst_mat)]
stopifnot(length(hox_ids) == length(hox_genes))
hox <- data.frame(dbgap_rnaseq_sample = colnames(vst_mat),
                  HOX_score = colMeans(vst_mat[hox_ids, ]))
d <- merge(d, as.data.table(hox), by = "dbgap_rnaseq_sample")
an <- left_join(an, hox, by = "dbgap_rnaseq_sample")
hox_desc <- an %>% group_by(NPM1_status) %>%
  summarise(n = n(), median_HOX = median(HOX_score), sd_HOX = sd(HOX_score), .groups = "drop")
cat("\nHOX score in the analysis set:\n"); print(hox_desc)
cat("Correlation with NPM1 status (point-biserial r):",
    round(cor(an$HOX_score, as.integer(an$NPM1_status == "mut")), 2),
    "| with monocyte-like (Spearman):",
    round(cor(an$HOX_score, an$Monocyte_like, method = "spearman"), 2), "\n")

res7 <- rbindlist(lapply(keep2, function(dr) {
  dd <- d[inhibitor == dr]
  h <- fit_drug(dd, auc ~ NPM1_status + HOX_score, "HOX_score")
  n <- fit_drug(dd, auc ~ NPM1_status + HOX_score, "NPM1_statusmut")
  hm <- fit_drug(dd, auc ~ NPM1_status + HOX_score + Monocyte_like, "HOX_score")
  data.table(inhibitor = dr, n = nrow(dd),
             hox_coef = h$estimate, ci_low = h$ci_low, ci_high = h$ci_high, p = h$p,
             npm1_coef_with_hox = n$estimate, npm1_p_with_hox = n$p,
             hox_coef_adj_monocyte = hm$estimate, p_adj_monocyte = hm$p)
}))
res7[, fdr := p.adjust(p, method = "BH")]
res7[, fdr_adj_monocyte := p.adjust(p_adj_monocyte, method = "BH")]
res7 <- merge(res7, drugs[, .(inhibitor, drug_class)], by = "inhibitor")
setorder(res7, p)
fwrite(res7, "results/drug_hox_score.csv")
cat("Step 7: drugs where HOX score tracks AUC beyond NPM1 status, FDR <", fdr_cut, ":",
    sum(res7$fdr < fdr_cut), "(also with monocyte-like:",
    sum(res7$fdr_adj_monocyte < fdr_cut), ")\n")
print(head(res7[, .(inhibitor, n, hox_coef, p, fdr, npm1_coef_with_hox,
                    hox_coef_adj_monocyte, fdr_adj_monocyte)], 15), digits = 3)

top7 <- head(res7[fdr < fdr_cut], 6)$inhibitor
if (length(top7) > 0) {
  hp <- d[inhibitor %in% top7]
  hp[, drug := factor(sub(" \\(.*\\)$", "", inhibitor),
                      levels = sub(" \\(.*\\)$", "", top7))]
  hp[, NPM1_lab := factor(ifelse(NPM1_status == "mut", "NPM1-mutated", "NPM1 wild-type"),
                          levels = c("NPM1 wild-type", "NPM1-mutated"))]
  p7 <- ggplot(hp, aes(HOX_score, auc, colour = NPM1_lab)) +
    geom_point(size = 0.9, alpha = 0.7) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.7) +
    facet_wrap(~ drug, scales = "free_y") +
    scale_colour_manual(values = lab_cols, name = NULL) +
    labs(title = "Drugs whose AUC tracks the HOX score beyond NPM1 status",
         subtitle = "HOX score = mean VST of HOXA9, HOXA10, MEIS1, PBX3; top drugs by p (FDR < 0.05)",
         x = "HOX score (VST)", y = "AUC (lower = more sensitive)") +
    theme_bw(base_size = 10) + theme(legend.position = "bottom")
  save_fig(p7, "09_hox_score_drugs.png", w = 9, h = 6)
}

fwrite(an %>% select(dbgap_rnaseq_sample, NPM1, block_group, FLT3_ITD,
                     Monocyte_like, HSC_like, HOX_score),
       "results/drug_analysis_samples.csv")

writeLines(capture.output(sessionInfo()), "results/sessionInfo_09_drug_response.txt")
cat("\nSaved results/sessionInfo_09_drug_response.txt\nDone.\n")
