# 06_subgroup.R
# Session 5, step 1: what is the NEK10 / FTO / USP44-high group inside
# NPM1-mutated AML? Report only; nothing is filtered and the DESeq2 model is
# not touched.
#
#  1a. Co-mutations from the Beat AML exome (WES) calls, plus FLT3-ITD from the
#      clinical table, high vs low group. Count samples with NPM1 + FLT3-ITD +
#      (TET2 or IDH2) together.
#  1b. Compare with two published NPM1 subgroups:
#      - Arza-Apalategi et al., Leukemia 2023 (doi:10.1038/s41375-023-01844-w):
#        an intragenic lncRNA inside FTO (IFEX9), high in about a third of
#        NPM1-mutated AML, which gene-level counts read as "high FTO". Linked
#        to FLT3-ITD, TET2, IDH2 and an immature profile.
#      - Mer et al., Nat Commun 2021 (doi:10.1038/s41467-021-21233-0):
#        "primitive" vs "committed" NPM1-mutated AML. We score their marker
#        genes and ask whether our high group lines up with "primitive".
#  1c. Which annotated genes sit inside FTO in our gene annotation (Ensembl 75),
#      to show that gene-level counts can't separate FTO from IFEX9.
#  1d. (Session 8) TET2, IDH1 or IDH2, high vs low. A consistency check, not an
#      independent test: the grouping came from Session 7 (TCGA) after looking
#      at the data, and partly from Beat AML itself (Session 5).
#
# Session 8: the high/low split is now k-means on the three-gene NEK10/FTO/USP44
# score (WDR86 dropped), from the rerun of 05_enrichment.R.
#
# Run from the repo root:  Rscript scripts/06_subgroup.R   (about 2 minutes)
# Needs results/cohort.csv (02), results/block_nek10_fto_scores.csv (05) and
# data/ (01_download.R).

suppressPackageStartupMessages({
  library(DESeq2)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(jsonlite)
})

set.seed(2026)
dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

counts_file <- "data/beataml_waves1to4_counts_dbgap.txt"
wes_file    <- "data/beataml_wes_wv1to4_mutations_dbgap.txt"
fto_file    <- "results/fto_locus_genes_grch37.csv"

save_fig <- function(p, file, w = 7, h = 5) {
  ggsave(file.path("figures", file), p, width = w, height = h, dpi = 300)
  cat("Saved figures/", file, "\n", sep = "")
}

# The high/low split comes from 05_enrichment.R (k-means on the three-gene
# NEK10/FTO/USP44 core score since Session 8).
block <- read.csv("results/block_nek10_fto_scores.csv")
cohort <- read.csv("results/cohort.csv")
mut <- cohort %>% filter(main_cohort, group == "NPM1_mut") %>%
  inner_join(block[, c("dbgap_rnaseq_sample", "core_score", "block_group")],
             by = "dbgap_rnaseq_sample") %>%
  mutate(block_group = factor(block_group, levels = c("low", "high")))
stopifnot(nrow(mut) == 100)
cat("NPM1-mutated samples:", nrow(mut), "\n"); print(table(mut$block_group))

# ---- 1a. Co-mutations --------------------------------------------------------
# WES calls are already filtered to coding, non-silent variants. A sample with
# no exome (no dbgap_dnaseq_sample, or absent from the WES file) is NA, not
# "wild-type". FLT3-ITD isn't called from exomes, so it comes from cohort.csv
# (clinical lab result). FLT3-TKD = missense / in-frame changes at D835 or I836.

wes <- fread(wes_file)
has_wes <- !is.na(mut$dbgap_dnaseq_sample) &
  mut$dbgap_dnaseq_sample %in% wes$dbgap_sample_id
cat("\nSamples with exome calls:", sum(has_wes), "of", nrow(mut), "\n")

mut_genes <- c("NPM1", "DNMT3A", "TET2", "IDH1", "IDH2", "FLT3_TKD",
               "NRAS", "KRAS", "PTPN11", "cohesin")
cohesin <- c("SMC1A", "SMC3", "RAD21", "STAG2")
carrier <- function(g) {
  hit <- switch(g,
    FLT3_TKD = wes$symbol == "FLT3" & grepl("^p\\.(D835|I836)", wes$hgvsp_short),
    cohesin  = wes$symbol %in% cohesin,
    wes$symbol == g)
  x <- mut$dbgap_dnaseq_sample %in% wes$dbgap_sample_id[hit]
  ifelse(has_wes, x, NA)
}
for (g in mut_genes) mut[[paste0("mut_", g)]] <- carrier(g)
mut$FLT3_ITD_pos <- mut$FLT3_ITD == "positive"

# NPM1 is detected in the exome for only some samples (small insertions in a
# hard-to-capture exon). Report it as a sensitivity check of the WES calls.
cat("NPM1 mutation seen in the exome:", sum(mut$mut_NPM1, na.rm = TRUE),
    "of", sum(has_wes), "with exome data\n")

# NPM1 VAF from the exome (rough purity proxy; highest NPM1 call per sample)
npm1_vaf <- wes[symbol == "NPM1", .(npm1_vaf = max(t_vaf)), by = dbgap_sample_id]
mut$npm1_vaf <- npm1_vaf$npm1_vaf[match(mut$dbgap_dnaseq_sample,
                                        npm1_vaf$dbgap_sample_id)]

feat <- c("FLT3_ITD_pos", paste0("mut_", setdiff(mut_genes, "NPM1")))
comut <- do.call(rbind, lapply(feat, function(f) {
  x <- mut[[f]]; ok <- !is.na(x)
  tab <- table(factor(x[ok], levels = c(FALSE, TRUE)), mut$block_group[ok])
  data.frame(feature = sub("^mut_", "", f),
             n_tested = sum(ok),
             low_n  = tab["TRUE", "low"],  low_of  = sum(tab[, "low"]),
             high_n = tab["TRUE", "high"], high_of = sum(tab[, "high"]),
             low_pct  = round(100 * tab["TRUE", "low"]  / sum(tab[, "low"]), 1),
             high_pct = round(100 * tab["TRUE", "high"] / sum(tab[, "high"]), 1),
             odds_ratio = unname(fisher.test(tab)$estimate),
             fisher_p = fisher.test(tab)$p.value)
}))
comut$fisher_padj <- p.adjust(comut$fisher_p, method = "BH")
write.csv(comut, "results/subgroup_comutations.csv", row.names = FALSE)
cat("\nCo-mutations, NEK10/FTO/USP44 low vs high (Fisher, BH across",
    nrow(comut), "features):\n")
print(comut, row.names = FALSE, digits = 3)

vaf_p <- wilcox.test(npm1_vaf ~ block_group, data = mut)$p.value
cat(sprintf("\nNPM1 VAF (exome), median low %.2f vs high %.2f, Wilcoxon p = %.2g\n",
            median(mut$npm1_vaf[mut$block_group == "low"], na.rm = TRUE),
            median(mut$npm1_vaf[mut$block_group == "high"], na.rm = TRUE), vaf_p))

# The triple combination from Arza-Apalategi et al.: NPM1 + FLT3-ITD + (TET2 or IDH2)
mut$tet2_or_idh2 <- ifelse(has_wes, mut$mut_TET2 | mut$mut_IDH2, NA)
mut$triple <- ifelse(has_wes, mut$FLT3_ITD_pos & mut$tet2_or_idh2, NA)
triple_tab <- table(triple = mut$triple, block = mut$block_group, useNA = "ifany")
cat("\nNPM1 + FLT3-ITD + (TET2 or IDH2):\n"); print(triple_tab)
# Denominator = the 97 samples with exome data, so every row uses the same set.
combos <- list(
  "NPM1 + FLT3-ITD"                  = mut$FLT3_ITD_pos & has_wes,
  "NPM1 + (TET2 or IDH2)"            = mut$tet2_or_idh2,
  "NPM1 + FLT3-ITD + TET2"           = mut$FLT3_ITD_pos & mut$mut_TET2,
  "NPM1 + FLT3-ITD + IDH2"           = mut$FLT3_ITD_pos & mut$mut_IDH2,
  "NPM1 + FLT3-ITD + (TET2 or IDH2)" = mut$triple)
lo <- has_wes & mut$block_group == "low"
hi <- has_wes & mut$block_group == "high"
triple <- do.call(rbind, lapply(names(combos), function(nm) {
  x <- combos[[nm]] & has_wes
  data.frame(combination = nm, n_all = sum(x[has_wes]), of_all = sum(has_wes),
             n_low = sum(x[lo]), of_low = sum(lo),
             n_high = sum(x[hi]), of_high = sum(hi),
             fisher_p = fisher.test(table(factor(x[has_wes], c(FALSE, TRUE)),
                                          mut$block_group[has_wes]))$p.value)
}))
write.csv(triple, "results/subgroup_triple_mutation.csv", row.names = FALSE)
print(triple, row.names = FALSE, digits = 3)

# ---- 1d. TET2, IDH1 or IDH2 (Session 8) --------------------------------------
# Mutant IDH1 and IDH2 make 2-hydroxyglutarate, which blocks TET2, so the three
# share one pathway. In Session 7 this grouping stood out in TCGA (all 6
# three-gene-high samples). It is NOT an independent test here: the idea came
# partly from Beat AML (TET2 and IDH2 were higher in the high group in
# Session 5) and was formed after looking at TCGA. One pre-stated test, so no
# multiple-testing correction; read it as "is Beat AML consistent with it?".
mut$tet2_idh1_idh2 <- ifelse(has_wes, mut$mut_TET2 | mut$mut_IDH1 | mut$mut_IDH2, NA)
tii_tab <- table(factor(mut$tet2_idh1_idh2[has_wes], c(FALSE, TRUE)),
                 mut$block_group[has_wes])
tii_fisher <- fisher.test(tii_tab)
tii <- data.frame(
  feature = "TET2, IDH1 or IDH2",
  n_low = tii_tab["TRUE", "low"], of_low = sum(tii_tab[, "low"]),
  n_high = tii_tab["TRUE", "high"], of_high = sum(tii_tab[, "high"]),
  low_pct = round(100 * tii_tab["TRUE", "low"] / sum(tii_tab[, "low"]), 1),
  high_pct = round(100 * tii_tab["TRUE", "high"] / sum(tii_tab[, "high"]), 1),
  odds_ratio = unname(tii_fisher$estimate),
  or_ci_low = tii_fisher$conf.int[1], or_ci_high = tii_fisher$conf.int[2],
  fisher_p = tii_fisher$p.value,
  note = "Consistency check, not independent: hypothesis formed partly from Beat AML")
write.csv(tii, "results/subgroup_tet2_idh1_idh2.csv", row.names = FALSE)
cat("\nTET2, IDH1 or IDH2 (consistency check, not an independent test):\n")
print(tii[, 1:10], row.names = FALSE, digits = 3)

# ---- 1b. Primitive vs committed markers (Mer et al. 2021) --------------------
# Same VST as Session 4 step 5 (300 samples, blind, >= 10 reads in >= 100
# samples), except the marker genes are kept even if they fail the low-count
# filter, so a lowly expressed marker isn't silently dropped.

primitive <- c("CDH2", "GPR12", "MDFI", "ZNF521", "PLXNB1", "HHIP")
committed <- c("CD163", "C1QA", "CD14", "MARCO", "MSR1")
fto_like  <- c("FTO", "NEK10", "USP44")

vst_cohort <- cohort %>% filter(main_cohort | group == "KMT2A_r")
raw <- fread(counts_file, check.names = FALSE)
counts <- as.matrix(raw[, vst_cohort$dbgap_rnaseq_sample, with = FALSE])
rownames(counts) <- raw$stable_id
storage.mode(counts) <- "integer"
sym <- setNames(raw$display_label, raw$stable_id)
markers_id <- names(sym)[sym %in% c(primitive, committed)]
keep <- rowSums(counts >= 10) >= 100
failed <- sym[markers_id[!keep[markers_id]]]
cat("\nMarker genes that fail the low-count filter (kept anyway):",
    if (length(failed)) paste(failed, collapse = ", ") else "none", "\n")
keep[markers_id] <- TRUE
dds <- DESeqDataSetFromMatrix(counts[keep, ], colData = vst_cohort, design = ~ 1)
vst_mat <- assay(vst(dds, blind = TRUE))
gene_vst <- function(s) {
  id <- names(sym)[sym == s & names(sym) %in% rownames(vst_mat)]
  if (length(id) != 1) stop(s, " matches ", length(id), " genes")
  vst_mat[id, mut$dbgap_rnaseq_sample]
}

zmat <- function(genes) t(scale(t(do.call(rbind, lapply(genes, gene_vst)))))
mut$primitive_score <- colMeans(zmat(primitive))
mut$committed_score <- colMeans(zmat(committed))
mut$prim_minus_comm <- mut$primitive_score - mut$committed_score

# Each marker on its own, high vs low
marker_tab <- do.call(rbind, lapply(c(primitive, committed), function(g) {
  v <- gene_vst(g)
  data.frame(gene = g, set = ifelse(g %in% primitive, "primitive", "committed"),
             median_low = median(v[mut$block_group == "low"]),
             median_high = median(v[mut$block_group == "high"]),
             diff_high_minus_low = median(v[mut$block_group == "high"]) -
               median(v[mut$block_group == "low"]),
             wilcox_p = wilcox.test(v ~ mut$block_group)$p.value,
             spearman_vs_core = cor(v, mut$core_score, method = "spearman"))
}))
marker_tab$wilcox_padj <- p.adjust(marker_tab$wilcox_p, method = "BH")
write.csv(marker_tab, "results/subgroup_primitive_committed_markers.csv",
          row.names = FALSE)
cat("\nPrimitive / committed markers, high vs low group:\n")
print(marker_tab, row.names = FALSE, digits = 3)

# Scores vs group and vs the NEK10/FTO/USP44 core score
score_tab <- do.call(rbind, lapply(c("primitive_score", "committed_score",
                                     "prim_minus_comm"), function(s) {
  ct <- cor.test(mut[[s]], mut$core_score, method = "spearman", exact = FALSE)
  data.frame(score = s,
             median_low = median(mut[[s]][mut$block_group == "low"]),
             median_high = median(mut[[s]][mut$block_group == "high"]),
             wilcox_p = wilcox.test(mut[[s]] ~ mut$block_group)$p.value,
             spearman_vs_core = unname(ct$estimate), spearman_p = ct$p.value)
}))
write.csv(score_tab, "results/subgroup_primitive_committed_scores.csv",
          row.names = FALSE)
cat("\nScores (mean z within NPM1-mutated) vs group:\n")
print(score_tab, row.names = FALSE, digits = 3)

# A descriptive primitive / committed call: two-cluster k-means on
# primitive-minus-committed, same idea as the Session 4 block split. Mer et al.
# used their own clustering on many more genes, so this is only an approximation.
km <- kmeans(mut$prim_minus_comm, centers = 2, nstart = 25)
mut$pc_call <- factor(ifelse(km$cluster == which.max(km$centers),
                             "primitive-like", "committed-like"),
                      levels = c("committed-like", "primitive-like"))
pc_tab <- table(block = mut$block_group, call = mut$pc_call)
cat("\nBlock group vs primitive/committed call:\n"); print(pc_tab)
cat("Fisher p =", signif(fisher.test(pc_tab)$p.value, 3), "\n")

# ---- 1c. What sits inside FTO in our annotation? -----------------------------
# The counts are per gene (Ensembl 75). Any read inside FTO's exons that isn't
# assigned to a separately annotated gene counts towards FTO. IFEX9 has no gene
# ID of its own in this annotation, so its reads (if exonic to FTO) end up in
# "FTO" and the two can't be told apart without transcript-level data.

if (!file.exists(fto_file)) {
  cat("\nFetching genes in the FTO locus from grch37.rest.ensembl.org ...\n")
  d <- fromJSON(paste0("https://grch37.rest.ensembl.org/overlap/region/human/",
                       "16:53737875-54155853?feature=gene;content-type=application/json"))
  write.csv(d[, c("id", "external_name", "biotype", "start", "end", "strand")],
            fto_file, row.names = FALSE)
}
fto_locus <- read.csv(fto_file)
fto_locus$in_counts <- fto_locus$id %in% raw$stable_id
fto_locus$passes_filter <- fto_locus$id %in% rownames(vst_mat)
fto_locus$spearman_vs_FTO_in_NPM1mut <- sapply(fto_locus$id, function(id) {
  if (!id %in% rownames(vst_mat) || sym[id] == "FTO") return(NA)
  cor(vst_mat[id, mut$dbgap_rnaseq_sample], gene_vst("FTO"), method = "spearman")
})
fto_locus$median_raw_count_NPM1mut <- sapply(fto_locus$id, function(id)
  if (id %in% rownames(counts)) median(counts[id, mut$dbgap_rnaseq_sample]) else NA)
write.csv(fto_locus, fto_file, row.names = FALSE)
cat("\nGenes annotated inside the FTO locus (GRCh37, Ensembl 75):\n")
print(fto_locus, row.names = FALSE, digits = 2)

# ---- Save per-sample table and figure ----------------------------------------

out <- mut %>% select(dbgap_rnaseq_sample, dbgap_dnaseq_sample, block_group,
                      core_score, FLT3_ITD, starts_with("mut_"), tet2_or_idh2,
                      tet2_idh1_idh2,
                      triple, npm1_vaf, primitive_score, committed_score,
                      prim_minus_comm, pc_call)
write.csv(out, "results/subgroup_per_sample.csv", row.names = FALSE)

long <- rbind(
  data.frame(sample = mut$dbgap_rnaseq_sample, group = mut$block_group,
             what = "Primitive score\n(CDH2, GPR12, MDFI, ZNF521, PLXNB1, HHIP)",
             value = mut$primitive_score),
  data.frame(sample = mut$dbgap_rnaseq_sample, group = mut$block_group,
             what = "Committed score\n(CD163, C1QA, CD14, MARCO, MSR1)",
             value = mut$committed_score))
p_lab <- setNames(sprintf("Wilcoxon p = %.2g", score_tab$wilcox_p[1:2]),
                  unique(long$what))
p <- ggplot(long, aes(group, value, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.6) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.7) +
  facet_wrap(~ what, labeller = labeller(what = function(x) paste0(x, "\n", p_lab[x]))) +
  scale_fill_manual(values = c(low = "grey65", high = "#D55E00"), guide = "none") +
  labs(title = "Mer et al. 2021 markers in the NEK10/FTO/USP44 groups",
       subtitle = "NPM1-mutated AML (n = 100); score = mean z-score of the markers",
       x = "NEK10 / FTO / USP44 group (three-gene split)", y = "Score") +
  theme_bw(base_size = 10)
save_fig(p, "06_subgroup_primitive_committed.png", w = 7, h = 4.5)

writeLines(capture.output(sessionInfo()), "results/sessionInfo_06_subgroup.txt")
cat("\nSaved results/sessionInfo_06_subgroup.txt\nDone.\n")
