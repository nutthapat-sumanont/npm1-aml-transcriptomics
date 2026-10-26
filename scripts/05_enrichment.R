# 05_enrichment.R
# Gene-set enrichment for NPM1-mutated vs NPM1-wild-type AML, plus two
# follow-up checks.
#
#  1. Remove sex-chromosome artefacts from the ranked list: genes in the
#     male-specific part of chrY, plus XIST (approved after Session 3).
#  2. Rank the remaining genes by the DESeq2 Wald statistic and run preranked
#     GSEA (fgsea) against MSigDB Hallmark and GO Biological Process.
#  3. Positive controls: VERHAAK_AML_WITH_NPM1_MUTATED_UP and _DN (MSigDB C2),
#     with an enrichment plot for each.
#  4. Dot plot of the top 15 up and top 15 down pathways by NES.
#  5. Three-group plot of HOXA9 and MEIS1 (VST): NPM1-mutated,
#     KMT2A-rearranged, other AML.
#  6. Does the NEK10 / FTO / USP44 / WDR86 block mark a subgroup within
#     NPM1-mutated samples? Report only, no filtering.
#
# Run from the repo root:  Rscript scripts/05_enrichment.R   (about 5 minutes)
# Needs results/deseq2_npm1.csv (04_deseq2.R), results/cohort.csv (02_cohort.R)
# and data/ (01_download.R, for the VST in step 5 and 6).

suppressPackageStartupMessages({
  library(fgsea)
  library(msigdbr)
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
chrY_file   <- "results/chrY_genes_grch37.csv"
min_size    <- 15     # smallest gene set tested (genes present in our list)
max_size    <- 500    # largest gene set tested
padj_cut    <- 0.05
n_top       <- 15
group_cols  <- c("NPM1-mutated" = "#D55E00", "KMT2A-rearranged" = "#0072B2",
                 "Other AML" = "grey55")

save_fig <- function(p, file, w = 7, h = 5) {
  ggsave(file.path("figures", file), p, width = w, height = h, dpi = 300)
  cat("Saved figures/", file, "\n", sep = "")
}

# fgsea returns leadingEdge as a list column; turn it into gene symbols in one
# string so the table can be written to CSV.
tidy_fgsea <- function(res, id2sym) {
  res <- as.data.frame(res)
  res$leadingEdge <- vapply(res$leadingEdge,
                            function(g) paste(id2sym[g], collapse = ";"), "")
  res %>% arrange(padj, desc(abs(NES)))
}

# ---- 1. Remove chrY genes and XIST -------------------------------------------
# Beat AML counts use Ensembl 75 (GRCh37) gene IDs. The DESeq2 table has no
# chromosome column, so we ask the Ensembl GRCh37 REST server for every gene on
# chrY once and keep that list in results/ (small, committed), so later runs
# don't need the internet.
#
# Only the male-specific region of Y (MSY) is removed. The two pseudoautosomal
# regions (PAR1, PAR2) are identical on X and Y and their genes (CD99, CRLF2,
# IL3RA, ...) are expressed from both chromosomes in both sexes, so they are not
# sex artefacts. GRCh37 PAR coordinates on Y are from the GRC.

par1 <- c(10001, 2649520)
par2 <- c(59034050, 59363566)
chrY_len <- 59373566

if (!file.exists(chrY_file)) {
  cat("Fetching chrY genes from grch37.rest.ensembl.org ...\n")
  starts <- seq(1, chrY_len, by = 5e6)   # the server takes at most 5 Mb a call
  chrY <- do.call(rbind, lapply(starts, function(s) {
    u <- sprintf(paste0("https://grch37.rest.ensembl.org/overlap/region/human/",
                        "Y:%d-%d?feature=gene;content-type=application/json"),
                 s, min(s + 5e6 - 1, chrY_len))
    d <- fromJSON(u)
    if (length(d) == 0) NULL else
      d[, c("id", "external_name", "biotype", "start", "end")]
  }))
  chrY <- unique(chrY)
  chrY$region <- ifelse(chrY$end <= par1[2], "PAR1",
                 ifelse(chrY$start >= par2[1], "PAR2", "MSY"))
  write.csv(chrY, chrY_file, row.names = FALSE)
}
chrY <- read.csv(chrY_file)

de <- read.csv("results/deseq2_npm1.csv")
sex_genes <- de %>%
  filter(gene_id %in% chrY$id[chrY$region == "MSY"] | symbol == "XIST") %>%
  mutate(reason = ifelse(symbol == "XIST", "XIST", "chrY (male-specific)"))
par_kept <- de %>% filter(gene_id %in% chrY$id[chrY$region != "MSY"])

cat("\nGenes in the DESeq2 table:", nrow(de), "\n")
cat("Removed before enrichment:", nrow(sex_genes),
    "(", sum(sex_genes$reason != "XIST"), "male-specific chrY +",
    sum(sex_genes$reason == "XIST"), "XIST )\n")
cat("  of these, significant in Session 3 (padj < 0.05, |log2FC| > 1):",
    sum(sex_genes$sig), "->",
    paste(sex_genes$symbol[sex_genes$sig], collapse = ", "), "\n")
cat("Pseudoautosomal (PAR) genes kept:", nrow(par_kept), "\n")
write.csv(sex_genes[, c("gene_id", "symbol", "biotype", "reason", "log2FC",
                        "stat", "padj", "rank", "sig")],
          "results/enrichment_removed_sex_genes.csv", row.names = FALSE)

de <- de %>% filter(!gene_id %in% sex_genes$gene_id)
cat("Genes left for ranking:", nrow(de), "\n")

# ---- 2. Ranked list and gene sets --------------------------------------------
# The ranking is the Wald statistic (log2FC / its standard error, unshrunken).
# It carries both the direction and how sure we are, so a gene with a big but
# noisy change doesn't jump to the top. Genes are matched to gene sets by
# Ensembl ID: more of them match than by symbol (our symbols are from 2014)
# and there are no duplicate IDs.

ranks  <- sort(setNames(de$stat, de$gene_id), decreasing = TRUE)
id2sym <- setNames(de$symbol, de$gene_id)

get_sets <- function(...) {
  m <- msigdbr(species = "Homo sapiens", ...)
  split(m$ensembl_gene, m$gs_name) %>% lapply(unique)
}
sets <- list(
  hallmark = get_sets(collection = "H"),
  gobp     = get_sets(collection = "C5", subcollection = "GO:BP"),
  verhaak  = get_sets(collection = "C2", subcollection = "CGP")[
    c("VERHAAK_AML_WITH_NPM1_MUTATED_UP", "VERHAAK_AML_WITH_NPM1_MUTATED_DN")]
)
msig_version <- unique(msigdbr(species = "Homo sapiens", collection = "H")$db_version)
cat("\nMSigDB version:", msig_version, "\n")
for (nm in names(sets)) cat(sprintf("  %-8s %5d gene sets\n", nm, length(sets[[nm]])))

# ---- 3. GSEA -----------------------------------------------------------------
# fgseaMultilevel estimates very small p-values accurately (eps = 0 means no
# lower bound). Each collection is tested on its own, so the FDR (padj) is
# corrected within Hallmark and within GO BP separately.

run_gsea <- function(pw) {
  set.seed(2026)
  fgseaMultilevel(pw, ranks, minSize = min_size, maxSize = max_size, eps = 0,
                  BPPARAM = BiocParallel::SerialParam(progressbar = FALSE))
}
gsea <- lapply(sets, run_gsea)
gsea_tab <- lapply(gsea, tidy_fgsea, id2sym = id2sym)

write.csv(gsea_tab$hallmark, "results/gsea_hallmark.csv", row.names = FALSE)
write.csv(gsea_tab$gobp,     "results/gsea_gobp.csv",     row.names = FALSE)
write.csv(gsea_tab$verhaak,  "results/gsea_verhaak_controls.csv", row.names = FALSE)

gsea_summary <- do.call(rbind, lapply(names(gsea_tab), function(nm) {
  t <- gsea_tab[[nm]]
  data.frame(collection = nm, sets_tested = nrow(t),
             padj_lt_0.05 = sum(t$padj < padj_cut, na.rm = TRUE),
             up = sum(t$padj < padj_cut & t$NES > 0, na.rm = TRUE),
             down = sum(t$padj < padj_cut & t$NES < 0, na.rm = TRUE))
}))
write.csv(gsea_summary, "results/gsea_summary.csv", row.names = FALSE)
cat("\nGSEA summary (padj < 0.05):\n"); print(gsea_summary, row.names = FALSE)

show_top <- function(t, n = 10) {
  t <- t[!is.na(t$padj), ]
  up <- t %>% filter(NES > 0) %>% arrange(desc(NES)) %>% head(n)
  dn <- t %>% filter(NES < 0) %>% arrange(NES) %>% head(n)
  print(rbind(up, dn)[, c("pathway", "size", "NES", "padj")], row.names = FALSE)
}
cat("\nHallmark, top 10 each way by NES:\n"); show_top(gsea_tab$hallmark)

# ---- 4. Positive controls: Verhaak NPM1 sets ---------------------------------
# Verhaak et al. 2005 (Haematologica) compared NPM1-mutated with wild-type AML
# on microarrays. Their UP set should sit at the top of our ranking (NES > 0)
# and their DN set at the bottom (NES < 0). If not, something is wrong upstream.

vh <- gsea_tab$verhaak
cat("\nPositive controls:\n")
print(vh[, c("pathway", "size", "NES", "pval", "padj")], row.names = FALSE)
up_ok <- with(vh[vh$pathway == "VERHAAK_AML_WITH_NPM1_MUTATED_UP", ], NES > 0 & padj < padj_cut)
dn_ok <- with(vh[vh$pathway == "VERHAAK_AML_WITH_NPM1_MUTATED_DN", ], NES < 0 & padj < padj_cut)
if (!isTRUE(up_ok) || !isTRUE(dn_ok)) {
  stop("Verhaak positive controls are not enriched in the expected direction; ",
       "look for a bug before going further")
}
cat("Both Verhaak positive controls pass.\n")

for (pw in vh$pathway) {
  r <- vh[vh$pathway == pw, ]
  p <- plotEnrichment(sets$verhaak[[pw]], ranks, ticksSize = 0.3) +
    labs(title = pw,
         subtitle = sprintf("NES = %.2f, FDR = %.1e, %d genes in our list",
                            r$NES, r$padj, r$size),
         x = "Gene rank (1 = most up in NPM1-mutated, by Wald statistic)",
         y = "Running enrichment score") +
    theme_bw(base_size = 11)
  save_fig(p, paste0("05_gsea_", tolower(sub(".*MUTATED_", "verhaak_", pw)), ".png"),
           w = 7, h = 4.5)
}

# ---- 5. Dot plot: top 15 up and top 15 down ----------------------------------
# The 15 pathways with the highest NES and the 15 with the lowest. One plot
# per collection, since their FDRs are separate. Not pre-filtered on FDR,
# because Hallmark has no "up" set below 0.05 and an empty half would hide
# that; instead filled dots are FDR < 0.05 and open dots are not.

dot_plot <- function(t, prefix, title) {
  t <- t %>% filter(!is.na(padj))
  top <- rbind(t %>% filter(NES > 0) %>% arrange(desc(NES)) %>% head(n_top),
               t %>% filter(NES < 0) %>% arrange(NES) %>% head(n_top))
  n_pos <- sum(top$NES > 0)
  top$label <- sub(paste0("^", prefix, "_"), "", top$pathway)
  top$label <- gsub("_", " ", top$label)
  top$label <- ifelse(nchar(top$label) > 60,
                      paste0(substr(top$label, 1, 57), "..."), top$label)
  top$label <- factor(top$label, levels = top$label[order(top$NES)])
  top$neglog10_padj <- -log10(pmax(top$padj, 1e-50))
  top$fdr_ok <- factor(ifelse(top$padj < padj_cut, "< 0.05", ">= 0.05"),
                       levels = c("< 0.05", ">= 0.05"))
  ggplot(top, aes(NES, label, size = size, colour = neglog10_padj, shape = fdr_ok)) +
    geom_vline(xintercept = 0, colour = "grey70") +
    geom_point(stroke = 1) +
    scale_shape_manual(values = c("< 0.05" = 16, ">= 0.05" = 1), name = "FDR",
                       drop = FALSE) +
    scale_colour_viridis_c(name = "-log10 FDR") +
    scale_size_area(name = "Genes in set", max_size = 6) +
    labs(title = title,
         subtitle = paste0("Top ", n_top, " up and down by NES. ",
                           "NES > 0 = higher in NPM1-mutated",
                           if (n_pos == 0) "\nNo set in this collection has NES > 0" else ""),
         x = "Normalized enrichment score (NES)", y = NULL) +
    theme_bw(base_size = 10)
}
save_fig(dot_plot(gsea_tab$hallmark, "HALLMARK", "MSigDB Hallmark"),
         "05_gsea_dotplot_hallmark.png", w = 8, h = 6.5)
save_fig(dot_plot(gsea_tab$gobp, "GOBP", "GO Biological Process"),
         "05_gsea_dotplot_gobp.png", w = 10, h = 8)

# ---- 6. VST for the three-group plot and the block check ---------------------
# All main-cohort samples plus the KMT2A-rearranged cases (same filters
# otherwise: AML, initial diagnosis, de novo, not APL). VST is blind to group,
# with the same low-count filter as before (>= 10 reads in >= 100 samples).

cohort <- read.csv("results/cohort.csv") %>%
  filter(main_cohort | group == "KMT2A_r") %>%
  mutate(group3 = factor(case_when(group == "NPM1_mut" ~ "NPM1-mutated",
                                   group == "KMT2A_r"  ~ "KMT2A-rearranged",
                                   TRUE                ~ "Other AML"),
                         levels = names(group_cols)))
rownames(cohort) <- cohort$dbgap_rnaseq_sample
cat("\nSamples for VST:\n"); print(table(cohort$group3))

raw <- fread(counts_file, check.names = FALSE)
counts <- as.matrix(raw[, cohort$dbgap_rnaseq_sample, with = FALSE])
rownames(counts) <- raw$stable_id
if (!identical(colnames(counts), rownames(cohort))) {
  stop("Sample order in the count matrix does not match the cohort table")
}
storage.mode(counts) <- "integer"
keep <- rowSums(counts >= 10) >= 100
dds  <- DESeqDataSetFromMatrix(counts[keep, ], colData = cohort, design = ~ 1)
vst_mat <- assay(vst(dds, blind = TRUE))
sym <- setNames(raw$display_label, raw$stable_id)
gene_vst <- function(s) vst_mat[names(sym)[sym == s & names(sym) %in% rownames(vst_mat)], ]

# ---- 7. Three-group plot: HOXA9 and MEIS1 ------------------------------------
# KMT2A fusions and NPM1 mutations both keep menin/KMT2A on the HOXA/MEIS1
# genes, which is why both are menin-inhibitor targets. If that's right, both
# groups should be high and "other AML" lower.

three <- do.call(rbind, lapply(c("HOXA9", "MEIS1"), function(g) {
  data.frame(sample = colnames(vst_mat), gene = g, vst = gene_vst(g),
             group3 = cohort$group3, NPM1 = cohort$NPM1)
}))
write.csv(three, "results/three_group_hoxa9_meis1.csv", row.names = FALSE)

three_summary <- three %>% group_by(gene, group3) %>%
  summarise(n = n(), median_vst = median(vst),
            q25 = quantile(vst, 0.25), q75 = quantile(vst, 0.75),
            .groups = "drop")
pairs <- list(c("NPM1-mutated", "Other AML"), c("KMT2A-rearranged", "Other AML"),
              c("NPM1-mutated", "KMT2A-rearranged"))
three_tests <- do.call(rbind, lapply(c("HOXA9", "MEIS1"), function(g) {
  do.call(rbind, lapply(pairs, function(pr) {
    a <- three$vst[three$gene == g & three$group3 == pr[1]]
    b <- three$vst[three$gene == g & three$group3 == pr[2]]
    data.frame(gene = g, comparison = paste(pr, collapse = " vs "),
               median_diff = median(a) - median(b),
               wilcox_p = wilcox.test(a, b)$p.value)
  }))
}))
write.csv(three_summary, "results/three_group_summary.csv", row.names = FALSE)
write.csv(three_tests, "results/three_group_tests.csv", row.names = FALSE)
cat("\nThree-group medians (VST):\n"); print(as.data.frame(three_summary))
cat("\nWilcoxon tests:\n"); print(three_tests, row.names = FALSE)

n_lab <- three %>% filter(gene == "HOXA9") %>% count(group3) %>%
  mutate(lab = paste0(group3, "\n(n = ", n, ")"))
p <- ggplot(three, aes(group3, vst, fill = group3)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, width = 0.6) +
  geom_jitter(aes(shape = ifelse(NPM1 == "mut" & group3 == "KMT2A-rearranged",
                                 "yes", "no")),
              width = 0.15, size = 1.1, alpha = 0.7) +
  scale_shape_manual(values = c(no = 16, yes = 17), breaks = "yes",
                     labels = "KMT2A-r that also carries an NPM1 mutation",
                     name = NULL) +
  facet_wrap(~ gene, scales = "free_y") +
  scale_fill_manual(values = group_cols, guide = "none") +
  scale_x_discrete(labels = setNames(n_lab$lab, n_lab$group3)) +
  labs(title = "HOXA9 and MEIS1 in the two menin-dependent AML groups",
       subtitle = "De novo AML at diagnosis, Beat AML 2.0",
       x = NULL, y = "Expression (VST, ~log2)") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
save_fig(p, "05_three_group_hoxa9_meis1.png", w = 8, h = 5)

# ---- 8. NEK10 / FTO / USP44 / WDR86 block within NPM1-mutated ----------------
# All four are among the top "up" genes. If they rise and fall together inside
# the NPM1-mutated group, they might mark a subgroup. Block score = mean z-score
# of the four genes across the 100 NPM1-mutated samples. Then we ask whether
# the score lines up with anything in cohort.csv. Report only.

block_genes <- c("NEK10", "FTO", "USP44", "WDR86")
mut <- cohort %>% filter(group == "NPM1_mut") %>%
  mutate(specimen_type = ifelse(specimenType == "Bone Marrow Aspirate", "BM", "Blood"))
bm <- t(sapply(block_genes, function(g) gene_vst(g)[mut$dbgap_rnaseq_sample]))

cor_mut <- cor(t(bm), method = "spearman")
cat("\nSpearman correlation of the block genes within NPM1-mutated (n = 100):\n")
print(round(cor_mut, 2))
cor_wt <- cor(sapply(block_genes, function(g)
  gene_vst(g)[cohort$dbgap_rnaseq_sample[cohort$group == "NPM1_wt"]]),
  method = "spearman")
cat("Same in NPM1-wild-type (n = 176):\n"); print(round(cor_wt, 2))

mut$block_score <- colMeans(t(scale(t(bm))))
# Is there a natural split? Two-cluster k-means on the score (seeded).
set.seed(2026)
km <- kmeans(mut$block_score, centers = 2, nstart = 25)
hi <- which.max(km$centers)
mut$block_group <- factor(ifelse(km$cluster == hi, "high", "low"), levels = c("low", "high"))
cat("\nk-means split of the block score:", paste(names(table(mut$block_group)),
    table(mut$block_group), collapse = ", "), "\n")
# Bimodality coefficient (> 0.555 hints at two peaks; a rough guide only)
sk <- mean(scale(mut$block_score)^3); ku <- mean(scale(mut$block_score)^4) - 3
n <- nrow(mut)
bc <- (sk^2 + 1) / (ku + 3 * (n - 1)^2 / ((n - 2) * (n - 3)))
cat("Bimodality coefficient:", round(bc, 3), "\n")

# Score vs each cohort variable. Categorical: Kruskal-Wallis on the score and
# Fisher's test on the high/low split. Age: Spearman correlation.
cat_vars <- c("FLT3_ITD", "consensus_sex", "specimen_type", "cohort",
              "centerID", "ELN2017")
assoc <- do.call(rbind, lapply(cat_vars, function(v) {
  x <- factor(mut[[v]])
  ok <- !is.na(x)
  med <- tapply(mut$block_score[ok], x[ok], median)
  data.frame(variable = v, test = "Kruskal-Wallis (score); Fisher (high/low)",
             levels = paste(sprintf("%s: n=%d, median %.2f", names(med),
                                    as.integer(table(x[ok])), med), collapse = " | "),
             p_score = kruskal.test(mut$block_score[ok], x[ok])$p.value,
             p_split = fisher.test(table(x[ok], mut$block_group[ok]),
                                   simulate.p.value = nlevels(droplevels(x[ok])) > 2,
                                   B = 1e4)$p.value)
}))
ct <- cor.test(mut$block_score, mut$ageAtDiagnosis, method = "spearman", exact = FALSE)
assoc <- rbind(assoc, data.frame(
  variable = "ageAtDiagnosis", test = "Spearman correlation",
  levels = sprintf("rho = %.2f", ct$estimate), p_score = ct$p.value,
  p_split = wilcox.test(ageAtDiagnosis ~ block_group, data = mut)$p.value))
write.csv(assoc, "results/block_nek10_fto_associations.csv", row.names = FALSE)
cat("\nBlock score vs cohort variables (NPM1-mutated only):\n")
print(assoc[, c("variable", "p_score", "p_split")], row.names = FALSE)
cat("\nFLT3-ITD vs block group:\n"); print(table(FLT3_ITD = mut$FLT3_ITD, block = mut$block_group))

# WDR86 correlates only weakly with the other three, so also test each gene on
# its own and a three-gene core score (NEK10, FTO, USP44) against FLT3-ITD.
core <- c("NEK10", "FTO", "USP44")
mut$core_score <- colMeans(t(scale(t(bm[core, ]))))
flt3_pos <- mut$FLT3_ITD == "positive"
per_gene <- do.call(rbind, lapply(c(block_genes, "core_score", "block_score"), function(g) {
  v <- if (g %in% block_genes) bm[g, ] else mut[[g]]
  data.frame(feature = g,
             median_FLT3_ITD_neg = median(v[!flt3_pos]),
             median_FLT3_ITD_pos = median(v[flt3_pos]),
             wilcox_p = wilcox.test(v[flt3_pos], v[!flt3_pos])$p.value)
}))
write.csv(per_gene, "results/block_nek10_fto_flt3_by_gene.csv", row.names = FALSE)
cat("\nEach block gene vs FLT3-ITD within NPM1-mutated:\n")
print(per_gene, row.names = FALSE, digits = 3)

write.csv(cbind(mut[, c("dbgap_rnaseq_sample", "FLT3_ITD", "consensus_sex",
                        "specimen_type", "cohort", "ELN2017", "ageAtDiagnosis",
                        "block_score", "core_score", "block_group")], t(bm)),
          "results/block_nek10_fto_scores.csv", row.names = FALSE)

long <- data.frame(sample = rep(mut$dbgap_rnaseq_sample, each = 4),
                   gene = factor(rep(block_genes, nrow(mut)), levels = block_genes),
                   vst = as.vector(bm),
                   FLT3_ITD = rep(mut$FLT3_ITD, each = 4),
                   score = rep(mut$block_score, each = 4))
long$sample <- factor(long$sample, levels = mut$dbgap_rnaseq_sample[order(mut$block_score)])
p <- ggplot(long, aes(sample, vst, colour = FLT3_ITD)) +
  geom_point(size = 1) +
  facet_wrap(~ gene, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = c(negative = "grey55", positive = "#CC79A7")) +
  labs(title = "NEK10 / FTO / USP44 / WDR86 within NPM1-mutated AML",
       subtitle = "Samples ordered by block score (mean z-score of the four genes)",
       x = "NPM1-mutated samples (n = 100)", y = "VST") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
save_fig(p, "05_block_nek10_fto_usp44_wdr86.png", w = 7, h = 7)

# ---- 9. Session info ---------------------------------------------------------

writeLines(c(paste("MSigDB version:", msig_version), capture.output(sessionInfo())),
           "results/sessionInfo_05_enrichment.txt")
cat("\nSaved results/sessionInfo_05_enrichment.txt\nDone.\n")
