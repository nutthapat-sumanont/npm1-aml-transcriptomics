# 03_qc.R
# Quality control of the RNA-seq counts for the main cohort (277 patients,
# 100 NPM1-mutated, 177 wild-type) before any differential expression.
#
#  1. Build a DESeq2 dataset, with a hard check that the count matrix columns
#     and the cohort rows are in exactly the same order.
#  2. Drop genes with very low counts.
#  3. Library size per sample, with outliers flagged.
#  4. VST, then PCA coloured by NPM1 status, specimen type and wave.
#  5. Sex check: XIST and RPS4Y1 against the recorded sex.
#  6. HOXA9, MEIS1 and CD34 by NPM1 status as a first sanity check.
#
# Run from the repo root:  Rscript scripts/03_qc.R
# Needs data/ from scripts/01_download.R and results/cohort.csv from 02_cohort.R.

suppressPackageStartupMessages({
  library(DESeq2)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
})

set.seed(2026)
dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

counts_file <- "data/beataml_waves1to4_counts_dbgap.txt"
npm1_cols   <- c(wt = "grey55", mut = "#D55E00")

save_fig <- function(p, file, w = 7, h = 5) {
  ggsave(file.path("figures", file), p, width = w, height = h, dpi = 300)
  cat("Saved figures/", file, "\n", sep = "")
}

# ---- 1. DESeq2 dataset -------------------------------------------------------

# QC looks at the cohort as it was before any QC exclusion: the 276 main-cohort
# samples plus BA2294R, which this QC flagged (sex mismatch) and 02_cohort.R
# step 9 then excluded. Keeping it here means a rerun still shows the evidence
# for dropping it (Session 12; before that, a rerun silently used 276).
cohort <- read.csv("results/cohort.csv") %>%
  filter(main_cohort | group == "QC_excluded") %>%
  mutate(
    NPM1_status   = factor(NPM1, levels = c("wt", "mut")),   # wt = reference
    # Leukapheresis is blood collected by apheresis, so it counts as blood
    specimen_type = factor(ifelse(specimenType == "Bone Marrow Aspirate",
                                  "BM", "Blood"), levels = c("BM", "Blood")),
    specimen_detail = factor(specimenType),
    wave          = factor(cohort, levels = c("Waves1+2", "Waves3+4", "Both")),
    sex           = factor(consensus_sex)
  )
rownames(cohort) <- cohort$dbgap_rnaseq_sample
cat("QC cohort (main cohort + QC-excluded):", nrow(cohort), "samples\n")
print(table(cohort$NPM1_status))

raw <- fread(counts_file, check.names = FALSE)
gene_info <- as.data.frame(raw[, .(stable_id, display_label, biotype)])
missing <- setdiff(cohort$dbgap_rnaseq_sample, names(raw))
if (length(missing) > 0) stop("Samples missing from count matrix: ",
                              paste(missing, collapse = ", "))

# Pull the cohort samples in cohort order, then check the order explicitly
counts <- as.matrix(raw[, cohort$dbgap_rnaseq_sample, with = FALSE])
rownames(counts) <- gene_info$stable_id
if (!identical(colnames(counts), rownames(cohort))) {
  stop("Sample order in the count matrix does not match the cohort table")
}
if (any(counts < 0) || any(counts != round(counts))) {
  stop("Counts must be non-negative whole numbers")
}
storage.mode(counts) <- "integer"
cat("Count matrix:", nrow(counts), "genes x", ncol(counts), "samples;",
    "sample order matches cohort table\n")

# Planned design from CLAUDE.md. It is stored now but not fitted in this script;
# all QC below is blind to the design.
dds <- DESeqDataSetFromMatrix(counts, colData = cohort,
                              design = ~ wave + specimen_type + NPM1_status)
rowData(dds)$symbol  <- gene_info$display_label
rowData(dds)$biotype <- gene_info$biotype

# ---- 2. Low-count filter -----------------------------------------------------
# Keep a gene if it has at least 10 reads in at least 100 samples. 100 is the
# size of the smaller group (NPM1-mutated), so a gene switched on only in the
# mutated group can still pass. This is the rule suggested in the DESeq2 guide.

min_count   <- 10
min_samples <- min(table(cohort$NPM1_status))
keep <- rowSums(counts(dds) >= min_count) >= min_samples
cat("\nLow-count filter: >=", min_count, "reads in >=", min_samples, "samples\n")
cat("Genes before:", nrow(dds), "| kept:", sum(keep), "| dropped:", sum(!keep), "\n")
cat("Kept genes by biotype (top 5):\n")
print(head(sort(table(rowData(dds)$biotype[keep]), decreasing = TRUE), 5))
dds <- dds[keep, ]

controls <- c("HOXA9", "HOXA10", "MEIS1", "PBX3", "CD34", "XIST", "RPS4Y1")
cat("Control genes kept:", paste(controls[controls %in% rowData(dds)$symbol],
                                 collapse = ", "), "\n")

# ---- 3. Library size ---------------------------------------------------------
# Library size = total reads assigned to genes (all genes, before filtering).
# Outlier rule: more than 3 MADs from the median on the log10 scale.

lib <- cohort %>%
  select(sample = dbgap_rnaseq_sample, NPM1_status, specimen_type, wave) %>%
  mutate(lib_size = colSums(counts), log_lib = log10(lib_size),
         size_factor = estimateSizeFactors(dds)$sizeFactor)
med <- median(lib$log_lib); md <- mad(lib$log_lib)
lib$outlier <- abs(lib$log_lib - med) > 3 * md
cat("\nLibrary size (millions): median", round(median(lib$lib_size) / 1e6, 1),
    "| range", round(min(lib$lib_size) / 1e6, 1), "to",
    round(max(lib$lib_size) / 1e6, 1), "\n")
cat("Library size outliers (> 3 MAD on log10 scale):", sum(lib$outlier), "\n")
if (any(lib$outlier)) print(lib[lib$outlier, c("sample", "lib_size", "NPM1_status", "wave")])
write.csv(lib, "results/qc_library_sizes.csv", row.names = FALSE)

lib_plot <- lib %>% arrange(lib_size) %>% mutate(rank = row_number())
p <- ggplot(lib_plot, aes(rank, lib_size / 1e6, fill = NPM1_status)) +
  geom_col(width = 1) +
  geom_hline(yintercept = 10^(med + c(-3, 3) * md) / 1e6, linetype = "dashed") +
  geom_text_repel(data = filter(lib_plot, outlier),
                  aes(label = sample), size = 3, min.segment.length = 0) +
  scale_fill_manual(values = npm1_cols) +
  labs(x = "Samples, sorted by library size", y = "Reads assigned to genes (millions)",
       fill = "NPM1", title = "Library size per sample",
       subtitle = "Dashed lines: median +/- 3 MAD (log10 scale)") +
  theme_bw()
save_fig(p, "03_library_size.png", w = 8, h = 4.5)

p <- ggplot(lib, aes(wave, lib_size / 1e6, colour = NPM1_status)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2), size = 1, alpha = 0.6) +
  scale_colour_manual(values = npm1_cols) +
  labs(x = "Wave", y = "Reads assigned to genes (millions)", colour = "NPM1",
       title = "Library size by wave and NPM1 status") +
  theme_bw()
save_fig(p, "03_library_size_by_wave.png", w = 6, h = 4.5)

# ---- 4. VST and PCA ----------------------------------------------------------
# blind = TRUE: the transform ignores the design, which is what we want for QC.

vsd <- vst(dds, blind = TRUE)
pca <- plotPCA(vsd, intgroup = c("NPM1_status", "specimen_detail", "wave"),
               ntop = 500, returnData = TRUE)
pv  <- round(100 * attr(pca, "percentVar"), 1)
# PCs 3 and 4 too, from the same 500 most variable genes
top <- head(order(matrixStats::rowVars(assay(vsd)), decreasing = TRUE), 500)
pc  <- prcomp(t(assay(vsd)[top, ]))
pv4 <- round(100 * pc$sdev^2 / sum(pc$sdev^2), 1)
pcs <- data.frame(sample = colnames(vsd), PC1 = pc$x[, 1], PC2 = pc$x[, 2],
                  PC3 = pc$x[, 3], PC4 = pc$x[, 4], as.data.frame(colData(vsd))[
                    , c("NPM1_status", "specimen_detail", "specimen_type", "wave", "sex")])
cat("\nPCA on top 500 variable genes, % variance PC1-PC4:",
    paste(pv4[1:4], collapse = ", "), "\n")

# Samples far from the rest: more than 4 SD from the centre on PC1 or PC2
pcs$pca_outlier <- abs(as.vector(scale(pcs$PC1))) > 4 | abs(as.vector(scale(pcs$PC2))) > 4
cat("PCA outliers (> 4 SD on PC1 or PC2):", sum(pcs$pca_outlier), "\n")
if (any(pcs$pca_outlier)) print(pcs[pcs$pca_outlier, c("sample", "PC1", "PC2", "NPM1_status")])

# How much of each PC goes with each variable (R squared from a simple lm)
r2 <- sapply(c("NPM1_status", "specimen_type", "wave", "sex"), function(v)
  sapply(paste0("PC", 1:4), function(p) summary(lm(pcs[[p]] ~ pcs[[v]]))$r.squared))
cat("R squared of each PC against each variable:\n"); print(round(r2, 3))
write.csv(pcs, "results/qc_pca_scores.csv", row.names = FALSE)
write.csv(round(r2, 4), "results/qc_pca_r2.csv")

pca_plot <- function(colour_by, values, title, x = "PC1", y = "PC2") {
  ggplot(pcs, aes(.data[[x]], .data[[y]], colour = .data[[colour_by]])) +
    geom_point(size = 1.8, alpha = 0.8) +
    geom_text_repel(data = filter(pcs, pca_outlier), aes(label = sample),
                    size = 3, colour = "black") +
    scale_colour_manual(values = values) +
    labs(x = paste0(x, " (", pv4[as.integer(sub("PC", "", x))], "%)"),
         y = paste0(y, " (", pv4[as.integer(sub("PC", "", y))], "%)"),
         colour = NULL, title = title, subtitle = "VST, top 500 variable genes") +
    coord_fixed() + theme_bw()
}
spec_cols <- c("Bone Marrow Aspirate" = "#0072B2", "Peripheral Blood" = "#CC79A7",
               "Leukapheresis" = "#009E73")
wave_cols <- c("Waves1+2" = "#E69F00", "Waves3+4" = "#56B4E9", "Both" = "black")
save_fig(pca_plot("NPM1_status", npm1_cols, "PCA by NPM1 status"), "03_pca_npm1.png")
save_fig(pca_plot("specimen_detail", spec_cols, "PCA by specimen type"), "03_pca_specimen.png")
save_fig(pca_plot("wave", wave_cols, "PCA by wave"), "03_pca_wave.png")
save_fig(pca_plot("NPM1_status", npm1_cols, "PCA by NPM1 status, PC3 vs PC4",
                  x = "PC3", y = "PC4"), "03_pca_npm1_pc3_pc4.png")

# ---- 5. Sex check ------------------------------------------------------------
# XIST is made only from the inactive X, so it is high in females. RPS4Y1 sits
# on the Y chromosome, so it is high in males. Predicted sex comes from the gap
# between the two (RPS4Y1 minus XIST on the VST scale). Note: AML cells can
# lose the Y chromosome, which makes a true male look "female" on RPS4Y1.

sym <- rowData(vsd)$symbol
get_gene <- function(g) assay(vsd)[which(sym == g)[1], ]
sexdf <- data.frame(sample = colnames(vsd), recorded = as.character(colData(vsd)$sex),
                    XIST = get_gene("XIST"), RPS4Y1 = get_gene("RPS4Y1"))
sexdf$diff <- sexdf$RPS4Y1 - sexdf$XIST
# Two clusters on the difference, seeded from the extremes
km <- kmeans(sexdf$diff, centers = c(min(sexdf$diff), max(sexdf$diff)))
sexdf$predicted <- ifelse(km$cluster == which.max(km$centers), "Male", "Female")
sexdf$mismatch  <- sexdf$predicted != sexdf$recorded
cat("\nSex check, recorded vs predicted:\n")
print(table(recorded = sexdf$recorded, predicted = sexdf$predicted))
mm <- sexdf[sexdf$mismatch, ]
cat("Mismatches:", nrow(mm), "\n")
if (nrow(mm) > 0) print(mm %>% mutate(across(c(XIST, RPS4Y1, diff), ~ round(.x, 2))),
                        row.names = FALSE)
write.csv(sexdf, "results/qc_sex_check.csv", row.names = FALSE)

p <- ggplot(sexdf, aes(XIST, RPS4Y1, colour = recorded)) +
  geom_point(size = 1.8, alpha = 0.8) +
  geom_point(data = mm, shape = 21, size = 4, colour = "black", stroke = 0.8) +
  geom_text_repel(data = mm, aes(label = sample), size = 3, colour = "black") +
  scale_colour_manual(values = c(Female = "#CC79A7", Male = "#0072B2")) +
  labs(x = "XIST (VST)", y = "RPS4Y1 (VST)", colour = "Recorded sex",
       title = "Sex check", subtitle = "Circled: expression disagrees with recorded sex") +
  theme_bw()
save_fig(p, "03_sex_check.png", w = 6, h = 5)

# ---- 6. Sanity check: HOXA9, MEIS1, CD34 -------------------------------------
# Expected: HOXA9 and MEIS1 higher in NPM1-mutated, CD34 lower.

genes <- c("HOXA9", "MEIS1", "CD34")
long <- bind_rows(lapply(genes, function(g) data.frame(
  gene = g, NPM1_status = colData(vsd)$NPM1_status, expr = get_gene(g))))
long$gene <- factor(long$gene, levels = genes)
sanity <- long %>% group_by(gene) %>% summarise(
  median_wt  = median(expr[NPM1_status == "wt"]),
  median_mut = median(expr[NPM1_status == "mut"]),
  diff_mut_minus_wt = median_mut - median_wt,
  wilcox_p = wilcox.test(expr ~ NPM1_status)$p.value)
cat("\nSanity genes, median VST by group (VST is roughly log2):\n")
print(as.data.frame(sanity), digits = 3)
write.csv(sanity, "results/qc_sanity_genes.csv", row.names = FALSE)

p <- ggplot(long, aes(NPM1_status, expr, fill = NPM1_status)) +
  geom_violin(alpha = 0.4, colour = NA) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  geom_jitter(width = 0.12, size = 0.6, alpha = 0.5) +
  facet_wrap(~ gene, scales = "free_y") +
  scale_fill_manual(values = npm1_cols, guide = "none") +
  labs(x = "NPM1 status", y = "Expression (VST, ~log2)",
       title = "Positive controls: HOXA9 and MEIS1 up, CD34 down in NPM1-mutated") +
  theme_bw()
save_fig(p, "03_sanity_genes.png", w = 8, h = 4)

# ---- Save --------------------------------------------------------------------
# dds and vsd go to data/ (gitignored); re-create them by re-running this script
saveRDS(dds, "data/dds_filtered.rds")
saveRDS(vsd, "data/vsd_blind.rds")
writeLines(capture.output(sessionInfo()), "results/sessionInfo_03_qc.txt")
cat("\nDone.\n")
