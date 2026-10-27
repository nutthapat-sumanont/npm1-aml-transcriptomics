# 12_build_site.R
# Session 11: data files for the project web page in docs/ (GitHub Pages).
# Session 12: the gene lookup also carries the TCGA-LAML NPM1-mutated vs
# wild-type result (results/tcga_deseq2_npm1.csv), matched on Ensembl gene_id.
# The page itself (docs/*.html, docs/assets/) is hand-written plain HTML, CSS
# and JavaScript; this script only rebuilds what the page reads.
#
#  1. Per-gene expression summary for the gene lookup: quartiles and box-plot
#     whiskers of VST expression in NPM1-mutated, KMT2A-rearranged and other
#     AML. Same samples and VST as the Session 4 three-group plot
#     (05_enrichment.R step 6). Saved to results/site_gene_group_quartiles.csv.
#     This is the only step that needs data/ (the counts file); if the CSV
#     already exists it is reused, unless you set REBUILD_QUARTILES=1.
#  2. docs/data/genes.json: one row per gene with the three group summaries,
#     the main DESeq2 result (results/deseq2_npm1.csv), the TCGA-LAML result
#     (results/tcga_deseq2_npm1.csv, null = not tested in TCGA), plus short
#     notes for genes discussed in LEARNING.md.
#  3. Copy the key figures from figures/ into docs/figures/.
#  4. Check docs/ stays under 20 MB.
#
# Run from the repo root:  Rscript scripts/12_build_site.R
# Needs (step 1 only, if the CSV is missing): Rscript scripts/01_download.R

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(jsonlite)
})
if (!file.exists("CLAUDE.md")) stop("Run this from the repo root")
set.seed(2026)  # nothing here is random; set for consistency with the other scripts

quart_file  <- "results/site_gene_group_quartiles.csv"
counts_file <- "data/beataml_waves1to4_counts_dbgap.txt"
groups <- c("NPM1-mutated", "KMT2A-rearranged", "Other AML")

# ---- 1. Per-gene group summaries ---------------------------------------------
# Box-plot numbers per gene and group: quartiles (R's default quantile type,
# as in results/three_group_summary.csv) and Tukey whiskers (the most extreme
# sample still within 1.5 x IQR of the box). Individual samples are not
# stored, so the page loads fast and no per-patient values are published.

box_stats <- function(x) {
  q <- quantile(x, c(0.25, 0.5, 0.75), names = FALSE)
  iqr <- q[3] - q[1]
  c(n = length(x), min = min(x), lower = min(x[x >= q[1] - 1.5 * iqr]),
    q1 = q[1], median = q[2], q3 = q[3],
    upper = max(x[x <= q[3] + 1.5 * iqr]), max = max(x))
}

if (!file.exists(quart_file) || Sys.getenv("REBUILD_QUARTILES") == "1") {
  suppressPackageStartupMessages(library(DESeq2))
  if (!file.exists(counts_file)) stop("Missing ", counts_file, "; run scripts/01_download.R first")

  # Same cohort, filter and VST as 05_enrichment.R step 6
  cohort <- read.csv("results/cohort.csv") %>%
    filter(main_cohort | group == "KMT2A_r") %>%
    mutate(group3 = factor(case_when(group == "NPM1_mut" ~ "NPM1-mutated",
                                     group == "KMT2A_r"  ~ "KMT2A-rearranged",
                                     TRUE                ~ "Other AML"),
                           levels = groups))
  rownames(cohort) <- cohort$dbgap_rnaseq_sample
  cat("Samples:\n"); print(table(cohort$group3))

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
  cat("Genes after the low-count filter:", nrow(vst_mat), "\n")

  sym <- setNames(raw$display_label, raw$stable_id)
  quart <- rbindlist(lapply(groups, function(g) {
    m <- vst_mat[, cohort$group3 == g, drop = FALSE]
    s <- t(apply(m, 1, box_stats))
    data.table(gene_id = rownames(m), symbol = unname(sym[rownames(m)]),
               group = g, s)
  }))

  # Check against the Session 4 table: HOXA9 and MEIS1 must match exactly
  old <- read.csv("results/three_group_summary.csv")
  chk <- merge(old, quart, by.x = c("gene", "group3"), by.y = c("symbol", "group"))
  bad <- with(chk, abs(median_vst - median) > 1e-6 | abs(q25 - q1) > 1e-6 |
                   abs(q75 - q3) > 1e-6 | n.x != n.y)
  if (nrow(chk) != 6 || any(bad)) stop("Quartiles differ from results/three_group_summary.csv")
  cat("HOXA9 and MEIS1 quartiles match results/three_group_summary.csv\n")
  num <- c("min", "lower", "q1", "median", "q3", "upper", "max")
  quart[, (num) := lapply(.SD, round, 3), .SDcols = num]
  fwrite(quart, quart_file)
} else {
  cat("Reusing", quart_file, "(set REBUILD_QUARTILES=1 to recompute from data/)\n")
}

quart <- fread(quart_file)

# ---- 2. docs/data/genes.json -------------------------------------------------
de <- fread("results/deseq2_npm1.csv")
tcga <- fread("results/tcga_deseq2_npm1.csv")  # Session 7, ~ sex + NPM1_status, 111 samples
if (anyDuplicated(tcga$gene_id)) stop("Duplicate gene_id in results/tcga_deseq2_npm1.csv")
sex_removed <- fread("results/enrichment_removed_sex_genes.csv")$gene_id

# Short notes for genes discussed in LEARNING.md (wording from the log)
pos_up <- c("HOXA9", "HOXA10", "MEIS1", "PBX3", paste0("HOXB", 2:9))
notes <- c(
  setNames(rep("Positive control: expected higher in NPM1-mutated AML (HOX program).",
               length(pos_up)), pos_up),
  CD34  = "Positive control: expected lower in NPM1-mutated AML, which is usually CD34-negative.",
  FTO   = paste("FTO locus / IFEX9. In TCGA-LAML coverage the rise in the NEK10/FTO/USP44-high",
                "group is IFEX9, a transcript starting in FTO intron 8 and spliced onto FTO exon 9;",
                "FTO exons 1 to 8 are not higher. Beat AML's gene-level FTO count includes IFEX9",
                "(listed as FTO-007 in Ensembl 75)."),
  NEK10 = paste("Part of the NEK10/FTO/USP44 module: up mainly in the primitive-like",
                "(NEK10/FTO/USP44-high) subgroup of NPM1-mutated AML, close to zero in the committed-like one."),
  USP44 = paste("Part of the NEK10/FTO/USP44 module: up mainly in the primitive-like",
                "(NEK10/FTO/USP44-high) subgroup of NPM1-mutated AML, close to zero in the committed-like one."),
  WDR86 = "Tracks FLT3-ITD rather than the NEK10/FTO/USP44 module, so it was dropped from the subgroup score.",
  RPS4Y1 = paste("Male-specific chrY gene. Its 'down' call is an artefact of one female wild-type",
                 "sample (BA3206R) with partial Y expression; removed before enrichment."),
  KDM5D = paste("Male-specific chrY gene. Its 'down' call is an artefact of one female wild-type",
                "sample (BA3206R) with partial Y expression; removed before enrichment."),
  XIST  = paste("XIST is on in female cells, so its expression mostly shows sex (64% of NPM1-mutated",
                "vs 41% of wild-type patients are female). Removed from the ranked list before enrichment.")
)

ids <- union(de$gene_id, quart$gene_id)
sym <- c(setNames(de$symbol, de$gene_id), setNames(quart$symbol, quart$gene_id))
sym <- sym[!duplicated(names(sym))]
biotype <- setNames(de$biotype, de$gene_id)

# Expression: 3 groups x 5 numbers (lower whisker, q1, median, q3, upper whisker)
expr_cols <- c("lower", "q1", "median", "q3", "upper")
expr <- lapply(split(quart, by = "gene_id"), function(d) {
  d <- d[match(groups, d$group)]
  round(as.vector(t(as.matrix(d[, ..expr_cols]))), 2)
})
n_by_group <- quart[, .(n = n[1]), by = group][match(groups, group), n]

de_idx <- match(ids, de$gene_id)
# Same Ensembl stable ID in both annotations (Ensembl 75 and GENCODE v36, version
# suffix already removed in 08_tcga_replication.R). A gene with no TCGA row either
# has no matching ID in GENCODE v36 or failed TCGA's low-count filter
# (>= 10 reads in >= 31 samples); the page calls both "not tested in TCGA".
tcga_idx <- match(ids, tcga$gene_id)
rows <- lapply(seq_along(ids), function(i) {
  id <- ids[i]; j <- de_idx[i]; s <- unname(sym[id])
  d <- if (is.na(j)) NULL else list(
    signif(de$log2FC[j], 3), signif(de$lfcSE[j], 3),
    if (is.na(de$padj[j])) NULL else signif(de$padj[j], 3),
    de$rank[j], as.integer(isTRUE(de$sig[j])))
  k <- tcga_idx[i]
  tc <- if (is.na(k)) NULL else list(
    signif(tcga$log2FC[k], 3), signif(tcga$lfcSE[k], 3),
    if (is.na(tcga$padj[k])) NULL else signif(tcga$padj[k], 3),
    tcga$rank[k], as.integer(isTRUE(tcga$sig[k])),
    if (identical(tcga$symbol[k], s)) NULL else tcga$symbol[k])  # renamed since 2014
  note <- if (id %in% sex_removed && !(s %in% names(notes))) {
    paste("Male-specific chrY gene: removed from the ranked list before gene-set enrichment, because",
          "a sex-adjusted model still handles Y genes with near-zero counts in females badly.")
  } else if (s %in% names(notes)) unname(notes[s]) else NULL
  list(s, id, if (is.na(j)) NULL else biotype[[id]], expr[[id]], d, note, tc)
})

meta <- list(
  groups = groups, n = n_by_group,
  expr_fields = expr_cols,
  de_fields = c("log2FC", "lfcSE", "padj", "rank", "sig"),
  de_n_genes = nrow(de),
  tcga_fields = c("log2FC", "lfcSE", "padj", "rank", "sig", "symbol_if_different"),
  tcga_n_genes = nrow(tcga),
  built = format(Sys.Date()),
  expression = paste("DESeq2 VST (roughly log2), blind, on the 276 main-cohort samples plus",
                     "24 KMT2A-rearranged cases; genes with at least 10 reads in at least 100 samples."),
  deseq2 = paste("~ wave + specimen_type + sex + NPM1_status on the 276 main-cohort samples,",
                 "NPM1 wild-type as reference; apeglm-shrunken log2FC; significant = padj < 0.05",
                 "and |log2FC| > 1."),
  tcga = paste("TCGA-LAML (GDC STAR counts, GENCODE v36): ~ sex + NPM1_status on 111 main-cohort",
               "samples (31 NPM1-mutated, 80 wild-type), same cohort filters as Beat AML; apeglm-shrunken",
               "log2FC; significant = padj < 0.05 and |log2FC| > 1.")
)
dir.create("docs/data", recursive = TRUE, showWarnings = FALSE)
write_json(list(meta = meta, genes = rows), "docs/data/genes.json",
           auto_unbox = TRUE, null = "null", digits = NA)
cat("docs/data/genes.json:", length(rows), "genes,",
    round(file.size("docs/data/genes.json") / 1024^2, 2), "MB\n")
cat("  with expression:", sum(ids %in% quart$gene_id),
    " with a DESeq2 result:", sum(!is.na(de_idx)),
    " with a TCGA result:", sum(!is.na(tcga_idx)), "\n")

# ---- 3. Figures ----------------------------------------------------------------
site_figs <- c(
  "05_three_group_hoxa9_meis1.png", "04_volcano.png", "05_gsea_verhaak_up.png",
  "08_tcga_lfc_vs_beataml.png", "07_celltype_adjusted_lfc.png",
  "05_block_nek10_fto_usp44_wdr86.png", "06_subgroup_shared_vs_specific.png",
  "11_fto_coverage.png", "09_venetoclax.png", "10_km_npm1_flt3.png"
)
dir.create("docs/figures", showWarnings = FALSE)
old_figs <- setdiff(list.files("docs/figures"), site_figs)
if (length(old_figs)) file.remove(file.path("docs/figures", old_figs))
ok <- file.copy(file.path("figures", site_figs), "docs/figures", overwrite = TRUE)
if (!all(ok)) stop("Could not copy: ", paste(site_figs[!ok], collapse = ", "))
cat("Copied", length(site_figs), "figures to docs/figures/\n")

# ---- 4. Size check ---------------------------------------------------------------
docs_mb <- sum(file.size(list.files("docs", recursive = TRUE, full.names = TRUE))) / 1024^2
cat(sprintf("docs/ total: %.1f MB\n", docs_mb))
if (docs_mb > 20) stop("docs/ is over 20 MB")

writeLines(capture.output(sessionInfo()), "results/sessionInfo_12_build_site.txt")
