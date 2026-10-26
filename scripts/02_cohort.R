# 02_cohort.R
# Builds the analysis cohort: one AML sample per patient, taken at initial
# diagnosis, with RNA-seq counts and a known NPM1 status.
#
# Part 1 only looks at the clinical summary (no filtering).
# Part 2 (added after the filtering plan is approved) applies the filters.
#
# Run from the repo root:  Rscript scripts/02_cohort.R

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(data.table)
})

clin_file   <- "data/beataml_wv1to4_clinical.xlsx"
counts_file <- "data/beataml_waves1to4_counts_dbgap.txt"

# ---- Part 1: what is in the clinical summary? -------------------------------

clin <- read_excel(clin_file, sheet = "summary")
cat("Clinical summary:", nrow(clin), "rows (samples) x", ncol(clin), "columns\n")
cat("Unique patients:", n_distinct(clin$dbgap_subject_id), "\n\n")

# One line per column: type, number of distinct values, missing values,
# and a few example values
for (col in names(clin)) {
  v <- clin[[col]]
  u <- unique(v[!is.na(v)])
  ex <- paste(substr(as.character(head(u, 3)), 1, 40), collapse = " | ")
  cat(sprintf("%-42s %-9s distinct=%-4d NA=%-4d e.g. %s\n",
              col, class(v)[1], length(u), sum(is.na(v)), ex))
}

# The columns that matter for the filters, tabulated in full
cat("\n-- Columns used by the proposed filters --\n")
key_cols <- c("dxAtSpecimenAcquisition", "diseaseStageAtSpecimenCollection",
              "rnaSeq", "NPM1", "specimenType", "isDenovo", "isTransformed")
for (col in key_cols) {
  cat("\n", col, "\n", sep = "")
  print(table(clin[[col]], useNA = "ifany"))
}

# Which RNA-seq sample IDs actually have a column in the counts matrix?
counts_ids <- names(fread(counts_file, nrows = 0))
cat("\nRNA-seq IDs in clinical file:", sum(!is.na(clin$dbgap_rnaseq_sample)),
    "| of those, present in counts matrix:",
    sum(clin$dbgap_rnaseq_sample %in% counts_ids), "\n")
