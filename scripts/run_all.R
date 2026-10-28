# run_all.R
# Session 12: runs the whole pipeline from scratch, in order, from the repo root:
#   Rscript scripts/run_all.R
# Each script runs in its own fresh R process (like typing Rscript by hand), so
# nothing left over from one step can leak into the next. Stops at the first
# script that fails, because every later script reads what the earlier ones
# wrote. The time each script took goes to results/run_all_timings.csv.
#
# To pick up again after a failure without redoing the earlier steps:
#   Rscript scripts/run_all.R --from 08
#
# Took 37.5 minutes from scratch on a 4-core cloud machine in Session 12, about
# 60% of it in 08 and 11, which download TCGA-LAML data (GDC, recount3).
# Session 13 added 13_classifier.R. It runs before 12, out of number order,
# because 12 copies its ROC figure into docs/figures/.
# Session 14 added 14_paper_figures.R (paper figures and supplementary
# workbook from results/ only), last because it reads what 04 to 13 wrote.

scripts <- c(
  "scripts/01_download.R",           # Beat AML 2.0 files into data/
  "scripts/02_cohort.R",             # cohort.csv: 276 main-cohort patients
  "scripts/03_qc.R",                 # filtering, PCA, sex check
  "scripts/04_deseq2.R",             # main DESeq2 model + FLT3-ITD sensitivity
  "scripts/05_enrichment.R",         # GSEA, three-group plot, NEK10/FTO/USP44 group
  "scripts/06_subgroup.R",           # subgroup vs published NPM1 subtypes
  "scripts/06b_deseq2_subgroups.R",  # sensitivity DESeq2 with NPM1-high/low
  "scripts/07_cell_composition.R",   # cell-type scores, maturity sensitivity
  "scripts/08_tcga_replication.R",   # replication in TCGA-LAML (downloads ~650 MB)
  "scripts/09_drug_response.R",      # ex vivo drug response (AUC) by NPM1 status
  "scripts/10_survival.R",           # exploratory overall survival (Beat AML, TCGA)
  "scripts/11_fto_locus.R",          # FTO vs IFEX9 in TCGA coverage (recount3)
  "scripts/13_classifier.R",         # LASSO NPM1 classifier, nested CV, locked test in TCGA
  "scripts/12_build_site.R",         # data files for the web page in docs/ (copies 13's figure)
  "scripts/14_paper_figures.R"       # paper figures (figures/paper/) and supplementary tables (results/paper/)
)

if (!file.exists("CLAUDE.md")) stop("Run this from the repo root")

# Guard against a new numbered script that was never added to the list above
on_disk <- file.path("scripts", list.files("scripts", pattern = "^[0-9]+[a-z]?_.*\\.R$"))
missing <- setdiff(on_disk, scripts)
if (length(missing)) stop("Not in run_all.R: ", paste(missing, collapse = ", "))
gone <- scripts[!file.exists(scripts)]
if (length(gone)) stop("Listed in run_all.R but not found: ", paste(gone, collapse = ", "))

args <- commandArgs(trailingOnly = TRUE)
from <- if (length(args) >= 2 && args[1] == "--from") args[2] else NULL
if (!is.null(from)) {
  start <- grep(paste0("^scripts/", from), scripts)
  if (length(start) != 1) stop("--from ", from, " doesn't match exactly one script")
  scripts <- scripts[start:length(scripts)]
}

rscript <- file.path(R.home("bin"), "Rscript")
log_file <- "results/run_all_timings.csv"
dir.create("results", showWarnings = FALSE)
timings <- data.frame(script = character(), started = character(),
                      minutes = numeric(), status = character())
t_all <- Sys.time()

for (s in scripts) {
  t0 <- Sys.time()
  cat("\n==== ", s, " (", format(t0, "%H:%M:%S"), ") ====\n", sep = "")
  status <- system2(rscript, s)
  mins <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2)
  timings[nrow(timings) + 1, ] <- list(s, format(t0, "%Y-%m-%d %H:%M:%S"), mins,
                                       if (status == 0) "ok" else paste("failed, exit", status))
  write.csv(timings, log_file, row.names = FALSE)
  cat(sprintf("---- %s: %s in %.1f min\n", s, timings$status[nrow(timings)], mins))
  if (status != 0) stop(s, " failed with exit status ", status,
                        ". Fix it, then resume with: Rscript scripts/run_all.R --from ",
                        sub("^scripts/([0-9]+[a-z]?)_.*", "\\1", s))
}

total <- as.numeric(difftime(Sys.time(), t_all, units = "mins"))
cat(sprintf("\nAll %d scripts finished in %.1f minutes. Timings: %s\n",
            length(scripts), total, log_file))
print(timings[, c("script", "minutes")], row.names = FALSE)
