# run_all.R
# Runs the whole pipeline from scratch, in order, from the repo root:
#   Rscript scripts/run_all.R
# Each script runs in its own fresh R process (like typing Rscript by hand), so
# nothing left over from one step can leak into the next. Stops at the first
# script that fails. Takes roughly 30 minutes plus download time.

scripts <- c(
  "scripts/01_download.R",          # Beat AML 2.0 files into data/
  "scripts/02_cohort.R",            # cohort.csv: 276 main-cohort patients
  "scripts/03_qc.R",                # filtering, PCA, sex check
  "scripts/04_deseq2.R",            # main DESeq2 model + FLT3-ITD sensitivity
  "scripts/05_enrichment.R",        # GSEA, three-group plot, NEK10/FTO block
  "scripts/06_subgroup.R",          # block group vs published NPM1 subgroups
  "scripts/06b_deseq2_subgroups.R",  # sensitivity DESeq2 with NPM1-high/low
  "scripts/08_tcga_replication.R"    # replication in TCGA-LAML (downloads ~650 MB)
)

if (!file.exists("CLAUDE.md")) stop("Run this from the repo root")
rscript <- file.path(R.home("bin"), "Rscript")

for (s in scripts) {
  cat("\n==== ", s, " (", format(Sys.time(), "%H:%M:%S"), ") ====\n", sep = "")
  status <- system2(rscript, s)
  if (status != 0) stop(s, " failed with exit status ", status)
}
cat("\nAll", length(scripts), "scripts finished.\n")
