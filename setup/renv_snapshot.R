# renv_snapshot.R
# Session 12: writes renv.lock, the exact version of every R package the
# pipeline uses (and everything those packages need), so the analysis can be
# rebuilt later with the same versions:  Rscript -e 'renv::restore()'
#
# Run from the repo root after setup/cloud_setup.sh:  Rscript setup/renv_snapshot.R
# To only compare what's installed with renv.lock, without writing anything:
#   Rscript setup/renv_snapshot.R --check
#
# The package list is written out by hand rather than found by renv's code
# scan, because the scan misses apeglm (DESeq2 only calls it through
# lfcShrink(type = "apeglm")). renv adds every dependency of these packages.
# This doesn't switch the project to an renv library: the scripts keep using
# the normal R library that cloud_setup.sh fills.

if (!file.exists("CLAUDE.md")) stop("Run this from the repo root")

pkgs <- c(
  # Bioconductor
  "DESeq2", "apeglm", "fgsea", "BiocParallel", "GenomicRanges",
  "BiocFileCache", "recount3", "rtracklayer",
  # CRAN
  "msigdbr", "data.table", "dplyr", "tidyr", "readxl", "jsonlite", "httr",
  "ggplot2", "ggrepel", "pheatmap", "patchwork", "scales", "matrixStats",
  "pROC", "survival", "glmnet", "openxlsx",
  # the tools that restore the rest
  "renv", "BiocManager"
)

if ("--check" %in% commandArgs(trailingOnly = TRUE)) {
  locked <- vapply(renv::lockfile_read("renv.lock")$Packages, `[[`, "", "Version")
  # Compare as version numbers: packageVersion() writes "3.8-12" as "3.8.12"
  installed <- vapply(names(locked), function(p) tryCatch(as.character(packageVersion(p)),
                                                          error = function(e) NA_character_), "")
  diff <- vapply(names(locked), function(p) is.na(installed[[p]]) ||
                   package_version(installed[[p]]) != package_version(locked[[p]]), logical(1))
  if (!any(diff)) {
    cat("All", length(locked), "packages match renv.lock\n")
  } else {
    cat(sum(diff), "of", length(locked), "packages differ from renv.lock:\n")
    print(data.frame(locked = locked[diff], installed = installed[diff]))
  }
  quit(save = "no")
}

missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Not installed: ", paste(missing, collapse = ", "),
                          ". Run setup/cloud_setup.sh first.")

renv::snapshot(project = ".", library = .libPaths(), packages = pkgs,
               lockfile = "renv.lock", prompt = FALSE)
lock <- renv::lockfile_read("renv.lock")
cat("renv.lock:", length(lock$Packages), "packages, R", lock$R$Version,
    if (!is.null(lock$Bioconductor)) paste(", Bioconductor", lock$Bioconductor$Version), "\n")
