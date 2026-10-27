# 01_download.R
# Downloads the Beat AML 2.0 files we need into data/.
# Source: https://github.com/biodev/beataml2.0_data (CC-BY-4.0)
# Cite: Bottomly et al., Cancer Cell 2022; Tyner et al., Nature 2018.
#
# The files live in Git LFS. If a download goes wrong you can end up with a
# ~130 byte text "pointer" instead of the real file, so we check every file.
#
# Run from the repo root:  Rscript scripts/01_download.R

options(timeout = 3600)  # the counts file is large; give it time

files <- c(
  raw_counts      = "beataml_waves1to4_counts_dbgap.txt",
  clinical        = "beataml_wv1to4_clinical.xlsx",
  sample_mapping  = "beataml_waves1to4_sample_mapping.xlsx",
  mutations       = "beataml_wes_wv1to4_mutations_dbgap.txt",
  drug_auc        = "beataml_probit_curve_fits_v4_dbgap.txt",
  norm_expression = "beataml_waves1to4_norm_exp_dbgap.txt",
  # Malignant cell-type scores (Bottomly et al. 2022, van Galen et al. 2019
  # signatures), added in Session 6 for scripts/07_cell_composition.R
  celltype_scores = "beataml2_manuscript_vg_cts.xlsx",
  # Drug-to-family table (which inhibitor class each drug belongs to), added
  # in Session 8 for scripts/09_drug_response.R
  drug_families   = "beataml_drug_families.xlsx"
)

# github.com/raw redirects to the real LFS file. media.githubusercontent.com
# is where that redirect lands, and we use it as a fallback in case
# github.com itself is blocked (some cloud sandboxes do this).
# A few small files (the cell-type scores) are plain git files, not LFS, so
# the LFS media URL gives a 404 for them; raw.githubusercontent.com serves those.
urls_for <- function(f) c(
  paste0("https://github.com/biodev/beataml2.0_data/raw/main/", f),
  paste0("https://media.githubusercontent.com/media/biodev/beataml2.0_data/main/", f),
  paste0("https://raw.githubusercontent.com/biodev/beataml2.0_data/main/", f)
)

is_lfs_pointer <- function(path) {
  if (file.size(path) > 1024) return(FALSE)
  first <- readLines(path, n = 1, warn = FALSE)
  length(first) == 1 && grepl("git-lfs.github.com/spec", first, fixed = TRUE)
}

dir.create("data", showWarnings = FALSE)

for (nm in names(files)) {
  f <- files[[nm]]
  dest <- file.path("data", f)

  if (file.exists(dest) && !is_lfs_pointer(dest)) {
    message("Already have ", f, ", skipping download")
  } else {
    ok <- FALSE
    for (u in urls_for(f)) {
      message("Downloading ", f, " from ", u)
      res <- tryCatch(download.file(u, dest, mode = "wb", quiet = TRUE),
                      error = function(e) 1L, warning = function(w) 1L)
      if (res == 0 && file.exists(dest) && !is_lfs_pointer(dest)) { ok <- TRUE; break }
    }
    if (!ok) stop("Could not get real data for ", f)
  }
}

# Final check: size of each file and whether it is a pointer
report <- data.frame(
  dataset     = names(files),
  file        = unname(files),
  size_MB     = round(file.size(file.path("data", files)) / 1024^2, 2),
  lfs_pointer = vapply(file.path("data", files), is_lfs_pointer, logical(1)),
  row.names   = NULL
)
print(report)

if (any(report$lfs_pointer)) stop("Some files are LFS pointers, not real data")
cat("\nAll", nrow(report), "files downloaded and none is a Git LFS pointer.\n")
