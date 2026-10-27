# 11_fto_locus.R
# Session 9: is the "FTO-high" signal in NPM1-mutated AML FTO itself, or a
# hidden transcript inside the FTO gene (IFEX9, Arza-Apalategi et al.,
# Leukemia 2023)? Beat AML only gives one count per gene, so we look at read
# coverage along the FTO locus in TCGA-LAML, re-processed by recount3.
#
#  1. Find TCGA-LAML in recount3 and check that per-base coverage (bigWig) and
#     splice-junction files exist for our samples. Record the genome build,
#     annotation and read lengths (results/fto_recount3_availability.csv).
#
# Run from the repo root:  Rscript scripts/11_fto_locus.R
# Needs results/tcga_cohort.csv and results/tcga_block_scores.csv from
# 08_tcga_replication.R (Session 7). Downloads are cached in data/recount3/.

suppressPackageStartupMessages({
  library(recount3)
  library(dplyr)
})

set.seed(2026)
options(timeout = 3600)
dir.create("results", showWarnings = FALSE)
r3_dir <- "data/recount3"
dir.create(r3_dir, recursive = TRUE, showWarnings = FALSE)
# recount3 caches its own downloads with BiocFileCache; keep them in data/ too
bfc <- BiocFileCache::BiocFileCache(file.path(r3_dir, "cache"), ask = FALSE)

# ---- 1. What does recount3 have for TCGA-LAML? ---------------------------------
# recount3 (Wilks et al., Genome Biol 2021) re-aligned every TCGA RNA-seq run
# with one pipeline (STAR to GRCh38 / hg38). For each sample it stores a bigWig
# file: the number of reads covering every single base of the genome. For each
# project it stores a junction table: every place where a read was split
# across an intron (one end on one exon, the other on the next).

hp   <- available_projects(organism = "human", bfc = bfc)
proj <- hp %>% filter(project == "LAML", file_source == "tcga")
if (nrow(proj) != 1) stop("TCGA-LAML not found (or found twice) in recount3")
cat("recount3 project:", proj$project, "| source:", proj$file_source,
    "| samples:", proj$n_samples, "\n")

md <- read_metadata(file_retrieve(locate_url("LAML", proj$project_home,
                                             type = "metadata"), bfc = bfc))
md$patient <- substr(md$tcga.tcga_barcode, 1, 12)
cat("Metadata rows:", nrow(md), "| patients:", n_distinct(md$patient), "\n")

# Match to the Session 7 cohort by RNA-seq aliquot barcode (exact), not only by
# patient, so we are sure it is the same sequencing library.
cohort <- read.csv("results/tcga_cohort.csv")
blocks <- read.csv("results/tcga_block_scores.csv")
main   <- cohort %>% filter(main_cohort)
main$in_recount3 <- main$aliquot %in% md$tcga.tcga_barcode
cat("Main cohort (111) aliquots found in recount3:", sum(main$in_recount3), "of",
    nrow(main), "\n")
if (!all(main$in_recount3)) stop("Some main-cohort aliquots are missing in recount3")
cat("NPM1-mutated:", sum(main$NPM1 == "mut"), "| wild-type:", sum(main$NPM1 == "wt"),
    "| core split (Session 7):", paste(names(table(blocks$core_group)),
                                       table(blocks$core_group), collapse = ", "), "\n")

# Check that every bigWig we need is really there (HTTP HEAD: status and size,
# no download). A small file would mean an error page, not coverage.
main <- main %>% left_join(md %>% select(aliquot = tcga.tcga_barcode, external_id,
                                         rail_id), by = "aliquot")
main$bw_url <- vapply(main$external_id, function(s)
  unname(locate_url("LAML", proj$project_home, type = "bw", sample = s)), "")
head_one <- function(u) {
  h <- httr::HEAD(u)
  c(status = httr::status_code(h), bytes = as.numeric(httr::headers(h)[["content-length"]]))
}
hd <- t(vapply(main$bw_url, head_one, c(status = 0, bytes = 0)))
main$bw_status <- hd[, "status"]; main$bw_mb <- round(hd[, "bytes"] / 1e6)
cat("bigWig files: HTTP 200 for", sum(main$bw_status == 200), "of", nrow(main),
    "| size", min(main$bw_mb), "to", max(main$bw_mb), "MB each,",
    round(sum(main$bw_mb) / 1e3, 1), "GB in total\n")

jx_urls <- locate_url("LAML", proj$project_home, type = "jxn")
jx <- t(vapply(jx_urls, head_one, c(status = 0, bytes = 0)))
cat("Junction files (project level, all samples in one table):\n")
print(data.frame(file = basename(jx_urls), status = jx[, "status"],
                 MB = round(jx[, "bytes"] / 1e6, 1)), row.names = FALSE)

# Annotation recount3 uses for TCGA: GENCODE v26 (also v29) on GRCh38 / hg38.
gtf_url <- paste0(getOption("recount3_url", "http://duffel.rail.bio/recount3"),
                  "/human/annotations/exon_sums/human.exon_sums.G026.gtf.gz")
gtf_h <- head_one(gtf_url)
cat("GENCODE v26 exon annotation (recount3):", gtf_h[["status"]],
    round(gtf_h[["bytes"]] / 1e6, 1), "MB\n")

# Sequencing depth and read length for the depth adjustment later. bc_auc is
# the total coverage of the sample (sum of reads x bases over the genome).
qc <- md %>% transmute(aliquot = tcga.tcga_barcode,
                       read_length = `recount_qc.star.average_input_read_length`,
                       mapped_reads = `recount_qc.star.all_mapped_reads`,
                       auc_all = `recount_qc.bc_auc.all_reads_all_bases`,
                       junction_count = recount_qc.junction_count)
main <- main %>% left_join(qc, by = "aliquot")
cat("Average input read length:", paste(names(table(main$read_length)),
                                        table(main$read_length), sep = ": n=", collapse = "; "), "\n")
cat("Mapped reads (millions): median", round(median(main$mapped_reads) / 1e6, 1),
    "range", round(min(main$mapped_reads) / 1e6, 1), "to",
    round(max(main$mapped_reads) / 1e6, 1), "\n")

avail <- main %>%
  left_join(blocks %>% select(patient, core_group), by = "patient") %>%
  mutate(fto_group = case_when(NPM1 == "wt" ~ "NPM1_wt",
                               core_group == "high" ~ "FTO_high",
                               core_group == "low" ~ "FTO_low")) %>%
  select(patient, aliquot, NPM1, fto_group, FLT3_ITD, sex, plate, external_id,
         rail_id, bw_status, bw_mb, read_length, mapped_reads, auc_all, junction_count)
write.csv(avail, "results/fto_recount3_availability.csv", row.names = FALSE)
cat("Groups:", paste(names(table(avail$fto_group)), table(avail$fto_group),
                     collapse = ", "), "\n")
cat("Saved results/fto_recount3_availability.csv\n")
