# 08_tcga_replication.R
# Session 7: do the main Beat AML findings hold in an independent cohort,
# TCGA-LAML (The Cancer Genome Atlas Research Network, NEJM 2013)?
#
#  1. Download TCGA-LAML from the NCI Genomic Data Commons (GDC) API: STAR gene
#     counts (GENCODE v36, GRCh38), the clinical patient table and the open
#     (masked) somatic mutation MAFs. Plus the TCGA paper's curated
#     cytogenetics from cBioPortal (study laml_tcga_pub), because GDC's own
#     fusion calls are controlled access. Everything is cached in data/tcga/.
#  2. Cohort with filters as close to Beat AML as TCGA allows. Filter table in
#     results/tcga_cohort_flow.csv, per-patient table in results/tcga_cohort.csv.
#  3. Match genes to Beat AML by Ensembl ID without the version suffix
#     (results/tcga_gene_matching.csv).
#  4. DESeq2 in TCGA, ~ sex + NPM1_status (approved), after a PCA check that
#     sequencing plate doesn't split the samples. Compare with Beat AML: fold
#     change correlation, overlap of significant genes, positive controls,
#     Verhaak NPM1 gene sets, and a FLT3-ITD sensitivity model.
#  5. Score every TCGA sample on the Beat AML top 50 up / top 50 down genes and
#     ask how well the score separates NPM1-mutated from wild-type (ROC AUC).
#  6. Three-group HOXA9 / MEIS1 plot: NPM1-mutated, KMT2A-rearranged, other AML.
#  7. Within NPM1-mutated: do NEK10, FTO and USP44 rise and fall together again,
#     and is there again a high and a low group? Same statistics as Session 4,
#     plus the Mer et al. primitive / committed markers from Session 5.
#
# Run from the repo root:  Rscript scripts/08_tcga_replication.R
# First run downloads about 650 MB (mostly the 151 count files) and takes about
# 10 minutes; later runs about 5 minutes.
# Needs data/ from 01_download.R (Beat AML gene IDs) and results/deseq2_npm1.csv
# (04_deseq2.R) and results/chrY_genes_grch37.csv (05_enrichment.R).

suppressPackageStartupMessages({
  library(jsonlite)
  library(data.table)
  library(dplyr)
  library(DESeq2)
  library(BiocParallel)
  library(ggplot2)
  library(ggrepel)
  library(fgsea)
  library(msigdbr)
  library(pROC)
  library(patchwork)
})

set.seed(2026)
options(timeout = 3600)
dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

tcga_dir <- "data/tcga"
gdc      <- "https://api.gdc.cancer.gov"
cbio     <- "https://www.cbioportal.org/api/studies/laml_tcga_pub"
dir.create(file.path(tcga_dir, "star"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(tcga_dir, "maf"),  recursive = TRUE, showWarnings = FALSE)
register(MulticoreParam(max(1, min(4, parallel::detectCores()))))

padj_cut   <- 0.05
lfc_cut    <- 1
coef_name  <- "NPM1_status_mut_vs_wt"
npm1_cols  <- c(wt = "grey55", mut = "#D55E00")
group_cols <- c("NPM1-mutated" = "#D55E00", "KMT2A-rearranged" = "#0072B2",
                "Other AML" = "grey55")

save_fig <- function(p, file, w = 7, h = 5) {
  ggsave(file.path("figures", file), p, width = w, height = h, dpi = 300)
  cat("Saved figures/", file, "\n", sep = "")
}

# ---- 1. Download from GDC -----------------------------------------------------
# The GDC "files" endpoint takes a JSON filter and returns file metadata. Each
# file is then fetched from /data/<file_id>. GDC gives an md5 checksum for every
# file, so we check each download against it (the same idea as the Git LFS
# pointer check in 01_download.R: a download can "succeed" and still be wrong).

gdc_files <- function(filters, fields) {
  body <- toJSON(list(filters = filters, fields = paste(fields, collapse = ","),
                      size = 2000, format = "JSON"), auto_unbox = TRUE)
  r <- httr::POST(paste0(gdc, "/files"), body = body, httr::content_type_json())
  httr::stop_for_status(r)
  fromJSON(httr::content(r, as = "text", encoding = "UTF-8"))$data$hits
}
f_and <- function(...) list(op = "and", content = list(...))
f_eq  <- function(field, value) list(op = "=", content = list(field = field, value = value))

fetch_gdc <- function(file_id, md5, dest) {
  if (file.exists(dest) && unname(tools::md5sum(dest)) == md5) return(invisible(TRUE))
  # Six tries, waiting 2, 4, ..., 64 s (about 2 minutes in all). Session 12's clean
  # run hit a few minutes of empty replies from GDC for one file; the old four
  # tries (30 s) gave up too early, and the same file downloaded fine later.
  for (attempt in 1:6) {
    ok <- tryCatch({ download.file(paste0(gdc, "/data/", file_id), dest,
                                   mode = "wb", quiet = TRUE); TRUE },
                   error = function(e) FALSE, warning = function(w) FALSE)
    if (ok && unname(tools::md5sum(dest)) == md5) return(invisible(TRUE))
    Sys.sleep(2^attempt)
  }
  stop("Download failed or md5 mismatch for ", file_id)
}

laml <- f_eq("cases.project.project_id", "TCGA-LAML")

# 1a. STAR gene counts. One file per RNA-seq aliquot.
star_manifest_file <- file.path(tcga_dir, "star_manifest.csv")
if (!file.exists(star_manifest_file)) {
  hits <- gdc_files(f_and(laml, f_eq("analysis.workflow_type", "STAR - Counts"),
                          f_eq("data_type", "Gene Expression Quantification")),
                    c("file_id", "file_name", "md5sum", "file_size", "cases.submitter_id",
                      "cases.samples.submitter_id", "cases.samples.sample_type",
                      "associated_entities.entity_submitter_id"))
  star_manifest <- data.frame(
    file_id     = hits$file_id, md5sum = hits$md5sum, file_size = hits$file_size,
    patient     = sapply(hits$cases, function(x) x$submitter_id[1]),
    sample      = sapply(hits$cases, function(x) x$samples[[1]]$submitter_id[1]),
    sample_type = sapply(hits$cases, function(x) x$samples[[1]]$sample_type[1]),
    aliquot     = sapply(hits$associated_entities, function(x) x$entity_submitter_id[1]))
  write.csv(star_manifest, star_manifest_file, row.names = FALSE)
}
star_manifest <- read.csv(star_manifest_file)
cat("STAR count files on GDC:", nrow(star_manifest), "\n")
for (i in seq_len(nrow(star_manifest))) {
  fetch_gdc(star_manifest$file_id[i], star_manifest$md5sum[i],
            file.path(tcga_dir, "star", paste0(star_manifest$file_id[i], ".tsv")))
}
cat("All", nrow(star_manifest), "STAR count files present and md5-checked\n")

# 1b. Clinical patient table (BCR Biotab). Rows 2 and 3 of the file are extra
# header lines (CDE names), not patients.
clin_file <- file.path(tcga_dir, "clinical_patient_laml.txt")
if (!file.exists(clin_file)) {
  hits <- gdc_files(f_and(laml, f_eq("data_format", "BCR Biotab"),
                          f_eq("file_name", "nationwidechildrens.org_clinical_patient_laml.txt")),
                    c("file_id", "md5sum"))
  fetch_gdc(hits$file_id[1], hits$md5sum[1], clin_file)
}

# 1c. Masked somatic mutations: one MAF per tumour/normal pair, open access.
# "Masked" means germline variants and low-quality calls are removed.
maf_manifest_file <- file.path(tcga_dir, "maf_manifest.csv")
if (!file.exists(maf_manifest_file)) {
  hits <- gdc_files(f_and(laml, f_eq("data_type", "Masked Somatic Mutation"),
                          f_eq("access", "open")),
                    c("file_id", "file_name", "md5sum", "cases.submitter_id"))
  write.csv(data.frame(file_id = hits$file_id, md5sum = hits$md5sum,
                       file_name = hits$file_name,
                       patient = sapply(hits$cases, function(x) x$submitter_id[1])),
            maf_manifest_file, row.names = FALSE)
}
maf_manifest <- read.csv(maf_manifest_file)
for (i in seq_len(nrow(maf_manifest))) {
  fetch_gdc(maf_manifest$file_id[i], maf_manifest$md5sum[i],
            file.path(tcga_dir, "maf", maf_manifest$file_name[i]))
}
cat("All", nrow(maf_manifest), "MAF files present and md5-checked\n")

# 1d. TCGA paper's curated cytogenetics and rearrangements, from cBioPortal.
# GDC's RNA fusion calls (Arriba, STAR-Fusion) are controlled access, and the
# GDC clinical table only partly records KMT2A (MLL) results. The NEJM 2013
# paper combined karyotype, FISH and whole-genome / RNA-seq rearrangement
# calls; cBioPortal serves that table openly.
cbio_file <- file.path(tcga_dir, "cbioportal_laml_tcga_pub_clinical.csv")
if (!file.exists(cbio_file)) {
  get_cbio <- function(type) {
    fromJSON(paste0(cbio, "/clinical-data?clinicalDataType=", type,
                    "&projection=SUMMARY&pageSize=100000"))
  }
  long <- rbind(get_cbio("PATIENT")[, c("patientId", "clinicalAttributeId", "value")],
                get_cbio("SAMPLE")[, c("patientId", "clinicalAttributeId", "value")])
  wide <- dcast(as.data.table(long), patientId ~ clinicalAttributeId,
                value.var = "value", fun.aggregate = function(x) x[1], fill = NA)
  write.csv(wide, cbio_file, row.names = FALSE)
}

# 1e. The paper's curated mutation calls for a few genes, also from cBioPortal.
# GDC's masked MAFs catch only 10 NPM1 insertions in 151 cases: like the Beat
# AML exome, the automated callers miss most of the 4-base NPM1 insertion. The
# paper's calls were reviewed by hand and cover all 200 patients.
cbio_mut_file <- file.path(tcga_dir, "cbioportal_laml_tcga_pub_mutations.csv")
if (!file.exists(cbio_mut_file)) {
  genes <- c(NPM1 = 4869, FLT3 = 2322, DNMT3A = 1788, TET2 = 54790, IDH1 = 3417,
             IDH2 = 3418, NRAS = 4893, KRAS = 3845, PTPN11 = 5781, KMT2A = 4297)
  m <- do.call(rbind, lapply(names(genes), function(g) {
    d <- fromJSON(paste0("https://www.cbioportal.org/api/molecular-profiles/",
                         "laml_tcga_pub_mutations/mutations?sampleListId=",
                         "laml_tcga_pub_all&projection=SUMMARY&entrezGeneId=", genes[[g]]))
    if (length(d) == 0) return(NULL)
    data.frame(gene = g, patient = d$patientId, protein_change = d$proteinChange,
               mutation_type = d$mutationType)
  }))
  write.csv(m, cbio_mut_file, row.names = FALSE)
}
cat("Downloads done.\n")

# ---- 2. Cohort ----------------------------------------------------------------

# 2a. Count matrix (unstranded column; TCGA-LAML RNA-seq is unstranded, so the
# three count columns give nearly the same totals and "unstranded" is the
# right one). The first four rows are STAR summary rows (N_unmapped, ...).
counts_rds <- file.path(tcga_dir, "tcga_laml_star_counts.rds")
if (!file.exists(counts_rds)) {
  one <- function(id) fread(file.path(tcga_dir, "star", paste0(id, ".tsv")),
                            skip = 1, select = c("gene_id", "gene_name", "gene_type",
                                                 "unstranded"))[!grepl("^N_", gene_id)]
  first <- one(star_manifest$file_id[1])
  mat <- sapply(star_manifest$file_id, function(id) {
    x <- one(id)
    if (!identical(x$gene_id, first$gene_id)) stop("Gene order differs in ", id)
    x$unstranded
  })
  rownames(mat) <- first$gene_id
  colnames(mat) <- star_manifest$aliquot
  saveRDS(list(counts = mat, genes = as.data.frame(first[, 1:3])), counts_rds)
}
tc <- readRDS(counts_rds)
cat("TCGA count matrix:", nrow(tc$counts), "genes x", ncol(tc$counts), "samples\n")

# 2b. Clinical table
clin <- read.delim(clin_file, check.names = FALSE)[-(1:2), ]
clin <- clin[, !duplicated(names(clin))]   # basophils_count appears twice in the file
na_str <- c("[Not Available]", "[Not Evaluated]", "[Unknown]", "[Not Applicable]")
clin[] <- lapply(clin, function(x) ifelse(x %in% na_str, NA, x))
cbio_tab <- read.csv(cbio_file)

# NPM1, FLT3 and PML-RARA from the clinical lab panel. The column holds strings
# like "FLT3 Mutation Positive|...|NPMc Positive". "NPMc" = cytoplasmic NPM1,
# the clinical NPM1 test (the counterpart of Beat AML's clinical NPM1 call).
lab_call <- function(txt, test) {
  m <- regmatches(txt, regexpr(paste0(test, " (Positive|Negative)"), txt))
  out <- rep(NA_character_, length(txt))
  out[grepl(paste0(test, " (Positive|Negative)"), txt)] <- sub(".* ", "", m)
  out
}
clin$NPM1_clinical <- lab_call(clin$molecular_abnormality_results, "NPMc")
clin$FLT3_clinical <- lab_call(clin$molecular_abnormality_results, "FLT3 Mutation")
clin$PMLRARA_lab   <- lab_call(clin$molecular_abnormality_results, "PML-RAR")

# MLL (KMT2A) FISH: test names and results are two "|"-separated lists that
# line up position by position. The result is the % of cells with a split signal.
tests <- strsplit(ifelse(is.na(clin$test_performed_indicator), "", clin$test_performed_indicator), "\\|")
outs  <- strsplit(ifelse(is.na(clin$fish_performed_outcome), "", clin$fish_performed_outcome), "\\|")
clin$MLL_fish_pct <- mapply(function(t, o) {
  if (length(t) == length(o) && "MLL" %in% t) suppressWarnings(as.numeric(o[t == "MLL"])) else NA
}, tests, outs)

# 2c. Somatic mutations from the MAFs. Keep protein-changing calls only.
maf_rds <- file.path(tcga_dir, "tcga_laml_maf.rds")
if (!file.exists(maf_rds)) {
  cols <- c("Hugo_Symbol", "Variant_Classification", "Variant_Type", "HGVSp_Short",
            "Tumor_Sample_Barcode", "t_depth", "t_alt_count")
  maf <- rbindlist(lapply(file.path(tcga_dir, "maf", maf_manifest$file_name),
                          function(f) fread(text = readLines(gzfile(f)),
                                            skip = "Hugo_Symbol", select = cols)))
  saveRDS(maf, maf_rds)
}
maf <- readRDS(maf_rds)
maf$patient <- substr(maf$Tumor_Sample_Barcode, 1, 12)
coding <- c("Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Ins", "Frame_Shift_Del",
            "In_Frame_Ins", "In_Frame_Del", "Splice_Site", "Nonstop_Mutation",
            "Translation_Start_Site")
maf <- maf[Variant_Classification %in% coding]
maf_patients <- unique(substr(maf_manifest$patient, 1, 12))
has_mut <- function(gene, filter = TRUE) {
  hit <- unique(maf[Hugo_Symbol == gene & eval(filter)]$patient)
  ifelse(clin$bcr_patient_barcode %in% maf_patients,
         clin$bcr_patient_barcode %in% hit, NA)
}
clin$NPM1_maf <- has_mut("NPM1")
cat("\nNPM1 in MAF, variant types:\n")
print(maf[Hugo_Symbol == "NPM1", .N, by = .(Variant_Classification, HGVSp_Short)])

# 2d. One row per TCGA-LAML patient with everything the filters need
cbio_mut <- read.csv(cbio_mut_file)
paper_has <- function(gene, pattern = ".") {
  unique(cbio_mut$patient[cbio_mut$gene == gene & grepl(pattern, cbio_mut$protein_change)])
}
# FLT3-ITD: in-frame insertions/duplications in the juxtamembrane domain
# (codons ~570 to 615). FLT3 point mutations (D835 etc.) are TKD, not ITD.
flt3 <- cbio_mut[cbio_mut$gene == "FLT3", ]
flt3_codon <- as.numeric(sub("^[A-Z]([0-9]+).*", "\\1", flt3$protein_change))
flt3_itd <- unique(flt3$patient[grepl("ins|dup", flt3$protein_change) &
                                flt3_codon >= 570 & flt3_codon <= 615])

pt <- clin %>%
  transmute(patient = bcr_patient_barcode,
            sex = factor(tools::toTitleCase(tolower(gender)), levels = c("Female", "Male")),
            age = as.numeric(age_at_diagnosis), fab = fab_category,
            cyto_risk = cyto_risk_group,
            prior_malignancy = history_other_malignancy,
            prior_heme_disorder = history_hematologic_disorder,
            NPM1_clinical, NPM1_maf, FLT3_clinical, PMLRARA_lab, MLL_fish_pct,
            karyotype_gdc = cytogenetic_abnormality_type,
            karyotype_gdc_other = cyto_abnormality_type_other) %>%
  left_join(cbio_tab %>% transmute(patient = patientId,
                                   cyto_code = CYTOGENETIC_CODE_OTHER,
                                   cyto_subtype = HISTOLOGICAL_SUBTYPE,
                                   rearrangements = INFERRED_GENOMIC_REARRANGEMENT,
                                   karyotype = CYTOGENETICS),
            by = "patient") %>%
  mutate(NPM1_paper = patient %in% paper_has("NPM1", "fs"),   # frameshift = NPM1c
         FLT3_ITD   = ifelse(patient %in% flt3_itd, "positive", "negative"))

# NPM1 status (approved in Session 7). Main call: the paper's curated
# sequencing call (all 200 patients, reviewed by hand). Where the paper and the
# clinical NPMc test disagree, the patient stays only if GDC's masked MAF
# confirms the paper call; otherwise the call is "unconfirmed" and the patient
# leaves the main cohort. Patients with no clinical result stay.
pt$NPM1 <- ifelse(pt$NPM1_paper, "mut", "wt")
pt$NPM1_confirmed <- is.na(pt$NPM1_clinical) |
  pt$NPM1_paper == (pt$NPM1_clinical == "Positive") |
  (!is.na(pt$NPM1_maf) & pt$NPM1_maf == pt$NPM1_paper)

# APL: FAB M3, or PML-RARA / t(15;17) in the curated or GDC cytogenetics
pt$APL <- pt$fab %in% "M3" | pt$cyto_code %in% "PML-RARA" |
  pt$cyto_subtype %in% "PML-RARA" | grepl("t\\(15;17\\)", pt$karyotype_gdc) |
  grepl("t\\(15;17\\)", pt$rearrangements)

# KMT2A (MLL) rearrangement, any of:
#  - the paper's curated cytogenetic code says "MLL translocation"
#  - a translocation involving 11q23 in the paper's rearrangement calls
#    (11q23.3 is KMT2A; t(10;11)(p12;q14) is PICALM-MLLT10, not KMT2A)
#  - t(9;11) or t(4;11) in the GDC karyotype fields
#  - MLL break-apart FISH positive in >= 20% of cells (positives are 74 to 98%;
#    negatives 0 to 19%, which is background for a break-apart probe)
pt$KMT2A_r_paper <- grepl("MLL", pt$cyto_code) | grepl("MLL", pt$cyto_subtype) |
  grepl("t\\([0-9;]+\\)\\([^)]*q23", pt$rearrangements) &
  grepl("t\\((11;|[0-9]+;11)", pt$rearrangements)
pt$KMT2A_r_gdc <- grepl("t\\(9;11\\)|t\\(4;11\\)", paste(pt$karyotype_gdc, pt$karyotype_gdc_other)) |
  (!is.na(pt$MLL_fish_pct) & pt$MLL_fish_pct >= 20)
pt$KMT2A_r <- pt$KMT2A_r_paper | pt$KMT2A_r_gdc

# Conflicting APL evidence we can't verify (approved drop, Session 7): the lab
# panel says PML-RARA positive, but FAB M2, normal karyotype, no t(15;17) in
# the paper's calls, and an NPM1 mutation. Not called APL, but not kept either.
pt$apl_conflict <- pt$PMLRARA_lab %in% "Positive" & !pt$APL

# De novo: Beat AML calls a case "not de novo" after any earlier cancer or
# blood disorder. TCGA records both as history fields.
pt$de_novo <- !(pt$prior_malignancy %in% "Yes" | pt$prior_heme_disorder %in% "YES")

# RNA-seq: one aliquot per patient in TCGA-LAML (checked below)
if (any(duplicated(star_manifest$patient))) stop("A patient has more than one RNA-seq file")
pt <- pt %>% left_join(star_manifest %>% select(patient, aliquot, sample_type), by = "patient")
pt$has_rna <- !is.na(pt$aliquot)
pt$plate   <- substr(pt$aliquot, 22, 25)   # sequencing plate, the closest thing to a batch

# 2e. Filter flow, same steps as Beat AML (02_cohort.R)
flow <- list(); step <- function(label, keep, note = "") {
  d <- pt[keep, ]
  flow[[length(flow) + 1]] <<- data.frame(
    step = length(flow), filter = label, patients = nrow(d),
    NPM1_mut = sum(d$NPM1 == "mut"), NPM1_wt = sum(d$NPM1 == "wt"),
    pct_mut = round(100 * mean(d$NPM1 == "mut"), 1), note = note)
  keep
}
k <- step("All TCGA-LAML patients", rep(TRUE, nrow(pt)))
k <- step("RNA-seq counts on GDC", k & pt$has_rna)
k <- step("AML at sampling", k, "all TCGA-LAML cases are AML")
k <- step("Initial diagnosis", k & pt$sample_type == "Primary Blood Derived Cancer - Peripheral Blood",
          "all samples are diagnostic (TCGA sample type 03)")
k <- step("NPM1 known and confirmed", k & pt$NPM1_confirmed,
          "paper call; if it disagrees with the clinical test, GDC MAF must confirm it")
k_before_npm1 <- k
k <- step("One sample per patient", k, "TCGA-LAML has one RNA-seq aliquot per patient")
k <- step("De novo AML only", k & pt$de_novo, "no earlier cancer or blood disorder")
k_before_apl <- k
k <- step("Drop APL", k & !pt$APL)
k_before_kmt2a <- k
k <- step("Drop KMT2A-rearranged", k & !pt$KMT2A_r)
k <- step("QC: conflicting APL evidence", k & !pt$apl_conflict,
          "TCGA-AB-2869: lab PML-RARA positive, FAB M2, normal karyotype")
flow <- do.call(rbind, flow)

pt$main_cohort <- k
pt$group <- case_when(pt$main_cohort & pt$NPM1 == "mut" ~ "NPM1_mut",
                      pt$main_cohort ~ "NPM1_wt",
                      k_before_kmt2a & pt$KMT2A_r ~ "KMT2A_r",
                      k_before_kmt2a & pt$apl_conflict ~ "QC_excluded",
                      TRUE ~ "excluded")
pt$exclusion_reason <- case_when(
  pt$main_cohort ~ NA_character_,
  !pt$has_rna ~ "no RNA-seq",
  !pt$NPM1_confirmed ~ "NPM1 call: paper and clinical test disagree, not confirmed by GDC MAF",
  !pt$de_novo ~ "not de novo",
  pt$APL ~ "APL",
  pt$KMT2A_r ~ "KMT2A-rearranged (kept for three-group plot)",
  pt$apl_conflict ~ "conflicting APL evidence: lab PML-RARA positive but FAB M2, normal karyotype, NPM1-mutated")

cat("\nTCGA-LAML cohort filters:\n"); print(flow, row.names = FALSE)
write.csv(flow, "results/tcga_cohort_flow.csv", row.names = FALSE)

# 2f. How do the NPM1 sources agree (RNA-seq patients)?
r <- pt[pt$has_rna, ]
cat("\nNPM1, paper call vs clinical NPMc test (patients with RNA-seq):\n")
print(table(paper = ifelse(r$NPM1_paper, "mut", "wt"), clinical_NPMc = r$NPM1_clinical, useNA = "ifany"))
cat("NPM1, paper call vs GDC masked MAF (NA = no MAF for that patient):\n")
print(table(paper = ifelse(r$NPM1_paper, "mut", "wt"), GDC_MAF = r$NPM1_maf, useNA = "ifany"))
disc <- r$NPM1_paper != (r$NPM1_clinical %in% "Positive")
cat("Discordant paper vs clinical:\n")
print(r[disc, c("patient", "NPM1_paper", "NPM1_clinical", "NPM1_maf", "main_cohort")], row.names = FALSE)

cat("\nKMT2A-rearranged with RNA-seq (sources agree?):\n")
print(r[r$KMT2A_r, c("patient", "KMT2A_r_paper", "KMT2A_r_gdc", "MLL_fish_pct", "NPM1",
                     "de_novo", "rearrangements")], row.names = FALSE)
cat("\nAPL with RNA-seq:", sum(r$APL), "(FAB M3:", sum(r$fab %in% "M3"), ")\n")
cat("Lab panel says PML-RAR positive but not called APL:\n")
print(r[r$PMLRARA_lab %in% "Positive" & !r$APL, c("patient", "fab", "cyto_code", "karyotype", "NPM1")],
      row.names = FALSE)

# 2g. Sex check from expression, as in Session 2: XIST (female) and RPS4Y1 (male)
lcpm <- function(g) {
  id <- tc$genes$gene_id[tc$genes$gene_name == g]
  log2(tc$counts[id, ] / colSums(tc$counts) * 1e6 + 1)
}
sx <- data.frame(aliquot = colnames(tc$counts), XIST = lcpm("XIST"), RPS4Y1 = lcpm("RPS4Y1"))
pt <- pt %>% left_join(sx, by = "aliquot")
pt$sex_expr <- ifelse(is.na(pt$XIST), NA, ifelse(pt$RPS4Y1 > pt$XIST, "Male", "Female"))
cat("\nRecorded sex vs expression sex (RNA-seq patients):\n")
print(table(recorded = pt$sex[pt$has_rna], expression = pt$sex_expr[pt$has_rna]))
mm <- pt$has_rna & pt$sex != pt$sex_expr
if (any(mm)) print(pt[mm, c("patient", "sex", "XIST", "RPS4Y1", "group")], row.names = FALSE)

cat("\nMain cohort:", sum(pt$main_cohort), "patients,", sum(pt$group == "NPM1_mut"),
    "NPM1-mutated (", round(100 * mean(pt$NPM1[pt$main_cohort] == "mut"), 1), "% ),",
    sum(pt$group == "NPM1_wt"), "wild-type; KMT2A-r kept for the three-group plot:",
    sum(pt$group == "KMT2A_r"), "\n")
cat("Sex in the main cohort (% female): NPM1-mut",
    round(100 * mean(pt$sex[pt$group == "NPM1_mut"] == "Female")), "| wt",
    round(100 * mean(pt$sex[pt$group == "NPM1_wt"] == "Female")), "\n")
cat("FLT3-ITD in the main cohort: NPM1-mut",
    round(100 * mean(pt$FLT3_ITD[pt$group == "NPM1_mut"] == "positive")), "% | wt",
    round(100 * mean(pt$FLT3_ITD[pt$group == "NPM1_wt"] == "positive")), "%\n")
cat("Sequencing plate by group (main cohort):\n")
print(table(plate = pt$plate[pt$main_cohort], NPM1 = pt$NPM1[pt$main_cohort]))

write.csv(pt %>% select(patient, aliquot, main_cohort, group, exclusion_reason, NPM1,
                        NPM1_paper, NPM1_clinical, NPM1_maf, NPM1_confirmed, FLT3_ITD, FLT3_clinical, sex, sex_expr, XIST, RPS4Y1, age,
                        fab, cyto_risk, cyto_code, de_novo, prior_malignancy,
                        prior_heme_disorder, APL, PMLRARA_lab, apl_conflict, KMT2A_r, KMT2A_r_paper,
                        KMT2A_r_gdc, MLL_fish_pct, rearrangements, karyotype, plate),
          "results/tcga_cohort.csv", row.names = FALSE)
cat("Saved results/tcga_cohort.csv and results/tcga_cohort_flow.csv\n")

# ---- 3. Match genes to Beat AML by Ensembl ID ---------------------------------
# TCGA uses GENCODE v36 on GRCh38, Beat AML uses Ensembl 75 on GRCh37. The
# Ensembl gene ID (ENSG...) stays the same for a gene across genome builds;
# only the ".15"-style version suffix changes when the gene model is edited,
# so we drop the suffix and match on the stable part. GENCODE also lists the
# chrY copies of pseudoautosomal genes as separate "_PAR_Y" rows (all zero in
# STAR counts); those are dropped.

tc$genes$ensg <- sub("\\..*$", "", tc$genes$gene_id)
par_y <- grepl("_PAR_Y$", tc$genes$gene_id)
tc$counts <- tc$counts[!par_y, ]; tc$genes <- tc$genes[!par_y, ]
if (any(duplicated(tc$genes$ensg))) stop("Duplicate Ensembl IDs after removing versions")
rownames(tc$counts) <- tc$genes$ensg

beat_all <- fread("data/beataml_waves1to4_counts_dbgap.txt",
                  select = c("stable_id", "display_label", "biotype"))
beat_de  <- read.csv("results/deseq2_npm1.csv")
in_tcga  <- function(ids) ids %in% tc$genes$ensg
gm <- data.frame(
  gene_set = c("All Beat AML genes (Ensembl 75)", "  protein coding",
               "Beat AML genes tested in DESeq2", "  protein coding",
               "Beat AML significant genes (padj < 0.05, |log2FC| > 1)",
               "Beat AML top 50 up + top 50 down"),
  n = c(nrow(beat_all), sum(beat_all$biotype == "protein_coding"), nrow(beat_de),
        sum(beat_de$biotype == "protein_coding"), sum(beat_de$sig), 100))
top_up   <- head(beat_de$gene_id[beat_de$log2FC > 0], 50)
top_down <- head(beat_de$gene_id[beat_de$log2FC < 0], 50)
gm$matched <- c(sum(in_tcga(beat_all$stable_id)),
                sum(in_tcga(beat_all$stable_id[beat_all$biotype == "protein_coding"])),
                sum(in_tcga(beat_de$gene_id)),
                sum(in_tcga(beat_de$gene_id[beat_de$biotype == "protein_coding"])),
                sum(in_tcga(beat_de$gene_id[beat_de$sig])),
                sum(in_tcga(c(top_up, top_down))))
gm$pct <- round(100 * gm$matched / gm$n, 1)
cat("\nGenes matched to TCGA (GENCODE v36) by Ensembl ID without version:\n")
print(gm, row.names = FALSE)
write.csv(gm, "results/tcga_gene_matching.csv", row.names = FALSE)

# Do matched IDs carry the same symbol? (Symbols get renamed over the years;
# the ID is the stable key, so a different symbol is not an error.)
m <- merge(beat_de[, c("gene_id", "symbol")], tc$genes[, c("ensg", "gene_name")],
           by.x = "gene_id", by.y = "ensg")
cat("Tested genes matched:", nrow(m), "| same symbol:", sum(m$symbol == m$gene_name),
    "| symbol changed:", sum(m$symbol != m$gene_name), "\n")
key_genes <- c("HOXA9", "HOXA10", "MEIS1", "PBX3", "CD34", "NEK10", "FTO", "USP44", "WDR86",
               "CDH2", "GPR12", "MDFI", "ZNF521", "PLXNB1", "HHIP",
               "CD163", "C1QA", "CD14", "MARCO", "MSR1")
key_ids <- beat_all$stable_id[match(key_genes, beat_all$display_label)]
cat("Key genes (controls, NEK10 block, Mer et al. markers) found in TCGA:",
    sum(in_tcga(key_ids)), "of", length(key_genes), "\n")
if (!all(in_tcga(key_ids))) cat("  missing:", key_genes[!in_tcga(key_ids)], "\n")
unmatched_sig <- beat_de[beat_de$sig & !in_tcga(beat_de$gene_id), c("gene_id", "symbol", "biotype")]
cat("Significant Beat AML genes with no TCGA match, by biotype:\n")
print(table(unmatched_sig$biotype))

# ---- 4. DESeq2 in TCGA ----------------------------------------------------------
# Approved design: ~ sex + NPM1_status (wt = reference). TCGA-LAML has no
# sequencing wave and no specimen-type split (all bone marrow, one centre), so
# sex is the only covariate. Sequencing plate is not in the model; the PCA below
# checks that plate doesn't split the samples.

main <- pt %>% filter(main_cohort) %>%
  mutate(NPM1_status = factor(NPM1, levels = c("wt", "mut")),
         FLT3_ITD    = factor(FLT3_ITD, levels = c("negative", "positive")),
         plate       = factor(plate))
rownames(main) <- main$aliquot
counts <- tc$counts[, main$aliquot]
if (!identical(colnames(counts), rownames(main))) stop("Sample order mismatch")
if (any(counts < 0) || any(counts != round(counts))) stop("Counts must be whole numbers")
storage.mode(counts) <- "integer"

dds <- DESeqDataSetFromMatrix(counts, colData = main, design = ~ sex + NPM1_status)
rowData(dds)$symbol  <- tc$genes$gene_name
rowData(dds)$biotype <- tc$genes$gene_type
# Same low-count rule as Beat AML: >= 10 reads in at least as many samples as
# the smaller group has (31 NPM1-mutated here, 100 in Beat AML)
n_small <- min(table(main$NPM1_status))
keep <- rowSums(counts(dds) >= 10) >= n_small
dds  <- dds[keep, ]
cat("\nTCGA main cohort:", ncol(dds), "samples; genes kept (>= 10 reads in >=", n_small,
    "samples):", nrow(dds), "of", length(keep), "\n")
cat("Beat AML tested genes among them:", sum(rownames(dds) %in% beat_de$gene_id), "\n")

# 4a. PCA (blind VST, top 500 variable genes), coloured by plate and by NPM1.
# R squared = share of each PC's variance explained by the variable.
vsd_qc <- vst(dds, blind = TRUE)
pv  <- apply(assay(vsd_qc), 1, var)
pca <- prcomp(t(assay(vsd_qc)[order(pv, decreasing = TRUE)[1:500], ]))
pct <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
pcs <- data.frame(main, pca$x[, 1:4])
r2 <- sapply(c("plate", "NPM1_status", "sex", "FLT3_ITD"), function(v)
  sapply(paste0("PC", 1:4), function(pc) summary(lm(pcs[[pc]] ~ pcs[[v]]))$r.squared))
cat("\nPCA variance explained (%):", paste0("PC", 1:4, " ", pct[1:4], collapse = ", "), "\n")
cat("R squared of each variable with PC1 to PC4:\n"); print(round(r2, 3))
write.csv(data.frame(PC = paste0("PC", 1:4), pct_var = pct[1:4], round(r2, 4)),
          "results/tcga_pca_r2.csv", row.names = FALSE)
plate_cols <- setNames(c("#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#0072B2",
                         "#D55E00", "#F0E442", "black")[seq_len(nlevels(main$plate))],
                       levels(main$plate))
p <- ggplot(pcs, aes(PC1, PC2, colour = plate, shape = NPM1_status)) +
  geom_point(size = 2.2, alpha = 0.85) +
  scale_colour_manual(values = plate_cols) +
  scale_shape_manual(values = c(wt = 16, mut = 17)) +
  labs(x = paste0("PC1 (", pct[1], "%)"), y = paste0("PC2 (", pct[2], "%)"),
       colour = "Plate", shape = "NPM1",
       title = paste0("TCGA-LAML main cohort (n = ", ncol(dds), "), blind VST, top 500 genes"),
       subtitle = sprintf("R squared with PC1 / PC2: plate %.2f / %.2f, NPM1 %.2f / %.2f",
                          r2["PC1", "plate"], r2["PC2", "plate"],
                          r2["PC1", "NPM1_status"], r2["PC2", "NPM1_status"])) +
  theme_bw()
save_fig(p, "08_tcga_pca_plate.png", w = 7.5, h = 5.5)

# 4b. Fit ~ sex + NPM1_status and shrink the NPM1 fold changes with apeglm,
# exactly as in 04_deseq2.R.
cat("\nFitting ~ sex + NPM1_status ...\n")
dds <- DESeq(dds, parallel = TRUE, quiet = TRUE)
res <- results(dds, name = coef_name, alpha = padj_cut)
shr <- lfcShrink(dds, coef = coef_name, type = "apeglm", res = res,
                 parallel = TRUE, quiet = TRUE)
tt <- data.frame(gene_id = rownames(res), symbol = rowData(dds)$symbol,
                 biotype = rowData(dds)$biotype, baseMean = res$baseMean,
                 log2FC = shr$log2FoldChange, lfcSE = shr$lfcSE,
                 log2FC_unshrunk = res$log2FoldChange, stat = res$stat,
                 pvalue = res$pvalue, padj = res$padj) %>%
  arrange(desc(abs(stat)))
tt$rank <- seq_len(nrow(tt))
tt$sig  <- !is.na(tt$padj) & tt$padj < padj_cut & abs(tt$log2FC) > lfc_cut
tt$direction <- ifelse(!tt$sig, "ns", ifelse(tt$log2FC > 0, "up", "down"))
tt$tested_in_beataml <- tt$gene_id %in% beat_de$gene_id
write.csv(tt, "results/tcga_deseq2_npm1.csv", row.names = FALSE)
cat("Saved results/tcga_deseq2_npm1.csv (", nrow(tt), " genes)\n", sep = "")
cat("TCGA: padj < 0.05:", sum(tt$padj < padj_cut, na.rm = TRUE),
    "| padj < 0.05 and |log2FC| > 1:", sum(tt$sig), "(up", sum(tt$direction == "up"),
    ", down", sum(tt$direction == "down"), ")\n")
cat("Top 10 in TCGA:", paste0(head(tt$symbol, 10), "(", head(tt$direction, 10), ")",
                               collapse = ", "), "\n")

# 4c. Compare with Beat AML on genes tested in both cohorts
cmp <- inner_join(beat_de %>% select(gene_id, symbol, log2FC, stat, padj, sig, direction, rank),
                  tt %>% select(gene_id, symbol_tcga = symbol, log2FC, stat, padj, sig,
                                direction, rank),
                  by = "gene_id", suffix = c("_beat", "_tcga"))
cmp$sig_tcga_05 <- !is.na(cmp$padj_tcga) & cmp$padj_tcga < padj_cut
cmp$same_sign   <- sign(cmp$log2FC_beat) == sign(cmp$log2FC_tcga)
beat_sig <- cmp$sig_beat; tcga_sig <- cmp$sig_tcga
ov <- table(beat = beat_sig, tcga = tcga_sig)
fisher_ov <- fisher.test(ov)
top50_beat <- cmp$rank_beat <= 50
cmp_summary <- data.frame(
  metric = c("Genes tested in both cohorts",
             "Pearson r, shrunken log2FC, all shared genes",
             "Spearman rho, shrunken log2FC, all shared genes",
             "Pearson r, shrunken log2FC, Beat AML significant genes",
             "Spearman rho, Wald statistic, all shared genes",
             "Beat AML significant genes (shared)",
             "TCGA significant genes (shared)",
             "Significant in both (same direction)",
             "Significant in both (opposite direction)",
             "Overlap odds ratio (Fisher)", "Overlap Fisher p",
             "Beat AML significant: same direction in TCGA (%)",
             "Beat AML significant: TCGA padj < 0.05, same direction",
             "Beat AML top 50 (shared): TCGA padj < 0.05, same direction"),
  value = c(nrow(cmp),
            cor(cmp$log2FC_beat, cmp$log2FC_tcga),
            cor(cmp$log2FC_beat, cmp$log2FC_tcga, method = "spearman"),
            cor(cmp$log2FC_beat[beat_sig], cmp$log2FC_tcga[beat_sig]),
            cor(cmp$stat_beat, cmp$stat_tcga, method = "spearman", use = "complete.obs"),
            sum(beat_sig), sum(tcga_sig),
            sum(beat_sig & tcga_sig & cmp$same_sign),
            sum(beat_sig & tcga_sig & !cmp$same_sign),
            unname(fisher_ov$estimate), fisher_ov$p.value,
            100 * mean(cmp$same_sign[beat_sig]),
            sum(beat_sig & cmp$sig_tcga_05 & cmp$same_sign),
            sum(top50_beat & cmp$sig_tcga_05 & cmp$same_sign)))
cat("\nBeat AML vs TCGA:\n"); print(cmp_summary, row.names = FALSE, digits = 3)
cat("Beat AML top 50 genes present in the TCGA test:", sum(top50_beat), "\n")
write.csv(cmp_summary, "results/tcga_vs_beataml_summary.csv", row.names = FALSE)

cmp$cat <- factor(case_when(beat_sig & tcga_sig ~ "Significant in both",
                            beat_sig ~ "Beat AML only", tcga_sig ~ "TCGA only",
                            TRUE ~ "Neither"),
                  levels = c("Significant in both", "Beat AML only", "TCGA only", "Neither"))
lab_genes <- c("HOXA9", "HOXA10", "MEIS1", "PBX3", "CD34", "NEK10", "FTO", "USP44",
               "WDR86", "TRH", "BAALC", "APP")
p <- ggplot(cmp %>% arrange(desc(cat)), aes(log2FC_beat, log2FC_tcga, colour = cat)) +
  geom_hline(yintercept = 0, colour = "grey60") + geom_vline(xintercept = 0, colour = "grey60") +
  geom_abline(linetype = "dashed", colour = "grey40") +
  geom_point(size = 0.7, alpha = 0.6) +
  geom_label_repel(data = filter(cmp, symbol %in% lab_genes), aes(label = symbol),
                   size = 2.8, colour = "black", fill = alpha("white", 0.8),
                   min.segment.length = 0, max.overlaps = Inf) +
  scale_colour_manual(values = c("Significant in both" = "#D55E00", "Beat AML only" = "#0072B2",
                                 "TCGA only" = "#009E73", "Neither" = "grey80")) +
  labs(x = "log2FC NPM1-mut vs wt, Beat AML (n = 276)",
       y = "log2FC NPM1-mut vs wt, TCGA-LAML (n = 111)", colour = "padj < 0.05, |log2FC| > 1",
       title = "NPM1 effect in Beat AML vs TCGA-LAML (apeglm-shrunken)",
       subtitle = sprintf("%d genes tested in both; Pearson r = %.2f, Spearman rho = %.2f",
                          nrow(cmp), cmp_summary$value[2], cmp_summary$value[3])) +
  coord_equal() + theme_bw()
save_fig(p, "08_tcga_lfc_vs_beataml.png", w = 8, h = 7)

# 4d. Positive controls in TCGA
hoxb <- sort(grep("^HOXB[0-9]+$", tt$symbol, value = TRUE))
pc_genes <- c("HOXA9", "HOXA10", hoxb, "MEIS1", "PBX3", "CD34")
pc <- tt %>% filter(symbol %in% pc_genes) %>%
  mutate(expected = ifelse(symbol == "CD34", "down", "up"),
         as_expected = ifelse(expected == "up", log2FC > 0, log2FC < 0) &
                       !is.na(padj) & padj < padj_cut) %>%
  left_join(beat_de %>% select(gene_id, log2FC_beataml = log2FC, padj_beataml = padj),
            by = "gene_id") %>%
  select(symbol, gene_id, baseMean, log2FC, padj, rank, expected, as_expected,
         log2FC_beataml, padj_beataml) %>%
  arrange(match(symbol, pc_genes))
cat("\nPositive controls in TCGA (HOXB genes that passed the TCGA low-count filter):\n")
print(pc %>% mutate(across(c(log2FC, log2FC_beataml), ~ round(.x, 2)),
                    across(c(padj, padj_beataml), ~ signif(.x, 2)), baseMean = round(baseMean)),
      row.names = FALSE)
write.csv(pc, "results/tcga_positive_controls.csv", row.names = FALSE)
core_controls <- c("HOXA9", "HOXA10", "MEIS1", "PBX3", "CD34")
if (!all(pc$as_expected[pc$symbol %in% core_controls])) {
  warning("A core positive control does not go the expected way in TCGA; look for a bug")
}

# 4e. Verhaak NPM1 gene sets on the TCGA ranking (Wald statistic), after
# removing male-specific chrY genes and XIST, as approved after Session 3.
# The chrY list is by Ensembl ID, which is the same in GRCh37 and GRCh38.
chrY <- read.csv("results/chrY_genes_grch37.csv")
sexrm <- tt$gene_id %in% chrY$id[chrY$region == "MSY"] | tt$symbol == "XIST"
cat("\nRemoved before GSEA:", sum(sexrm), "genes (male-specific chrY + XIST)\n")
ranks <- sort(setNames(tt$stat, tt$gene_id)[!sexrm & !is.na(tt$stat)], decreasing = TRUE)
vh_sets <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CGP") %>%
  filter(gs_name %in% c("VERHAAK_AML_WITH_NPM1_MUTATED_UP", "VERHAAK_AML_WITH_NPM1_MUTATED_DN"))
vh_sets <- split(vh_sets$ensembl_gene, vh_sets$gs_name) %>% lapply(unique)
set.seed(2026)
vh <- fgsea(vh_sets, ranks, minSize = 15, maxSize = 500, eps = 0)
vh$leadingEdge <- vapply(vh$leadingEdge, function(g) paste(tt$symbol[match(g, tt$gene_id)],
                                                           collapse = ";"), "")
cat("Verhaak NPM1 gene sets in TCGA:\n")
print(as.data.frame(vh)[, c("pathway", "size", "NES", "pval", "padj")], row.names = FALSE)
write.csv(as.data.frame(vh), "results/tcga_gsea_verhaak.csv", row.names = FALSE)

# 4f. FLT3-ITD sensitivity: ~ sex + FLT3_ITD + NPM1_status
cat("\nFLT3-ITD by NPM1 status (TCGA main cohort):\n")
print(table(FLT3_ITD = main$FLT3_ITD, NPM1 = main$NPM1_status))
dds_f <- dds
design(dds_f) <- ~ sex + FLT3_ITD + NPM1_status
dds_f <- DESeq(dds_f, parallel = TRUE, quiet = TRUE)
res_f <- results(dds_f, name = coef_name, alpha = padj_cut)
shr_f <- lfcShrink(dds_f, coef = coef_name, type = "apeglm", res = res_f,
                   parallel = TRUE, quiet = TRUE)
tt_f <- data.frame(gene_id = rownames(res_f), log2FC = shr_f$log2FoldChange,
                   stat = res_f$stat, padj = res_f$padj) %>%
  arrange(desc(abs(stat))) %>% mutate(rank = row_number())
tt_f$sig <- !is.na(tt_f$padj) & tt_f$padj < padj_cut & abs(tt_f$log2FC) > lfc_cut
top50_t <- head(tt$gene_id, 50)
both_f <- inner_join(tt[, c("gene_id", "log2FC")], tt_f[, c("gene_id", "log2FC")],
                     by = "gene_id", suffix = c("", "_flt3"))
sens <- data.frame(
  top50_kept_in_top50 = length(intersect(top50_t, head(tt_f$gene_id, 50))),
  top50_still_sig = sum(top50_t %in% tt_f$gene_id[tt_f$sig]),
  sig_main = sum(tt$sig), sig_flt3 = sum(tt_f$sig),
  sig_shared = length(intersect(tt$gene_id[tt$sig], tt_f$gene_id[tt_f$sig])),
  lfc_correlation = cor(both_f$log2FC, both_f$log2FC_flt3))
cat("FLT3-ITD sensitivity (TCGA):\n"); print(sens, row.names = FALSE, digits = 3)
cat("TCGA top 50 genes that leave the top 50:",
    paste(tt$symbol[match(setdiff(top50_t, head(tt_f$gene_id, 50)), tt$gene_id)],
          collapse = ", "), "\n")
write.csv(sens, "results/tcga_deseq2_sensitivity_flt3_summary.csv", row.names = FALSE)

# ---- 5. Beat AML signature score in TCGA ---------------------------------------
# The Beat AML top 50 up and top 50 down genes (ranked by |Wald statistic|, as
# everywhere else in this project), kept if the Ensembl ID exists in TCGA.
# One blind VST on the 111 main-cohort samples plus the 9 KMT2A-r cases is used
# for steps 5 to 7. Same low-count rule as the DESeq2 step, but the score genes,
# the NEK10 block and the Mer et al. markers are kept even if they fail it.

vs <- pt %>% filter(main_cohort | group == "KMT2A_r") %>%
  mutate(group3 = factor(case_when(group == "NPM1_mut" ~ "NPM1-mutated",
                                   group == "KMT2A_r"  ~ "KMT2A-rearranged",
                                   TRUE                ~ "Other AML"),
                         levels = names(group_cols)))
rownames(vs) <- vs$aliquot
cnt <- tc$counts[, vs$aliquot]; storage.mode(cnt) <- "integer"
beat_up   <- head(beat_de$gene_id[beat_de$stat > 0], 50)
beat_down <- head(beat_de$gene_id[beat_de$stat < 0], 50)
primitive <- c("CDH2", "GPR12", "MDFI", "ZNF521", "PLXNB1", "HHIP")
committed <- c("CD163", "C1QA", "CD14", "MARCO", "MSR1")
block_genes <- c("NEK10", "FTO", "USP44", "WDR86")
force_ids <- intersect(c(beat_up, beat_down, key_ids), rownames(cnt))
keep_v <- rowSums(cnt >= 10) >= n_small
keep_v[force_ids] <- TRUE
vst_mat <- assay(vst(DESeqDataSetFromMatrix(cnt[keep_v, ], colData = vs, design = ~ 1),
                     blind = TRUE))
cat("\nVST for steps 5 to 7:", ncol(vst_mat), "samples,", nrow(vst_mat), "genes\n")
# Look a gene up by its Beat AML symbol -> Beat AML Ensembl ID -> TCGA row
gene_vst <- function(s) {
  id <- beat_all$stable_id[beat_all$display_label == s]
  id <- intersect(id, rownames(vst_mat))
  if (length(id) != 1) stop(s, " matches ", length(id), " genes")
  vst_mat[id, ]
}

up_ok <- intersect(beat_up, rownames(vst_mat)); dn_ok <- intersect(beat_down, rownames(vst_mat))
cat("Score genes found in TCGA: up", length(up_ok), "of 50, down", length(dn_ok), "of 50\n")
cat("  missing:", beat_de$symbol[match(setdiff(c(beat_up, beat_down), c(up_ok, dn_ok)),
                                        beat_de$gene_id)], "\n")
# z-score each gene across the 111 main-cohort samples (no labels used), then
# score = mean z of the up genes minus mean z of the down genes.
is_main <- vs$main_cohort
zscore <- function(ids) {
  m <- vst_mat[ids, , drop = FALSE]
  (m - rowMeans(m[, is_main])) / apply(m[, is_main], 1, sd)
}
vs$score_up   <- colMeans(zscore(up_ok))
vs$score_down <- colMeans(zscore(dn_ok))
vs$score      <- vs$score_up - vs$score_down

mm <- vs[is_main, ]
roc_one <- function(x) roc(mm$NPM1, x, levels = c("wt", "mut"), direction = "<", quiet = TRUE)
roc_main <- roc_one(mm$score)
auc_tab <- do.call(rbind, lapply(c("score", "score_up", "score_down"), function(s) {
  r <- if (s == "score_down") roc(mm$NPM1, -mm[[s]], levels = c("wt", "mut"),
                                  direction = "<", quiet = TRUE) else roc_one(mm[[s]])
  ci <- ci.auc(r, method = "delong")
  data.frame(score = s, genes = c(score = length(up_ok) + length(dn_ok),
                                  score_up = length(up_ok), score_down = length(dn_ok))[[s]],
             AUC = as.numeric(auc(r)), CI_low = ci[1], CI_high = ci[3],
             wilcox_p = wilcox.test(mm[[s]] ~ mm$NPM1)$p.value)
}))
auc_tab$note <- c("up minus down", "up genes only", "down genes only (sign flipped)")
cat("\nROC AUC, NPM1-mutated vs wild-type (TCGA main cohort, n = 111):\n")
print(auc_tab, row.names = FALSE, digits = 3)
yj <- coords(roc_main, "best", best.method = "youden", ret = c("threshold", "sensitivity", "specificity"))
cat("Best cut-off (Youden, optimistic because chosen on the same data):",
    round(yj$threshold[1], 2), "| sensitivity", round(yj$sensitivity[1], 2),
    "| specificity", round(yj$specificity[1], 2), "\n")
write.csv(auc_tab, "results/tcga_signature_auc.csv", row.names = FALSE)
odd <- mm %>% filter((NPM1 == "wt" & score > yj$threshold[1]) |
                     (NPM1 == "mut" & score < yj$threshold[1])) %>%
  select(patient, NPM1, NPM1_clinical, score, FLT3_ITD, cyto_code, fab) %>% arrange(desc(score))
cat("Samples on the 'wrong' side of that cut-off:\n"); print(odd, row.names = FALSE, digits = 2)
write.csv(vs %>% select(patient, aliquot, group3, NPM1, score, score_up, score_down),
          "results/tcga_signature_scores.csv", row.names = FALSE)

n_lab3 <- vs %>% dplyr::count(group3) %>% mutate(lab = paste0(group3, "\n(n = ", n, ")"))
pa <- ggplot(vs, aes(group3, score, fill = group3)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, width = 0.6) +
  geom_jitter(width = 0.15, size = 1.1, alpha = 0.7) +
  scale_fill_manual(values = group_cols, guide = "none") +
  scale_x_discrete(labels = setNames(n_lab3$lab, n_lab3$group3)) +
  labs(x = NULL, y = "Score (mean z up - mean z down)",
       title = "Beat AML NPM1 signature in TCGA-LAML") +
  theme_bw(base_size = 11)
rd <- data.frame(fpr = 1 - roc_main$specificities, tpr = roc_main$sensitivities) %>% arrange(fpr, tpr)
pb <- ggplot(rd, aes(fpr, tpr)) +
  geom_abline(linetype = "dashed", colour = "grey50") +
  geom_step(colour = "#D55E00", linewidth = 1) +
  annotate("text", x = 0.55, y = 0.15, hjust = 0,
           label = sprintf("AUC = %.3f\n95%% CI %.3f to %.3f", auc_tab$AUC[1],
                           auc_tab$CI_low[1], auc_tab$CI_high[1])) +
  coord_equal() +
  labs(x = "1 - specificity", y = "Sensitivity",
       title = sprintf("ROC, NPM1-mut vs wt (n = %d)", nrow(mm))) +
  theme_bw(base_size = 11)
ggsave("figures/08_tcga_signature_score_roc.png", pa + pb + plot_layout(widths = c(1.2, 1)) +
         plot_annotation(subtitle = sprintf(
           "Beat AML top %d up and top %d down genes (matched by Ensembl ID), blind VST",
           length(up_ok), length(dn_ok))),
       width = 11, height = 5, dpi = 300)
cat("Saved figures/08_tcga_signature_score_roc.png\n")

# ---- 6. Three-group plot: HOXA9 and MEIS1 -------------------------------------
three <- do.call(rbind, lapply(c("HOXA9", "MEIS1"), function(g)
  data.frame(sample = colnames(vst_mat), gene = g, vst = gene_vst(g), group3 = vs$group3)))
three_summary <- three %>% group_by(gene, group3) %>%
  summarise(n = n(), median_vst = median(vst), q25 = quantile(vst, 0.25),
            q75 = quantile(vst, 0.75), .groups = "drop")
pairs <- list(c("NPM1-mutated", "Other AML"), c("KMT2A-rearranged", "Other AML"),
              c("NPM1-mutated", "KMT2A-rearranged"))
three_tests <- do.call(rbind, lapply(c("HOXA9", "MEIS1"), function(g) {
  do.call(rbind, lapply(pairs, function(pr) {
    a <- three$vst[three$gene == g & three$group3 == pr[1]]
    b <- three$vst[three$gene == g & three$group3 == pr[2]]
    data.frame(gene = g, comparison = paste(pr, collapse = " vs "),
               median_diff = median(a) - median(b), wilcox_p = wilcox.test(a, b)$p.value)
  }))
}))
write.csv(three_summary, "results/tcga_three_group_summary.csv", row.names = FALSE)
write.csv(three_tests, "results/tcga_three_group_tests.csv", row.names = FALSE)
cat("\nThree-group medians (VST), TCGA:\n"); print(as.data.frame(three_summary), digits = 3)
cat("Wilcoxon tests:\n"); print(three_tests, row.names = FALSE, digits = 3)
p <- ggplot(three, aes(group3, vst, fill = group3)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, width = 0.6) +
  geom_jitter(width = 0.15, size = 1.1, alpha = 0.7) +
  facet_wrap(~ gene, scales = "free_y") +
  scale_fill_manual(values = group_cols, guide = "none") +
  scale_x_discrete(labels = setNames(n_lab3$lab, n_lab3$group3)) +
  labs(title = "HOXA9 and MEIS1 in the two menin-dependent AML groups",
       subtitle = "De novo AML at diagnosis, TCGA-LAML", x = NULL,
       y = "Expression (VST, ~log2)") +
  theme_bw(base_size = 11)
save_fig(p, "08_tcga_three_group_hoxa9_meis1.png", w = 8, h = 5)

# ---- 7. NEK10 / FTO / USP44 (/ WDR86) within NPM1-mutated -----------------------
# Same statistics as Session 4 (05_enrichment.R step 8): Spearman correlations
# in mutated and wild-type samples, block score = mean z-score of the four
# genes, two-cluster k-means on it (seed 2026), the bimodality coefficient, the
# score against cohort variables, and each gene against FLT3-ITD.

mut <- vs %>% filter(group == "NPM1_mut")
wt  <- vs %>% filter(group == "NPM1_wt")
bm  <- t(sapply(block_genes, function(g) gene_vst(g)[mut$aliquot]))
cor_mut <- cor(t(bm), method = "spearman")
cor_wt  <- cor(sapply(block_genes, function(g) gene_vst(g)[wt$aliquot]), method = "spearman")
cat("\nSpearman correlation of the block genes within NPM1-mutated (n =", nrow(mut), "):\n")
print(round(cor_mut, 2))
cat("Same in NPM1-wild-type (n =", nrow(wt), "):\n"); print(round(cor_wt, 2))
# Is the NEK10-FTO-USP44 correlation stronger than chance? cor.test per pair
core <- c("NEK10", "FTO", "USP44")
pair_p <- combn(core, 2, function(pr) {
  ct <- cor.test(bm[pr[1], ], bm[pr[2], ], method = "spearman", exact = FALSE)
  data.frame(pair = paste(pr, collapse = "-"), rho_mut = unname(ct$estimate), p_mut = ct$p.value,
             rho_wt = cor(gene_vst(pr[1])[wt$aliquot], gene_vst(pr[2])[wt$aliquot],
                          method = "spearman"))
}, simplify = FALSE) %>% bind_rows()
cat("Core pairs:\n"); print(pair_p, row.names = FALSE, digits = 3)

bimodality <- function(x) {
  sk <- mean(scale(x)^3); ku <- mean(scale(x)^4) - 3; n <- length(x)
  (sk^2 + 1) / (ku + 3 * (n - 1)^2 / ((n - 2) * (n - 3)))
}
mut$block_score <- colMeans(t(scale(t(bm))))
mut$core_score  <- colMeans(t(scale(t(bm[core, ]))))
set.seed(2026)
km <- kmeans(mut$block_score, centers = 2, nstart = 25)
mut$block_group <- factor(ifelse(km$cluster == which.max(km$centers), "high", "low"),
                          levels = c("low", "high"))
set.seed(2026)
km_core <- kmeans(mut$core_score, centers = 2, nstart = 25)
mut$core_group <- factor(ifelse(km_core$cluster == which.max(km_core$centers), "high", "low"),
                         levels = c("low", "high"))
bc_block <- bimodality(mut$block_score); bc_core <- bimodality(mut$core_score)
cat("\nk-means split of the block score (4 genes): low", sum(mut$block_group == "low"),
    "| high", sum(mut$block_group == "high"),
    sprintf("(%.0f%% high)", 100 * mean(mut$block_group == "high")), "\n")
cat("k-means split of the core score (NEK10/FTO/USP44): low", sum(mut$core_group == "low"),
    "| high", sum(mut$core_group == "high"), "; agreement with block split:",
    sum(mut$core_group == mut$block_group), "of", nrow(mut), "\n")
cat("Bimodality coefficient: block score", round(bc_block, 3), "| core score",
    round(bc_core, 3), "(> 0.555 hints at two groups)\n")
# How well separated are the two k-means groups? Gap between groups in SD units
cat("Block score, gap between the highest 'low' and lowest 'high' sample:",
    round(min(mut$block_score[mut$block_group == "high"]) -
            max(mut$block_score[mut$block_group == "low"]), 2), "\n")

cat_vars <- c("FLT3_ITD", "sex", "plate", "cyto_risk")
assoc <- do.call(rbind, lapply(cat_vars, function(v) {
  x <- factor(mut[[v]]); ok <- !is.na(x)
  med <- tapply(mut$block_score[ok], x[ok], median)
  data.frame(variable = v, test = "Kruskal-Wallis (score); Fisher (high/low)",
             levels = paste(sprintf("%s: n=%d, median %.2f", names(med),
                                    as.integer(table(x[ok])), med), collapse = " | "),
             p_score = kruskal.test(mut$block_score[ok], x[ok])$p.value,
             p_split = fisher.test(table(x[ok], mut$block_group[ok]),
                                   simulate.p.value = nlevels(droplevels(x[ok])) > 2,
                                   B = 1e4)$p.value)
}))
ct <- cor.test(mut$block_score, mut$age, method = "spearman", exact = FALSE)
assoc <- rbind(assoc, data.frame(variable = "age", test = "Spearman correlation",
                                 levels = sprintf("rho = %.2f", ct$estimate), p_score = ct$p.value,
                                 p_split = wilcox.test(age ~ block_group, data = mut)$p.value))
cat("\nBlock score vs cohort variables (NPM1-mutated only):\n")
print(assoc[, c("variable", "p_score", "p_split")], row.names = FALSE, digits = 3)
write.csv(assoc, "results/tcga_block_associations.csv", row.names = FALSE)

flt3_pos <- mut$FLT3_ITD == "positive"
per_gene <- do.call(rbind, lapply(c(block_genes, "core_score", "block_score"), function(g) {
  v <- if (g %in% block_genes) bm[g, ] else mut[[g]]
  data.frame(feature = g, median_FLT3_ITD_neg = median(v[!flt3_pos]),
             median_FLT3_ITD_pos = median(v[flt3_pos]),
             wilcox_p = wilcox.test(v[flt3_pos], v[!flt3_pos])$p.value)
}))
cat("\nEach block gene vs FLT3-ITD within NPM1-mutated (FLT3-ITD+", sum(flt3_pos), "of",
    nrow(mut), "):\n")
print(per_gene, row.names = FALSE, digits = 3)
write.csv(per_gene, "results/tcga_block_flt3_by_gene.csv", row.names = FALSE)

# 7b. Mer et al. 2021 primitive and committed markers (as in 06_subgroup.R)
zmat <- function(genes) t(scale(t(do.call(rbind, lapply(genes, function(g) gene_vst(g)[mut$aliquot])))))
mut$primitive_score <- colMeans(zmat(primitive))
mut$committed_score <- colMeans(zmat(committed))
mut$prim_minus_comm <- mut$primitive_score - mut$committed_score
marker_tab <- do.call(rbind, lapply(c(primitive, committed), function(g) {
  v <- gene_vst(g)[mut$aliquot]
  data.frame(gene = g, set = ifelse(g %in% primitive, "primitive", "committed"),
             median_low = median(v[mut$block_group == "low"]),
             median_high = median(v[mut$block_group == "high"]),
             diff_high_minus_low = median(v[mut$block_group == "high"]) -
               median(v[mut$block_group == "low"]),
             wilcox_p = wilcox.test(v ~ mut$block_group)$p.value,
             spearman_vs_core = cor(v, mut$core_score, method = "spearman"))
}))
marker_tab$wilcox_padj <- p.adjust(marker_tab$wilcox_p, method = "BH")
cat("\nPrimitive / committed markers, high vs low block group (TCGA):\n")
print(marker_tab, row.names = FALSE, digits = 3)
write.csv(marker_tab, "results/tcga_primitive_committed_markers.csv", row.names = FALSE)
score_tab <- do.call(rbind, lapply(c("primitive_score", "committed_score", "prim_minus_comm"),
                                   function(s) {
  ct <- cor.test(mut[[s]], mut$core_score, method = "spearman", exact = FALSE)
  data.frame(score = s, median_low = median(mut[[s]][mut$block_group == "low"]),
             median_high = median(mut[[s]][mut$block_group == "high"]),
             wilcox_p = wilcox.test(mut[[s]] ~ mut$block_group)$p.value,
             spearman_vs_core = unname(ct$estimate), spearman_p = ct$p.value)
}))
cat("Scores vs block group:\n"); print(score_tab, row.names = FALSE, digits = 3)
write.csv(score_tab, "results/tcga_primitive_committed_scores.csv", row.names = FALSE)
set.seed(2026)
km_pc <- kmeans(mut$prim_minus_comm, centers = 2, nstart = 25)
mut$pc_call <- factor(ifelse(km_pc$cluster == which.max(km_pc$centers),
                             "primitive-like", "committed-like"),
                      levels = c("committed-like", "primitive-like"))
pc_tab <- table(block = mut$block_group, call = mut$pc_call)
cat("Block group vs primitive/committed call:\n"); print(pc_tab)
cat("Fisher p =", signif(fisher.test(pc_tab)$p.value, 3), "\n")

# 7c. Extra: co-mutations, high vs low (the paper's calls; all patients have
# mutation data). Session 5 saw TET2/IDH2 in the high group and RAS-pathway
# mutations and DNMT3A in the low group. Small numbers, descriptive only.
comut_genes <- c("DNMT3A", "TET2", "IDH1", "IDH2", "NRAS", "KRAS", "PTPN11")
comut <- do.call(rbind, lapply(c("FLT3_ITD", comut_genes), function(g) {
  has <- if (g == "FLT3_ITD") flt3_pos else
    mut$patient %in% cbio_mut$patient[cbio_mut$gene == g & cbio_mut$mutation_type != "Silent"]
  tb <- table(factor(has, levels = c(FALSE, TRUE)), mut$block_group)
  data.frame(feature = g, n_low = sum(mut$block_group == "low"), mut_low = tb["TRUE", "low"],
             n_high = sum(mut$block_group == "high"), mut_high = tb["TRUE", "high"],
             fisher_p = fisher.test(tb)$p.value)
}))
tet_idh <- mut$patient %in% cbio_mut$patient[cbio_mut$gene %in% c("TET2", "IDH2")]
comut <- rbind(comut, data.frame(feature = "TET2 or IDH2", n_low = sum(mut$block_group == "low"),
  mut_low = sum(tet_idh & mut$block_group == "low"), n_high = sum(mut$block_group == "high"),
  mut_high = sum(tet_idh & mut$block_group == "high"),
  fisher_p = fisher.test(table(tet_idh, mut$block_group))$p.value))
# Post hoc (added after seeing the IDH1 result): IDH1 and IDH2 mutants both make
# 2-hydroxyglutarate, which blocks TET2, so "TET2, IDH1 or IDH2" is one pathway.
tet_idh12 <- mut$patient %in% cbio_mut$patient[cbio_mut$gene %in% c("TET2", "IDH1", "IDH2")]
for (grp in c("block_group", "core_group")) {
  tb <- table(factor(tet_idh12, levels = c(FALSE, TRUE)), mut[[grp]])
  comut <- rbind(comut, data.frame(
    feature = paste0("TET2, IDH1 or IDH2 (post hoc", if (grp == "core_group") ", core split" else "", ")"),
    n_low = sum(mut[[grp]] == "low"), mut_low = tb["TRUE", "low"],
    n_high = sum(mut[[grp]] == "high"), mut_high = tb["TRUE", "high"],
    fisher_p = fisher.test(tb)$p.value))
}
cat("\nCo-mutations by block group (TCGA, NPM1-mutated):\n"); print(comut, row.names = FALSE, digits = 3)
write.csv(comut, "results/tcga_block_comutations.csv", row.names = FALSE)

# 7d. Check with the three-gene core split. In TCGA, WDR86 doesn't follow
# NEK10/FTO/USP44 (as in Beat AML, it tracks FLT3-ITD instead), so the
# four-gene score can pull in samples that are high for WDR86 only. Repeat the
# key 7b/7c tests with the NEK10/FTO/USP44 k-means split.
core_checks <- rbind(
  data.frame(test = "primitive_score, Wilcoxon high vs low",
             value = median(mut$primitive_score[mut$core_group == "high"]) -
               median(mut$primitive_score[mut$core_group == "low"]),
             p = wilcox.test(mut$primitive_score ~ mut$core_group)$p.value),
  data.frame(test = "committed_score, Wilcoxon high vs low",
             value = median(mut$committed_score[mut$core_group == "high"]) -
               median(mut$committed_score[mut$core_group == "low"]),
             p = wilcox.test(mut$committed_score ~ mut$core_group)$p.value),
  data.frame(test = "primitive-like calls in core-high (value = count of high)",
             value = sum(mut$pc_call == "primitive-like" & mut$core_group == "high"),
             p = fisher.test(table(mut$core_group, mut$pc_call))$p.value),
  do.call(rbind, lapply(c("FLT3_ITD", comut_genes), function(g) {
    has <- if (g == "FLT3_ITD") flt3_pos else
      mut$patient %in% cbio_mut$patient[cbio_mut$gene == g & cbio_mut$mutation_type != "Silent"]
    data.frame(test = paste0(g, ": mutated in core-high (value) of ", sum(mut$core_group == "high"),
                             "; in core-low ", sum(has & mut$core_group == "low"), " of ",
                             sum(mut$core_group == "low")),
               value = sum(has & mut$core_group == "high"),
               p = fisher.test(table(factor(has, levels = c(FALSE, TRUE)), mut$core_group))$p.value)
  })))
cat("\nChecks with the NEK10/FTO/USP44 (core) split, high", sum(mut$core_group == "high"),
    "vs low", sum(mut$core_group == "low"), ":\n")
print(core_checks, row.names = FALSE, digits = 3)
write.csv(core_checks, "results/tcga_block_core_split_checks.csv", row.names = FALSE)

write.csv(cbind(mut %>% select(patient, aliquot, FLT3_ITD, sex, age, plate, block_score,
                               core_score, block_group, core_group, primitive_score,
                               committed_score, prim_minus_comm, pc_call), t(bm)),
          "results/tcga_block_scores.csv", row.names = FALSE)
block_summary <- data.frame(
  cohort = "TCGA-LAML", n_NPM1_mut = nrow(mut),
  rho_NEK10_FTO = cor_mut["NEK10", "FTO"], rho_NEK10_USP44 = cor_mut["NEK10", "USP44"],
  rho_FTO_USP44 = cor_mut["FTO", "USP44"],
  rho_WDR86_min = min(cor_mut["WDR86", core]), rho_WDR86_max = max(cor_mut["WDR86", core]),
  n_high = sum(mut$block_group == "high"), n_low = sum(mut$block_group == "low"),
  bimodality_block = bc_block, bimodality_core = bc_core,
  rho_primitive_vs_core = score_tab$spearman_vs_core[1],
  fisher_p_block_vs_pc_call = fisher.test(pc_tab)$p.value)
write.csv(block_summary, "results/tcga_block_summary.csv", row.names = FALSE)

long <- data.frame(sample = rep(mut$aliquot, each = 4),
                   gene = factor(rep(block_genes, nrow(mut)), levels = block_genes),
                   vst = as.vector(bm), FLT3_ITD = rep(mut$FLT3_ITD, each = 4),
                   group = rep(mut$block_group, each = 4))
long$sample <- factor(long$sample, levels = mut$aliquot[order(mut$block_score)])
p <- ggplot(long, aes(sample, vst, colour = FLT3_ITD, shape = group)) +
  geom_point(size = 1.6) +
  facet_wrap(~ gene, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = c(negative = "grey55", positive = "#CC79A7")) +
  scale_shape_manual(values = c(low = 16, high = 17), name = "k-means group") +
  labs(title = "NEK10 / FTO / USP44 / WDR86 within NPM1-mutated AML, TCGA-LAML",
       subtitle = "Samples ordered by block score (mean z-score of the four genes)",
       x = paste0("NPM1-mutated samples (n = ", nrow(mut), ")"), y = "VST") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
save_fig(p, "08_tcga_block_nek10_fto_usp44_wdr86.png", w = 7, h = 7)

pl <- rbind(
  data.frame(group = mut$block_group, value = mut$primitive_score,
             what = "Primitive score\n(CDH2, GPR12, MDFI, ZNF521, PLXNB1, HHIP)"),
  data.frame(group = mut$block_group, value = mut$committed_score,
             what = "Committed score\n(CD163, C1QA, CD14, MARCO, MSR1)"))
p_lab <- setNames(sprintf("Wilcoxon p = %.2g", score_tab$wilcox_p[1:2]), unique(pl$what))
p <- ggplot(pl, aes(group, value, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.6) +
  geom_jitter(width = 0.15, size = 1.2, alpha = 0.7) +
  facet_wrap(~ what, labeller = labeller(what = function(x) paste0(x, "\n", p_lab[x]))) +
  scale_fill_manual(values = c(low = "grey65", high = "#D55E00"), guide = "none") +
  labs(title = "Mer et al. 2021 markers in the NEK10/FTO/USP44 groups, TCGA-LAML",
       subtitle = paste0("NPM1-mutated AML (n = ", nrow(mut), "); score = mean z-score of the markers"),
       x = "NEK10 / FTO / USP44 / WDR86 group (k-means)", y = "Score") +
  theme_bw(base_size = 10)
save_fig(p, "08_tcga_primitive_committed.png", w = 7, h = 4.5)

writeLines(c(paste("GDC data release:", fromJSON(paste0(gdc, "/status"))$data_release),
             paste("MSigDB version:", unique(msigdbr(species = "Homo sapiens",
                                                       collection = "H")$db_version)),
             capture.output(sessionInfo())),
           "results/sessionInfo_08_tcga_replication.txt")
cat("\nDone.\n")
