# 08_tcga_replication.R
# Session 6: do the main Beat AML findings hold in an independent cohort,
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
#  Steps 4 to 7 (DESeq2 and comparison, signature score, three-group plot,
#  NEK10/FTO/USP44 block) follow once the cohort is approved.
#
# Run from the repo root:  Rscript scripts/08_tcga_replication.R
# First run downloads about 650 MB (mostly the 151 count files).
# Needs results/deseq2_npm1.csv (04_deseq2.R) for the later steps.

suppressPackageStartupMessages({
  library(jsonlite)
  library(data.table)
  library(dplyr)
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
  for (attempt in 1:4) {
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

# NPM1 status. Main call: the paper's curated sequencing call (all 200
# patients, reviewed by hand). The clinical NPMc test is kept for comparison.
pt$NPM1 <- ifelse(pt$NPM1_paper, "mut", "wt")

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
k <- step("NPM1 known", k & !is.na(pt$NPM1_paper), "paper call exists for all 200")
k <- step("One sample per patient", k, "TCGA-LAML has one RNA-seq aliquot per patient")
k <- step("De novo AML only", k & pt$de_novo, "no earlier cancer or blood disorder")
k_before_apl <- k
k <- step("Drop APL", k & !pt$APL)
k_before_kmt2a <- k
k <- step("Drop KMT2A-rearranged", k & !pt$KMT2A_r)
flow <- do.call(rbind, flow)

pt$main_cohort <- k
pt$group <- case_when(pt$main_cohort & pt$NPM1 == "mut" ~ "NPM1_mut",
                      pt$main_cohort ~ "NPM1_wt",
                      k_before_kmt2a & pt$KMT2A_r ~ "KMT2A_r",
                      TRUE ~ "excluded")

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

write.csv(pt %>% select(patient, aliquot, main_cohort, group, NPM1, NPM1_paper, NPM1_clinical,
                        NPM1_maf, FLT3_ITD, FLT3_clinical, sex, sex_expr, XIST, RPS4Y1, age,
                        fab, cyto_risk, cyto_code, de_novo, prior_malignancy,
                        prior_heme_disorder, APL, PMLRARA_lab, KMT2A_r, KMT2A_r_paper,
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
