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

# ---- Part 2: apply the approved filters --------------------------------------
# Plan approved 2026-10-26. Each step keeps a subset of the previous one and we
# record how many samples and patients are left, so the cohort is traceable.

flow <- list()
log_step <- function(step, description, df) {
  flow[[length(flow) + 1]] <<- data.frame(
    step = step, description = description,
    samples = nrow(df), patients = n_distinct(df$dbgap_subject_id),
    npm1_positive = sum(df$NPM1 == "positive", na.rm = TRUE),
    npm1_negative = sum(df$NPM1 == "negative", na.rm = TRUE)
  )
  df
}

s0 <- log_step(0, "All samples in clinical summary", clin)

# 1. Has an RNA-seq sample that made it into the counts matrix (passed QC)
s1 <- s0 %>%
  filter(!is.na(dbgap_rnaseq_sample), dbgap_rnaseq_sample %in% counts_ids) %>%
  log_step(1, "Has RNA-seq and is in counts matrix", .)

# 2. Diagnosis at the time of sampling is AML (drops MDS, MPN, mixed lineage)
s2 <- s1 %>%
  filter(dxAtSpecimenAcquisition ==
           "ACUTE MYELOID LEUKAEMIA (AML) AND RELATED PRECURSOR NEOPLASMS") %>%
  log_step(2, "AML at specimen acquisition", .)

# 3. Sample taken at initial diagnosis (not relapse, remission, residual)
s3 <- s2 %>%
  filter(diseaseStageAtSpecimenCollection == "Initial Diagnosis") %>%
  log_step(3, "Initial diagnosis sample", .)

# 4. NPM1 status known (clinical lab call)
s4 <- s3 %>%
  filter(NPM1 %in% c("positive", "negative")) %>%
  log_step(4, "NPM1 status known", .)

# 5. One sample per patient: earliest collection time, bone marrow first on ties
s5 <- s4 %>%
  mutate(is_bm = specimenType == "Bone Marrow Aspirate") %>%
  arrange(dbgap_subject_id, timeOfSampleCollectionRelativeToInclusion, desc(is_bm)) %>%
  distinct(dbgap_subject_id, .keep_all = TRUE) %>%
  select(-is_bm) %>%
  log_step(5, "One sample per patient (earliest, bone marrow first)", .)

stopifnot(!any(duplicated(s5$dbgap_subject_id)))

# Flags used by steps 6 to 8. consensusAMLFusions comes from manual review of
# clinical cytogenetics plus RNA-seq fusion calls. "KMT2A_re" (partner not
# named) and "MLLT3-KMT2A" (t(9;11)) are KMT2A fusions, not KMT2A-PTD.
s5 <- s5 %>%
  mutate(KMT2A_rearranged = grepl("KMT2A", consensusAMLFusions),
         APL = grepl("PML-RARA", consensusAMLFusions))

# 6. De novo AML only. isDenovo is FALSE for AML that transformed from MDS/MPN,
#    therapy-related AML, and anyone with a previous cancer. 1 "unknown" dropped.
s6 <- s5 %>%
  filter(isDenovo == "TRUE") %>%
  log_step(6, "De novo AML only", .)

# 7. Drop APL (PML-RARA): very different biology and treatment
s7 <- s6 %>%
  filter(!APL) %>%
  log_step(7, "Drop APL (PML-RARA)", .)

# 8. Drop KMT2A fusions from the main comparison. They are menin-dependent
#    like NPM1-mutated AML, so they stay in cohort.csv as a separate group.
s8 <- s7 %>%
  filter(!KMT2A_rearranged) %>%
  log_step(8, "Drop KMT2A fusions (kept in cohort.csv, main_cohort = FALSE)", .)

flow <- bind_rows(flow) %>%
  mutate(npm1_mut_pct = round(100 * npm1_positive / (npm1_positive + npm1_negative), 1))
cat("\n-- Cohort flow --\n")
print(flow, row.names = FALSE)

# cohort.csv holds the main cohort plus the KMT2A-fusion cases from step 7.
# main_cohort marks the 2-group NPM1 comparison; group is for plotting.
cohort <- s7 %>%
  transmute(
    dbgap_subject_id, dbgap_rnaseq_sample, dbgap_dnaseq_sample,
    NPM1 = factor(ifelse(NPM1 == "positive", "mut", "wt"), levels = c("wt", "mut")),
    main_cohort = !KMT2A_rearranged,
    group = ifelse(KMT2A_rearranged, "KMT2A_r",
                   ifelse(NPM1 == "mut", "NPM1_mut", "NPM1_wt")),
    FLT3_ITD = `FLT3-ITD`,
    specimenType, cohort, centerID, consensus_sex, ageAtDiagnosis,
    isDenovo, specificDxAtAcquisition, consensusAMLFusions,
    KMT2A_rearranged, ELN2017
  )

cat("\nGroups:\n");                          print(table(cohort$group))
cat("\nKMT2A fusions by NPM1 status:\n");    print(table(cohort$consensusAMLFusions[cohort$KMT2A_rearranged],
                                                         cohort$NPM1[cohort$KMT2A_rearranged]))
main <- filter(cohort, main_cohort)
cat("\nMain cohort, specimen type by NPM1:\n"); print(table(main$specimenType, main$NPM1))

dir.create("results", showWarnings = FALSE)
write.csv(flow, "results/cohort_flow.csv", row.names = FALSE)
write.csv(cohort, "results/cohort.csv", row.names = FALSE)
cat("\nSaved results/cohort_flow.csv and results/cohort.csv\n")
