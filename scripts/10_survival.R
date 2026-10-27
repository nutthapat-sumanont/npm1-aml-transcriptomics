# 10_survival.R
# Session 9: an exploratory look at overall survival. Treatment varied a lot
# across Beat AML patients (intensive chemotherapy, azacitidine or decitabine,
# trial drugs, transplant or not), and we can't model that properly, so every
# result here is exploratory: a direction and a rough size, not a verdict.
#
#  1. Find and explain the survival columns in the Beat AML clinical summary:
#     what time zero is, how censoring works, how long follow-up is. Describe
#     treatment per group (report only, not in any model).
#  2. Main cohort: Kaplan-Meier and Cox, NPM1-mutated vs wild-type, adjusted
#     for age. Then the four NPM1 x FLT3-ITD groups. Does NPM1-mutated without
#     FLT3-ITD (ELN 2022 favorable) do better? Checked again after taking the
#     other favorable-risk AML (CBF, CEBPA) out of the comparison group.
#  3. Within NPM1-mutated: NEK10/FTO/USP44-high vs -low (three-gene split,
#     Session 8: 37 high, 63 low), Cox adjusted for age and FLT3-ITD. Hazard
#     ratio, 95% CI, deaths per group, and what that means for power.
#  4. Same as step 3 in TCGA-LAML (three-gene split from Session 7, 6 high /
#     25 low), with the TCGA paper's overall survival from cBioPortal.
#
# Run from the repo root:  Rscript scripts/10_survival.R   (under a minute)
# Needs data/ (01_download.R), results/cohort.csv (02),
# results/block_nek10_fto_scores.csv (05), results/subgroup_per_sample.csv
# (06), results/tcga_block_scores.csv (08). Downloads the TCGA paper's clinical
# table from cBioPortal into data/tcga/ if 08 hasn't already.

suppressPackageStartupMessages({
  library(survival)
  library(data.table)
  library(dplyr)
  library(readxl)
  library(jsonlite)
  library(ggplot2)
  library(patchwork)
})

set.seed(2026)   # nothing here is random; set for consistency with the rest
dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
dir.create("data/tcga", recursive = TRUE, showWarnings = FALSE)

clin_file  <- "data/beataml_wv1to4_clinical.xlsx"
cbio_file  <- "data/tcga/cbioportal_laml_tcga_pub_clinical.csv"   # same cache as 08
cbio       <- "https://www.cbioportal.org/api/studies/laml_tcga_pub"
days_per_month <- 365.25 / 12
npm1_cols  <- c(wt = "grey55", mut = "#D55E00")
grp4_cols  <- c("NPM1-wt / ITD-neg" = "grey55", "NPM1-wt / ITD-pos" = "black",
                "NPM1-mut / ITD-neg" = "#D55E00", "NPM1-mut / ITD-pos" = "#0072B2")
block_cols <- c(low = "#009E73", high = "#D55E00")

save_fig <- function(p, file, w = 7, h = 5) {
  ggsave(file.path("figures", file), p, width = w, height = h, dpi = 300)
  cat("Saved figures/", file, "\n", sep = "")
}

# One row per term of a Cox model: HR, 95% CI, Wald p, plus n and deaths.
cox_rows <- function(fit, model, cohort) {
  s  <- summary(fit)
  ci <- s$conf.int
  data.frame(cohort = cohort, model = model, term = rownames(ci),
             HR = ci[, "exp(coef)"], CI_low = ci[, "lower .95"],
             CI_high = ci[, "upper .95"], p = s$coefficients[, "Pr(>|z|)"],
             n = fit$n, deaths = fit$nevent, row.names = NULL)
}

# Proportional hazards check (Schoenfeld residuals): p < 0.05 means the hazard
# ratio seems to change over time, so a single HR is a poor summary.
ph_rows <- function(fit, model, cohort) {
  z <- cox.zph(fit)$table
  data.frame(cohort = cohort, model = model, term = rownames(z),
             chisq = z[, "chisq"], p = z[, "p"], row.names = NULL)
}

# Per-group Kaplan-Meier summary: n, deaths, median OS (months) with 95% CI,
# and survival at 1, 2 and 3 years.
km_rows <- function(fit, cohort, comparison) {
  tb <- summary(fit)$table
  if (is.null(dim(tb))) tb <- t(as.matrix(tb))
  s  <- summary(fit, times = c(12, 24, 36), extend = TRUE)
  grp <- if (is.null(s$strata)) factor(rep("all", length(s$time))) else s$strata
  at  <- function(t) tapply(s$surv[s$time == t], grp[s$time == t], identity)
  out <- data.frame(cohort = cohort, comparison = comparison,
                    group = sub("^[^=]*=", "", rownames(tb) %||% "all"),
                    n = tb[, "records"], deaths = tb[, "events"],
                    median_OS_months = tb[, "median"],
                    median_CI_low = tb[, "0.95LCL"], median_CI_high = tb[, "0.95UCL"],
                    row.names = NULL)
  for (t in c(12, 24, 36)) out[[paste0("OS_", t / 12, "yr")]] <- as.numeric(at(t))
  out
}
`%||%` <- function(a, b) if (is.null(a)) b else a

# Kaplan-Meier plot with numbers at risk underneath (ggplot only; survminer
# isn't installed and isn't needed).
km_plot <- function(fit, cols, title, legend_title, xmax = 60, by = 12) {
  s <- summary(fit, censored = TRUE)
  strata <- if (is.null(fit$strata)) rep("all", length(fit$time)) else
    rep(sub("^[^=]*=", "", names(fit$strata)), fit$strata)
  d <- data.frame(time = fit$time, surv = fit$surv, n.censor = fit$n.censor,
                  group = strata)
  # start every curve at time 0, survival 1
  d <- rbind(data.frame(time = 0, surv = 1, n.censor = 0, group = unique(strata)), d)
  d$group <- factor(d$group, levels = names(cols))
  p <- ggplot(d, aes(time, surv, colour = group)) +
    geom_step(linewidth = 0.7) +
    geom_point(data = d[d$n.censor > 0, ], shape = 3, size = 1.4) +
    scale_colour_manual(values = cols, name = legend_title, drop = FALSE) +
    scale_x_continuous(breaks = seq(0, xmax, by)) +
    scale_y_continuous(labels = scales::percent) +
    coord_cartesian(xlim = c(0, xmax), ylim = c(0, 1)) +
    guides(colour = guide_legend(nrow = 2)) +
    labs(x = "Months from time zero", y = "Overall survival", title = title) +
    theme_bw(base_size = 10) + theme(legend.position = "bottom")
  times <- seq(0, xmax, by)
  r <- summary(fit, times = times, extend = TRUE)
  rgrp <- if (is.null(r$strata)) "all" else sub("^[^=]*=", "", as.character(r$strata))
  risk <- data.frame(time = r$time, n.risk = r$n.risk, group = rgrp)
  risk$group <- factor(risk$group, levels = rev(names(cols)))
  tab <- ggplot(risk, aes(time, group, label = n.risk, colour = group)) +
    geom_text(size = 2.8, show.legend = FALSE) +
    scale_colour_manual(values = cols) +
    scale_x_continuous(breaks = times) +
    coord_cartesian(xlim = c(0, xmax)) +
    labs(x = NULL, y = NULL, subtitle = "Number at risk") +
    theme_minimal(base_size = 9) +
    theme(panel.grid = element_blank(), axis.text.x = element_blank())
  p / tab + plot_layout(heights = c(4, 1))
}

# Power of a two-group Cox / log-rank test (Schoenfeld 1983): with D deaths
# and a fraction p of patients in one group, the test has power
#   Phi( sqrt(D p (1 - p)) * |log HR| - z_(1 - alpha/2) )
# and the smallest HR it can pick up with 80% power is
#   exp( (z_(1 - alpha/2) + z_0.8) / sqrt(D p (1 - p)) ).
# Covariate adjustment changes this only a little, so it's a fair rough guide.
power_rows <- function(D, p, cohort, comparison, hrs = c(1.5, 2, 3), alpha = 0.05) {
  za <- qnorm(1 - alpha / 2)
  k  <- sqrt(D * p * (1 - p))
  pw <- pnorm(k * log(hrs) - za)
  needed <- function(hr) ceiling((za + qnorm(0.8))^2 / (p * (1 - p) * log(hr)^2))
  data.frame(cohort = cohort, comparison = comparison, deaths = D,
             fraction_in_group = round(p, 3),
             power_HR_1.5 = pw[1], power_HR_2 = pw[2], power_HR_3 = pw[3],
             min_HR_80pct_power = exp((za + qnorm(0.8)) / k),
             deaths_needed_HR_1.5 = needed(1.5), deaths_needed_HR_2 = needed(2))
}

cox_all <- list(); ph_all <- list(); km_all <- list(); pow_all <- list()

# ---- 1. Survival columns in the Beat AML clinical summary --------------------
# overallSurvival: days. vitalStatus: Alive / Dead / Unknown. causeOfDeath:
# Alive, Dead-Disease, Dead-Treatment, Dead-Other, Dead-Unknown.
#
# Time zero. The data dictionary isn't reachable from here and Bottomly et al.
# 2022 only say "post-diagnosis survival, measured in days". The data
# themselves pin it down: overallSurvival is a property of the PATIENT, not of
# the sample (it's identical across all samples of a patient), and it's never
# shorter than the gap between study inclusion and any of the patient's
# samples (apart from a handful of rounding cases). So time zero is the day
# the patient entered the study. For our main cohort that's the same day as
# the diagnostic sample (checked below), so in practice time zero = diagnosis.
#
# Censoring. vitalStatus "Dead": overallSurvival is days to death (event).
# "Alive": overallSurvival is days to last contact; we know they lived at
# least that long, but not what happened afterwards. That's a censored
# observation: it counts towards the "at risk" group up to that day and then
# drops out without counting as a death. "Unknown": no follow-up at all
# (overallSurvival 0), so the patient can't be used.

clin <- as.data.frame(read_excel(clin_file))
by_pt <- split(clin$overallSurvival, clin$dbgap_subject_id)
n_multi <- sum(lengths(by_pt) > 1)
n_same  <- sum(sapply(by_pt[lengths(by_pt) > 1], function(x) length(unique(x)) == 1))
n_short <- sum(clin$overallSurvival < clin$timeOfSampleCollectionRelativeToInclusion,
               na.rm = TRUE)
cat("Patients with more than one sample:", n_multi, "| same overallSurvival on all:",
    n_same, "\n")
cat("Samples where overallSurvival < days from inclusion to that sample:", n_short,
    "of", nrow(clin), "\n")
print(table(vitalStatus = clin$vitalStatus, causeOfDeath = clin$causeOfDeath,
            useNA = "ifany"))

cohort <- read.csv("results/cohort.csv") %>% filter(main_cohort)
stopifnot(nrow(cohort) == 276)
surv_cols <- c("dbgap_rnaseq_sample", "overallSurvival", "vitalStatus", "causeOfDeath",
               "timeOfSampleCollectionRelativeToInclusion", "ageAtSpecimenAcquisition",
               "typeInductionTx", "responseToInductionTx",
               "cumulativeTreatmentTypes", "cumulativeTreatmentRegimens")
ba <- cohort %>% left_join(clin[, surv_cols], by = "dbgap_rnaseq_sample")
cat("\nMain cohort, days from study inclusion to the RNA-seq sample:\n")
print(table(ba$timeOfSampleCollectionRelativeToInclusion))
cat("Main cohort vital status by NPM1:\n")
print(table(ba$vitalStatus, ba$NPM1))

# Treatment, for description only. "Intensive" = ever had an anthracycline /
# cytarabine-based intensive regimen (7+3 and variants, HiDAC/MiDAC, FLAG,
# CLAG, MEC, CPX-351, ...). "Transplant" = any bone marrow transplant recorded.
# Transplant happens months after time zero, so putting it in a Cox model as a
# plain yes/no would give the transplanted patients credit for the time they
# had to survive to get there (immortal time bias). We only describe it.
intensive_re <- paste0("7\\+3|3\\+7|5\\+2|HiDAC|HIDAC|MiDAC|MIDAC|IDAC|FLAG|CLAG|",
                       "GCLAC|MEC \\(|HAM \\(|ADE \\(|CPX-351|Vyxeos|COG-AAML|",
                       "SWOG S1203|CALGB|Daunorubicin|Idarubicin")
ba <- ba %>%
  mutate(
    os_days  = overallSurvival,
    os_months = overallSurvival / days_per_month,
    death    = as.integer(vitalStatus == "Dead"),
    NPM1_status = factor(NPM1, levels = c("wt", "mut")),
    FLT3_ITD = factor(FLT3_ITD, levels = c("negative", "positive")),
    age10    = ageAtDiagnosis / 10,
    intensive = grepl(intensive_re, cumulativeTreatmentRegimens),
    transplant = grepl("Bone Marrow Transplant", cumulativeTreatmentTypes),
    treatment_recorded = !is.na(cumulativeTreatmentRegimens),
    grp4 = factor(paste0("NPM1-", NPM1, " / ITD-",
                         ifelse(FLT3_ITD == "positive", "pos", "neg")),
                  levels = names(grp4_cols)),
    # ELN 2022 favorable via NPM1: NPM1-mutated without FLT3-ITD (and without
    # adverse cytogenetics; the 2 NPM1-mutated ELN2017-adverse cases are ITD+)
    npm1_fav = factor(ifelse(grp4 == "NPM1-mut / ITD-neg", "NPM1-mut / ITD-neg",
                             "everyone else"),
                      levels = c("everyone else", "NPM1-mut / ITD-neg")),
    # Other favorable-risk AML in the wild-type group: CBF (RUNX1-RUNX1T1,
    # CBFB-MYH11) and CEBPA, i.e. NPM1-wt with ELN2017 "Favorable".
    other_fav = NPM1 == "wt" & ELN2017 %in% "Favorable"
  )

ba_s <- ba %>% filter(vitalStatus %in% c("Alive", "Dead"))
cat("\nDropped for unknown vital status (no follow-up):",
    paste(setdiff(ba$dbgap_rnaseq_sample, ba_s$dbgap_rnaseq_sample), collapse = ", "),
    "\n")
cat("Survival set:", nrow(ba_s), "patients,", sum(ba_s$death), "deaths\n")

# Median follow-up by the reverse Kaplan-Meier method: flip the event
# indicator, so "censored alive" becomes the event. It answers "how long were
# patients followed", without the short times of the patients who died early.
fu <- survfit(Surv(os_months, 1 - death) ~ 1, data = ba_s)
fu_med <- summary(fu)$table["median"]
cat("Median follow-up (reverse KM):", round(fu_med, 1), "months\n")

surv_def <- data.frame(
  item = c("time column", "event column", "time zero",
           "event", "censored", "not usable",
           "patients with >1 sample, same OS on all", "samples with OS < inclusion-to-sample gap",
           "main cohort samples taken on the inclusion day",
           "main cohort, usable for survival", "deaths", "median follow-up (reverse KM), months"),
  value = c("overallSurvival (days)", "vitalStatus (Alive / Dead / Unknown)",
            "study inclusion; for the main cohort the same day as the diagnostic sample",
            "vitalStatus = Dead: overallSurvival is days to death",
            "vitalStatus = Alive: overallSurvival is days to last contact",
            "vitalStatus = Unknown (overallSurvival = 0)",
            paste(n_same, "of", n_multi), paste(n_short, "of", nrow(clin)),
            paste(sum(ba$timeOfSampleCollectionRelativeToInclusion == 0), "of", nrow(ba)),
            nrow(ba_s), sum(ba_s$death), round(fu_med, 1)))
write.csv(surv_def, "results/survival_definitions.csv", row.names = FALSE)

treat <- ba_s %>%
  mutate(group = as.character(grp4)) %>%
  bind_rows(ba_s %>% mutate(group = paste0("NPM1-", NPM1, " (all)"))) %>%
  group_by(group) %>%
  summarise(n = n(), median_age = median(ageAtDiagnosis, na.rm = TRUE),
            treatment_recorded = sum(treatment_recorded),
            intensive = sum(intensive), transplant = sum(transplant),
            CR_after_induction = sum(grepl("^Complete Response", responseToInductionTx)),
            refractory = sum(responseToInductionTx %in% "Refractory"),
            deaths = sum(death), .groups = "drop")
write.csv(treat, "results/survival_treatment_summary.csv", row.names = FALSE)
cat("\nTreatment and response by group (description only):\n")
print(as.data.frame(treat))

# ---- 2. NPM1-mutated vs wild-type, then NPM1 x FLT3-ITD ----------------------
# Kaplan-Meier: the share of patients still alive over time, recalculated at
# each death among those still being followed. Log-rank: are the curves
# different? Cox: hazard ratio (HR), the death rate in one group relative to
# the other at any given time; HR < 1 = fewer deaths. Age goes in as a
# straight line per 10 years, because age is the strongest known predictor of
# AML survival and NPM1-mutated patients are a bit older here.
# 1 patient has no age, so the age-adjusted models use 274.

km_npm1 <- survfit(Surv(os_months, death) ~ NPM1_status, data = ba_s)
lr_npm1 <- survdiff(Surv(os_months, death) ~ NPM1_status, data = ba_s)
cox_npm1_u <- coxph(Surv(os_months, death) ~ NPM1_status, data = ba_s)
cox_npm1 <- coxph(Surv(os_months, death) ~ NPM1_status + age10, data = ba_s)
cat("\nNPM1 log-rank p =", signif(1 - pchisq(lr_npm1$chisq, 1), 3), "\n")
print(summary(cox_npm1))

km_g4 <- survfit(Surv(os_months, death) ~ grp4, data = ba_s)
lr_g4 <- survdiff(Surv(os_months, death) ~ grp4, data = ba_s)
cox_g4 <- coxph(Surv(os_months, death) ~ grp4 + age10, data = ba_s)
cat("\nFour-group log-rank p =", signif(1 - pchisq(lr_g4$chisq, 3), 3), "\n")
print(summary(cox_g4))

# ELN 2022 favorable via NPM1 (NPM1-mut, no FLT3-ITD) vs everyone else.
cox_fav <- coxph(Surv(os_months, death) ~ npm1_fav + age10, data = ba_s)
# Same, after taking CBF and CEBPA AML (the other ELN favorable groups) out of
# the comparison group, so it's NPM1-mut/ITD-neg vs AML that isn't favorable
# by genetics. Report only; tells us how much the first HR is diluted.
ba_nofav <- ba_s %>% filter(!other_fav)
cox_fav2 <- coxph(Surv(os_months, death) ~ npm1_fav + age10, data = ba_nofav)
# And NPM1-mut/ITD-neg vs the other favorable AML directly (both "favorable").
ba_favs <- ba_s %>% filter(other_fav | npm1_fav == "NPM1-mut / ITD-neg") %>%
  mutate(fav_type = factor(ifelse(other_fav, "CBF / CEBPA", "NPM1-mut / ITD-neg"),
                           levels = c("CBF / CEBPA", "NPM1-mut / ITD-neg")))
cox_fav3 <- coxph(Surv(os_months, death) ~ fav_type + age10, data = ba_favs)
cat("\nOther favorable-risk AML in the wild-type group:",
    sum(ba_s$other_fav), "patients\n")
print(table(ba_s$consensusAMLFusions[ba_s$other_fav], useNA = "ifany"))

cox_all <- c(cox_all, list(
  cox_rows(cox_npm1_u, "NPM1 (unadjusted)", "Beat AML"),
  cox_rows(cox_npm1, "NPM1 + age", "Beat AML"),
  cox_rows(cox_g4, "NPM1 x FLT3-ITD (4 groups) + age; ref NPM1-wt / ITD-neg", "Beat AML"),
  cox_rows(cox_fav, "NPM1-mut / ITD-neg vs everyone else + age", "Beat AML"),
  cox_rows(cox_fav2, "NPM1-mut / ITD-neg vs everyone else, CBF/CEBPA removed + age", "Beat AML"),
  cox_rows(cox_fav3, "NPM1-mut / ITD-neg vs CBF/CEBPA + age", "Beat AML")))
ph_all <- c(ph_all, list(ph_rows(cox_npm1, "NPM1 + age", "Beat AML"),
                         ph_rows(cox_g4, "NPM1 x FLT3-ITD + age", "Beat AML")))
km_all <- c(km_all, list(km_rows(km_npm1, "Beat AML", "NPM1"),
                         km_rows(km_g4, "Beat AML", "NPM1 x FLT3-ITD")))

lab_hr <- function(fit, term) {
  r <- cox_rows(fit, "", "")
  r <- r[r$term == term, ]
  sprintf("HR %.2f (95%% CI %.2f to %.2f), p = %.2g", r$HR, r$CI_low, r$CI_high, r$p)
}
p1 <- km_plot(km_npm1, npm1_cols, "Beat AML main cohort: NPM1 status",
              "NPM1") +
  plot_annotation(subtitle = paste("Mutated vs wild-type, adjusted for age:",
                                   lab_hr(cox_npm1, "NPM1_statusmut")))
p2 <- km_plot(km_g4, grp4_cols, "Beat AML main cohort: NPM1 x FLT3-ITD", "Group") +
  plot_annotation(subtitle = paste("NPM1-mut / ITD-neg vs everyone else, adjusted for age:",
                                   lab_hr(cox_fav, "npm1_favNPM1-mut / ITD-neg")))
save_fig(wrap_elements(p1) | wrap_elements(p2), "10_km_npm1_flt3.png", w = 13, h = 6)

# ---- 3. NEK10/FTO/USP44-high vs -low within NPM1-mutated ---------------------
# The three-gene split from Session 8 (block_group: 37 high, 63 low). Cox with
# age and FLT3-ITD, because FLT3-ITD is the main known risk factor inside
# NPM1-mutated AML and could differ between the groups. The Mer et al.
# primitive / committed call (Session 5, pc_call) is run the same way as a
# report-only check, since that's the definition their survival result used.

block <- read.csv("results/block_nek10_fto_scores.csv")
pc    <- read.csv("results/subgroup_per_sample.csv")
mut <- ba_s %>% filter(NPM1 == "mut") %>%
  left_join(block[, c("dbgap_rnaseq_sample", "block_group", "core_score")],
            by = "dbgap_rnaseq_sample") %>%
  left_join(pc[, c("dbgap_rnaseq_sample", "pc_call")], by = "dbgap_rnaseq_sample") %>%
  mutate(block_group = factor(block_group, levels = c("low", "high")),
         pc_call = factor(pc_call, levels = c("committed-like", "primitive-like")))
stopifnot(nrow(mut) == 100, all(table(mut$block_group) == c(63, 37)))

km_blk  <- survfit(Surv(os_months, death) ~ block_group, data = mut)
lr_blk  <- survdiff(Surv(os_months, death) ~ block_group, data = mut)
cox_blk_u <- coxph(Surv(os_months, death) ~ block_group, data = mut)
cox_blk <- coxph(Surv(os_months, death) ~ block_group + age10 + FLT3_ITD, data = mut)
cat("\nNEK10/FTO/USP44 high vs low, log-rank p =",
    signif(1 - pchisq(lr_blk$chisq, 1), 3), "\n")
print(summary(cox_blk))
print(table(group = mut$block_group, died = mut$death))
cat("FLT3-ITD by group:\n"); print(table(mut$block_group, mut$FLT3_ITD))
# What else differs between the groups (helps read the unadjusted vs adjusted
# HR). Description only.
blk_chars <- mut %>% group_by(group = block_group) %>%
  summarise(cohort = "Beat AML", n = n(), median_age = median(ageAtDiagnosis),
            FLT3_ITD_pos = sum(FLT3_ITD == "positive"), intensive = sum(intensive),
            transplant = sum(transplant),
            CR_after_induction = sum(grepl("^Complete Response", responseToInductionTx)),
            deaths = sum(death),
            deaths_by_60_days = sum(death == 1 & os_days <= 60),
            ITD_pos_among_early_deaths = sum(death == 1 & os_days <= 60 &
                                             FLT3_ITD == "positive"),
            .groups = "drop")

# Report-only extras: continuous score, Mer et al. call, intensive-only subset.
cox_blk_cont <- coxph(Surv(os_months, death) ~ core_score + age10 + FLT3_ITD, data = mut)
cox_pc <- coxph(Surv(os_months, death) ~ pc_call + age10 + FLT3_ITD, data = mut)
mut_int <- mut %>% filter(intensive)
cox_blk_int <- coxph(Surv(os_months, death) ~ block_group + age10 + FLT3_ITD,
                     data = mut_int)

cox_all <- c(cox_all, list(
  cox_rows(cox_blk_u, "NPM1-mut: block high vs low (unadjusted)", "Beat AML"),
  cox_rows(cox_blk, "NPM1-mut: block high vs low + age + FLT3-ITD", "Beat AML"),
  cox_rows(cox_blk_cont, "NPM1-mut: core score (per unit of mean z)  + age + FLT3-ITD (report only)", "Beat AML"),
  cox_rows(cox_pc, "NPM1-mut: Mer primitive vs committed + age + FLT3-ITD (report only)", "Beat AML"),
  cox_rows(cox_blk_int, "NPM1-mut, intensive chemo only: block high vs low + age + FLT3-ITD (report only)", "Beat AML")))
ph_all <- c(ph_all, list(ph_rows(cox_blk, "NPM1-mut: block + age + FLT3-ITD", "Beat AML")))
km_all <- c(km_all, list(km_rows(km_blk, "Beat AML", "NPM1-mut: NEK10/FTO/USP44 group")))
pow_all <- c(pow_all, list(power_rows(sum(mut$death), mean(mut$block_group == "high"),
                                      "Beat AML", "NPM1-mut: block high vs low")))

# ---- 4. TCGA-LAML ------------------------------------------------------------
# The TCGA paper's clinical table on cBioPortal (study laml_tcga_pub) has
# OS_MONTHS ("overall survival in months since initial diagnosis") and
# OS_STATUS ("1:DECEASED" = event, "0:LIVING" = censored at last follow-up).
# Same table 08_tcga_replication.R caches; fetched here if it's missing.
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
tcga_clin <- read.csv(cbio_file)
tb <- read.csv("results/tcga_block_scores.csv")
tm <- tb %>%
  left_join(tcga_clin[, c("patientId", "OS_MONTHS", "OS_STATUS", "INDUCTION",
                          "TRANSPLANT_TYPE")], by = c("patient" = "patientId")) %>%
  mutate(os_months = as.numeric(OS_MONTHS),
         death = as.integer(OS_STATUS == "1:DECEASED"),
         core_group = factor(core_group, levels = c("low", "high")),
         FLT3_ITD = factor(FLT3_ITD, levels = c("negative", "positive")),
         age10 = age / 10)
stopifnot(nrow(tm) == 31, !anyNA(tm$os_months), !anyNA(tm$death))
cat("\nTCGA NPM1-mutated:", nrow(tm), "patients,", sum(tm$death), "deaths\n")
print(table(group = tm$core_group, died = tm$death))
fu_t <- survfit(Surv(os_months, 1 - death) ~ 1, data = tm)
cat("TCGA median follow-up (reverse KM):", round(summary(fu_t)$table["median"], 1),
    "months\n")

blk_chars <- bind_rows(blk_chars, tm %>% group_by(group = core_group) %>%
  summarise(cohort = "TCGA", n = n(), median_age = median(age),
            FLT3_ITD_pos = sum(FLT3_ITD == "positive"),
            intensive = sum(grepl("7\\+3", INDUCTION)),
            transplant = sum(!TRANSPLANT_TYPE %in% c("0", NA)),
            deaths = sum(death),
            deaths_by_60_days = sum(death == 1 & os_months * days_per_month <= 60),
            ITD_pos_among_early_deaths = sum(death == 1 & os_months * days_per_month <= 60 &
                                             FLT3_ITD == "positive"),
            .groups = "drop"))
write.csv(blk_chars, "results/survival_block_characteristics.csv", row.names = FALSE)
cat("\nNEK10/FTO/USP44 groups, characteristics (TCGA intensive = 7+3-based induction):\n")
print(as.data.frame(blk_chars))

km_t  <- survfit(Surv(os_months, death) ~ core_group, data = tm)
lr_t  <- survdiff(Surv(os_months, death) ~ core_group, data = tm)
cox_t_u <- coxph(Surv(os_months, death) ~ core_group, data = tm)
cox_t <- coxph(Surv(os_months, death) ~ core_group + age10 + FLT3_ITD, data = tm)
cat("TCGA high vs low, log-rank p =", signif(1 - pchisq(lr_t$chisq, 1), 3), "\n")
print(summary(cox_t))

cox_all <- c(cox_all, list(
  cox_rows(cox_t_u, "NPM1-mut: core high vs low (unadjusted)", "TCGA"),
  cox_rows(cox_t, "NPM1-mut: core high vs low + age + FLT3-ITD", "TCGA")))
ph_all <- c(ph_all, list(ph_rows(cox_t, "NPM1-mut: core + age + FLT3-ITD", "TCGA")))
km_all <- c(km_all, list(km_rows(km_t, "TCGA", "NPM1-mut: NEK10/FTO/USP44 group")))
pow_all <- c(pow_all, list(power_rows(sum(tm$death), mean(tm$core_group == "high"),
                                      "TCGA", "NPM1-mut: core high vs low")))

p3 <- km_plot(km_blk, block_cols, "Beat AML, NPM1-mutated (n = 100)",
              "NEK10/FTO/USP44") +
  plot_annotation(subtitle = paste("High vs low, adjusted for age and FLT3-ITD:\n",
                                   lab_hr(cox_blk, "block_grouphigh")))
p4 <- km_plot(km_t, block_cols, "TCGA-LAML, NPM1-mutated (n = 31)",
              "NEK10/FTO/USP44", xmax = 96) +
  plot_annotation(subtitle = paste("High vs low, adjusted for age and FLT3-ITD:\n",
                                   lab_hr(cox_t, "core_grouphigh")))
save_fig(wrap_elements(p3) | wrap_elements(p4), "10_km_nek10_fto_usp44.png", w = 13, h = 6)

# ---- Save tables -------------------------------------------------------------
cox_tab <- bind_rows(cox_all)
write.csv(cox_tab, "results/survival_cox.csv", row.names = FALSE)
write.csv(bind_rows(ph_all), "results/survival_ph_checks.csv", row.names = FALSE)
write.csv(bind_rows(km_all), "results/survival_km_summary.csv", row.names = FALSE)
pow <- bind_rows(pow_all)
write.csv(pow, "results/survival_power.csv", row.names = FALSE)
write.csv(ba_s %>% select(dbgap_rnaseq_sample, NPM1, FLT3_ITD, grp4, ageAtDiagnosis,
                          ELN2017, os_days, os_months, death, causeOfDeath,
                          intensive, transplant),
          "results/survival_beataml.csv", row.names = FALSE)
write.csv(tm %>% select(patient, core_group, FLT3_ITD, age, os_months, death,
                        INDUCTION, TRANSPLANT_TYPE),
          "results/survival_tcga.csv", row.names = FALSE)

cat("\nCox models:\n")
print(cox_tab %>% mutate(across(c(HR, CI_low, CI_high), ~ round(.x, 2)),
                         p = signif(p, 2)), row.names = FALSE)
cat("\nKaplan-Meier summaries:\n"); print(bind_rows(km_all), row.names = FALSE, digits = 3)
cat("\nPower:\n"); print(pow, row.names = FALSE, digits = 3)
cat("\nProportional hazards checks:\n"); print(bind_rows(ph_all), row.names = FALSE, digits = 3)

writeLines(capture.output(sessionInfo()), "results/sessionInfo_10_survival.txt")
cat("\nSaved results/sessionInfo_10_survival.txt\nDone.\n")
