# 13_classifier.R
# Session 13 (optional): can a small gene panel predict NPM1 status from
# expression, and does it hold up in a different cohort?
#
#  1. Beat AML expression: the 276 main-cohort samples, the same low-count rule
#     as 04_deseq2.R, one blind VST. Candidate genes are those with an Ensembl
#     match in TCGA-LAML (so the locked model can be applied there), minus the
#     male-specific chrY genes and XIST (sex is unbalanced between the groups,
#     the same reason they're removed before enrichment).
#  2. Nested cross-validation of a LASSO logistic regression (glmnet):
#     outer 5-fold CV, repeated 10 times, gives the performance estimate; an
#     inner 10-fold CV inside each outer training set picks lambda (lambda.1se).
#     The outer test fold is never seen while the model is fitted or tuned.
#     Also: how many genes does it need (models capped at 1, 2, 3, 5, 10, 20
#     genes), and how often each gene is picked across the 50 outer fits.
#  3. The final model: the same procedure on all 276 samples. Saved with the
#     VST parameters it needs, then locked before anything is read from TCGA.
#     The simplest model, HOXA9 + CD34 (plain logistic regression), is cross-
#     validated on the same outer folds and locked the same way.
#  4. One test in TCGA-LAML (111 main-cohort patients, 31 NPM1-mutated): TCGA
#     counts go through Beat AML's own VST (its size-factor reference and
#     dispersion trend), then each gene is put on Beat AML's scale (mean and SD
#     matched per gene, using no labels). The rescaling was fixed before TCGA
#     was looked at; "no rescaling" is reported next to it to show why it's needed.
#  5. Figure: cross-validated ROC in Beat AML and locked-model ROC in TCGA.
#
# Run from the repo root:  Rscript scripts/13_classifier.R   (about 5 minutes)
# Needs data/ from 01_download.R (Beat AML counts), data/tcga/ from
# 08_tcga_replication.R (TCGA counts), results/cohort.csv (02),
# results/tcga_cohort.csv (08) and results/chrY_genes_grch37.csv (05).

suppressPackageStartupMessages({
  library(DESeq2)
  library(glmnet)
  library(pROC)
  library(data.table)
  library(dplyr)
  library(ggplot2)
})

set.seed(2026)
dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

counts_file  <- "data/beataml_waves1to4_counts_dbgap.txt"
tcga_rds     <- "data/tcga/tcga_laml_star_counts.rds"
locked_file  <- "results/classifier_locked_model.csv"
n_outer      <- 5     # outer folds (performance estimate)
n_repeats    <- 10    # outer CV repeated with different fold splits
n_inner      <- 10    # inner folds (choosing lambda)
size_caps    <- c(1, 2, 3, 5, 10, 20)
simple_genes <- c("HOXA9", "CD34")
cohort_cols  <- c("Beat AML" = "#0072B2", "TCGA-LAML" = "#D55E00")

for (f in c(counts_file, tcga_rds)) if (!file.exists(f)) {
  stop("Missing ", f, ": run scripts/01_download.R and scripts/08_tcga_replication.R")
}

# Folds with the same NPM1-mutated / wild-type mix in each (stratified)
make_folds <- function(y, k) {
  folds <- integer(length(y))
  for (cl in unique(y)) {
    idx <- which(y == cl)
    folds[idx] <- sample(rep(seq_len(k), length.out = length(idx)))
  }
  folds
}

# ---- 1. Beat AML expression ---------------------------------------------------
# The VST is "blind": it uses only the counts, never the NPM1 labels, so fitting
# it once on all 276 samples doesn't leak the answer into the cross-validation.
# Its two learned parts (the size-factor reference and the dispersion trend) are
# kept so exactly the same transformation can be applied to TCGA later.

cohort <- read.csv("results/cohort.csv") %>% filter(main_cohort)
if (nrow(cohort) != 276) stop("Expected 276 main-cohort samples; re-run 02_cohort.R")
y <- as.integer(cohort$NPM1 == "mut")             # 1 = NPM1-mutated
cat("Beat AML main cohort:", length(y), "samples,", sum(y), "NPM1-mutated\n")

raw <- fread(counts_file, check.names = FALSE)
genes_beat <- as.data.frame(raw[, .(stable_id, display_label, biotype)])
counts <- as.matrix(raw[, cohort$dbgap_rnaseq_sample, with = FALSE])
rownames(counts) <- genes_beat$stable_id
storage.mode(counts) <- "integer"

# Same low-count rule as 04_deseq2.R: >= 10 reads in >= 100 samples
keep <- rowSums(counts >= 10) >= 100

# Genes with an Ensembl match in TCGA (GENCODE v36, version suffix dropped)
tc <- readRDS(tcga_rds)
tc_ids <- sub("\\..*$", "", tc$genes$gene_id[!grepl("_PAR_Y$", tc$genes$gene_id)])
in_tcga <- rownames(counts) %in% tc_ids

# Male-specific chrY genes (not PAR1/PAR2) and XIST
chrY <- read.csv("results/chrY_genes_grch37.csv")
sex_ids <- c(chrY$id[chrY$region == "MSY"],
             genes_beat$stable_id[genes_beat$display_label == "XIST"])
is_sex <- rownames(counts) %in% sex_ids

cand <- keep & in_tcga & !is_sex
gene_flow <- data.frame(
  step = c("All Beat AML genes", "Pass the low-count rule (04_deseq2.R)",
           "  ... and have an Ensembl match in TCGA", "  ... and are not chrY (male-specific) or XIST"),
  genes = c(nrow(counts), sum(keep), sum(keep & in_tcga), sum(cand)))
print(gene_flow, row.names = FALSE)
write.csv(gene_flow, "results/classifier_gene_flow.csv", row.names = FALSE)

# Blind VST, done by hand so its parts can be reused:
#   size factors = median-of-ratios against each gene's geometric mean,
#   dispersion trend = parametric fit of gene-wise dispersions on mean (design ~ 1)
cnt_beat <- counts[cand, ]
log_geo  <- rowMeans(log(cnt_beat))               # -Inf if any zero
geo_ref  <- exp(log_geo)                          # 0 = gene ignored for size factors
dds_b <- DESeqDataSetFromMatrix(cnt_beat, colData = data.frame(row.names = colnames(cnt_beat),
                                                              x = seq_len(ncol(cnt_beat))),
                                design = ~ 1)
dds_b <- estimateSizeFactors(dds_b, geoMeans = geo_ref)
dds_b <- estimateDispersionsGeneEst(dds_b, quiet = TRUE)
dds_b <- estimateDispersionsFit(dds_b, fitType = "parametric", quiet = TRUE)
X <- t(assay(varianceStabilizingTransformation(dds_b, blind = FALSE)))  # samples x genes
sym <- setNames(genes_beat$display_label, genes_beat$stable_id)
cat("VST matrix:", nrow(X), "samples x", ncol(X), "genes\n")

simple_ids <- names(sym)[match(simple_genes, sym)]
if (any(is.na(simple_ids)) || !all(simple_ids %in% colnames(X))) stop("HOXA9 or CD34 missing")

# ---- 2. Nested cross-validation -----------------------------------------------
# For each repeat r and outer fold f:
#   training set = the other 4 folds; test set = fold f (about 55 samples)
#   inner: cv.glmnet with 10 stratified folds of the training set picks lambda
#   (lambda.1se: the sparsest model within one standard error of the best
#   cross-validated deviance, the usual choice when the aim is a small panel)
#   then predict the test fold with that model.
# glmnet standardizes each gene inside the training set (standardize = TRUE),
# so the test fold doesn't influence the scaling either.

fit_lasso <- function(Xtr, ytr) {
  cv.glmnet(Xtr, ytr, family = "binomial", alpha = 1, nfolds = n_inner,
            foldid = make_folds(ytr, n_inner), type.measure = "deviance")
}
fit_simple <- function(Xtr, ytr) {
  d <- data.frame(y = ytr, Xtr[, simple_ids, drop = FALSE])
  names(d) <- c("y", simple_genes)
  glm(y ~ HOXA9 + CD34, family = binomial, data = d)
}
pred_simple <- function(m, Xte) {
  d <- as.data.frame(Xte[, simple_ids, drop = FALSE]); names(d) <- simple_genes
  as.numeric(predict(m, d, type = "response"))
}
auc_of <- function(yy, p) as.numeric(auc(roc(yy, p, levels = c(0, 1), direction = "<", quiet = TRUE)))

oof <- list()          # out-of-fold predictions
fits <- list()         # per outer fit: lambda, genes kept
cap_rows <- list()     # fixed-size models
t0 <- Sys.time()
for (r in seq_len(n_repeats)) {
  set.seed(2026 + r)
  outer <- make_folds(y, n_outer)
  for (f in seq_len(n_outer)) {
    tr <- outer != f; te <- outer == f
    set.seed(2026 + 100 * r + f)
    cvf <- fit_lasso(X[tr, ], y[tr])
    p_1se <- as.numeric(predict(cvf, X[te, ], s = "lambda.1se", type = "response"))
    p_min <- as.numeric(predict(cvf, X[te, ], s = "lambda.min", type = "response"))
    m_s   <- fit_simple(X[tr, ], y[tr])
    oof[[length(oof) + 1]] <- data.frame(
      repeat_id = r, fold = f, sample = rownames(X)[te], y = y[te],
      lasso = p_1se, lasso_min = p_min, simple = pred_simple(m_s, X[te, ]))
    b <- coef(cvf, s = "lambda.1se")[-1, 1]
    fits[[length(fits) + 1]] <- data.frame(
      repeat_id = r, fold = f, lambda_1se = cvf$lambda.1se, lambda_min = cvf$lambda.min,
      n_genes_1se = sum(b != 0),
      n_genes_min = sum(coef(cvf, s = "lambda.min")[-1, 1] != 0),
      genes_1se = paste(sym[names(b)[b != 0]], collapse = ";"))
    # "How many genes does it need?": along the same LASSO path, the least
    # penalized model with at most k genes. k fixes lambda, so no inner tuning.
    path <- cvf$glmnet.fit
    for (k in size_caps) {
      j <- max(which(path$df <= k))
      p_k <- as.numeric(predict(path, X[te, ], s = path$lambda[j], type = "response"))
      cap_rows[[length(cap_rows) + 1]] <- data.frame(repeat_id = r, fold = f, cap = k,
                                                     n_genes = path$df[j], sample = rownames(X)[te],
                                                     y = y[te], p = p_k)
    }
  }
  cat(sprintf("Repeat %d of %d done (%.1f min)\n", r, n_repeats,
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
oof  <- bind_rows(oof)
fits <- bind_rows(fits)
caps <- bind_rows(cap_rows)

# AUC per repeat, from that repeat's pooled out-of-fold predictions (every
# sample predicted exactly once, by a model that never saw it)
cv_auc <- oof %>% group_by(repeat_id) %>%
  summarise(lasso = auc_of(y, lasso), lasso_min = auc_of(y, lasso_min),
            simple = auc_of(y, simple), .groups = "drop")
cv_summary <- data.frame(
  model = c("LASSO (lambda.1se)", "LASSO (lambda.min)", "HOXA9 + CD34"),
  auc_mean = sapply(cv_auc[, -1], mean), auc_min = sapply(cv_auc[, -1], min),
  auc_max = sapply(cv_auc[, -1], max),
  genes_median = c(median(fits$n_genes_1se), median(fits$n_genes_min), 2),
  genes_min = c(min(fits$n_genes_1se), min(fits$n_genes_min), 2),
  genes_max = c(max(fits$n_genes_1se), max(fits$n_genes_min), 2))
cat("\nNested CV in Beat AML (", n_repeats, " x ", n_outer, "-fold):\n", sep = "")
print(cv_summary, row.names = FALSE, digits = 3)
write.csv(cv_auc, "results/classifier_cv_auc_by_repeat.csv", row.names = FALSE)
write.csv(fits, "results/classifier_cv_fits.csv", row.names = FALSE)
write.csv(oof, "results/classifier_cv_predictions.csv", row.names = FALSE)

cap_summary <- caps %>% group_by(cap, repeat_id) %>%
  summarise(auc = auc_of(y, p), genes = median(n_genes), .groups = "drop") %>%
  group_by(max_genes = cap) %>%
  summarise(auc_mean = mean(auc), auc_min = min(auc), auc_max = max(auc),
            genes_median = median(genes), .groups = "drop")
cat("\nAUC with at most k genes (nested CV):\n"); print(as.data.frame(cap_summary), digits = 3)
write.csv(cap_summary, "results/classifier_cv_by_size.csv", row.names = FALSE)

# How often each gene is picked across the 50 outer fits (stability)
picked <- unlist(strsplit(fits$genes_1se, ";"))
sel_freq <- as.data.frame(table(symbol = picked), stringsAsFactors = FALSE) %>%
  mutate(pct_of_fits = round(100 * Freq / nrow(fits), 1)) %>%
  rename(n_fits = Freq) %>% arrange(desc(n_fits))
write.csv(sel_freq, "results/classifier_cv_selection_frequency.csv", row.names = FALSE)
cat("\nGenes picked in at least half of the", nrow(fits), "outer fits:\n")
print(sel_freq[sel_freq$pct_of_fits >= 50, ], row.names = FALSE)

# ---- 3. Final model on all 276, then lock ---------------------------------------
# Same recipe (10-fold CV for lambda, lambda.1se) on every Beat AML sample.
# Coefficients are on the VST scale (glmnet returns them unstandardized);
# "per_sd" multiplies by the gene's Beat AML SD so genes can be compared.

set.seed(2026)
final <- fit_lasso(X, y)
b <- coef(final, s = "lambda.1se")[, 1]
kept <- names(b)[-1][b[-1] != 0]
de <- read.csv("results/deseq2_npm1.csv")
model_tab <- data.frame(
  term = c("(Intercept)", kept), gene_id = c(NA, kept), symbol = c(NA, sym[kept]),
  coef = unname(b[c("(Intercept)", kept)]),
  beat_mean = c(NA, colMeans(X[, kept, drop = FALSE])),
  beat_sd = c(NA, apply(X[, kept, drop = FALSE], 2, sd)))
model_tab$per_sd <- model_tab$coef * model_tab$beat_sd
model_tab$direction <- ifelse(is.na(model_tab$gene_id), NA,
                              ifelse(model_tab$coef > 0, "higher in NPM1-mut", "lower in NPM1-mut"))
model_tab$deseq2_rank <- de$rank[match(model_tab$gene_id, de$gene_id)]
model_tab$cv_pct_of_fits <- sel_freq$pct_of_fits[match(model_tab$symbol, sel_freq$symbol)]
model_tab <- model_tab[c(1, 1 + order(-abs(model_tab$per_sd[-1]))), ]
cat("\nFinal LASSO model:", length(kept), "genes (lambda.1se =",
    signif(final$lambda.1se, 3), ")\n")
print(model_tab[, c("symbol", "coef", "per_sd", "direction", "deseq2_rank", "cv_pct_of_fits")],
      row.names = FALSE, digits = 3)

# HOXA9 + CD34 on all 276
final_simple <- fit_simple(X, y)
simple_tab <- data.frame(term = names(coef(final_simple)),
                         gene_id = c(NA, simple_ids), symbol = c(NA, simple_genes),
                         coef = unname(coef(final_simple)),
                         beat_mean = c(NA, colMeans(X[, simple_ids])),
                         beat_sd = c(NA, apply(X[, simple_ids], 2, sd)))
simple_tab$per_sd <- simple_tab$coef * simple_tab$beat_sd
cat("\nHOXA9 + CD34 model:\n"); print(simple_tab[, c("term", "coef", "per_sd")], row.names = FALSE, digits = 3)

# Would the model use sex genes if they were allowed in? (A check only, same
# seed and recipe, on all candidate genes plus chrY and XIST.)
sex_in <- keep & in_tcga & is_sex
if (any(sex_in)) {
  dds_s <- DESeqDataSetFromMatrix(counts[keep & in_tcga, ], colData = colData(dds_b), design = ~ 1)
  dds_s <- estimateSizeFactors(dds_s, geoMeans = exp(rowMeans(log(counts[keep & in_tcga, ]))))
  dds_s <- estimateDispersionsFit(estimateDispersionsGeneEst(dds_s, quiet = TRUE),
                                  fitType = "parametric", quiet = TRUE)
  Xs <- t(assay(varianceStabilizingTransformation(dds_s, blind = FALSE)))
  set.seed(2026)
  fs <- fit_lasso(Xs, y)
  bs <- coef(fs, s = "lambda.1se")[-1, 1]
  sex_check <- data.frame(sex_genes_offered = sum(sex_in),
                          genes_kept = sum(bs != 0),
                          sex_genes_kept = paste(sym[names(bs)[bs != 0 & names(bs) %in% sex_ids]],
                                                 collapse = ";"))
  cat("\nCheck with chrY + XIST allowed:\n"); print(sex_check, row.names = FALSE)
  write.csv(sex_check, "results/classifier_sex_gene_check.csv", row.names = FALSE)
  rm(Xs, dds_s)
}

# Lock: everything needed to apply both models to new samples, written before
# TCGA expression is read.
locked <- bind_rows(cbind(model = "lasso", model_tab[, c("term", "gene_id", "symbol", "coef",
                                                         "beat_mean", "beat_sd")]),
                    cbind(model = "hoxa9_cd34", simple_tab[, c("term", "gene_id", "symbol", "coef",
                                                               "beat_mean", "beat_sd")]))
write.csv(locked, locked_file, row.names = FALSE)
write.csv(model_tab, "results/classifier_final_model.csv", row.names = FALSE)
saveRDS(list(geo_ref = geo_ref, dispersion_function = dispersionFunction(dds_b),
             genes = colnames(X)), "data/classifier_vst_reference.rds")
cat("\nLocked:", locked_file, "and data/classifier_vst_reference.rds\n")

writeLines(capture.output(sessionInfo()), "results/sessionInfo_13_classifier.txt")
