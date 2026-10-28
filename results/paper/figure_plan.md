# Session 14: proposed figure layout for a short paper

Status: **proposal, waiting for approval.** Nothing is built yet.

Rules for the build (from the Session 14 prompt):
- Every panel is drawn from tables already in `results/`. No new analysis, no refits, no change to any existing file in `results/` or `figures/`.
- Where a panel has to redraw a curve from saved values (ROC curves from saved predictions, Kaplan-Meier curves from the saved survival table), the script recomputes the headline number and **stops with an error if it doesn't match LEARNING.md** (for example AUC 0.992, HR 0.41). So the figures can't drift from the reported results.
- One style everywhere: one font family (sans, Helvetica/Arial metrics), base size 8 pt, panel letters A, B, C in bold 10 pt, 180 mm wide (two-column page width), each figure saved as PDF and 300 dpi PNG in `figures/paper/`.

## Colours

| What | Colour | Where it comes from |
|---|---|---|
| NPM1-mutated | `#D55E00` (vermillion) | `npm1_cols`, `group_cols` in 03 to 10 |
| KMT2A-rearranged | `#0072B2` (blue) | `group_cols` in 05, 08 |
| Other AML / NPM1 wild-type | `grey55` | `npm1_cols`, `group_cols` |
| Up / down genes | `#D55E00` / `#0072B2`, non-significant `grey80` | `dir_cols` in 04, 07 |

Two choices I need you to make (questions 1 and 2 at the end): the colours for the two NPM1 subgroups, and how to tell Beat AML from TCGA on the same panel.

## Figure 1. A HOX-centred expression signature in NPM1-mutated AML (Beat AML, n = 276)

**Point:** NPM1 status is the largest structured signal in the data, and its signature (HOX/MEIS1/PBX3 up, CD34 down) reproduces the published Verhaak signature and is shared with KMT2A-rearranged AML.

| Panel | Content | Source table | Numbers it shows (LEARNING.md) |
|---|---|---|---|
| A | PCA, PC1 vs PC2, coloured by NPM1 status | `qc_pca_scores.csv` | NPM1 explains 43% of PC2 |
| B | Volcano plot, positive controls labelled (HOXA9, HOXA10, MEIS1, PBX3, CD34, plus NEK10 and FTO as the top up genes) | `deseq2_npm1.csv` | 1,019 genes: 307 up, 712 down |
| C | Gene-set enrichment, NES bars: Verhaak UP and DN, then the top 5 up and top 5 down *main* collapsed GO BP terms by NES; FDR printed next to each bar | `gsea_verhaak_controls.csv`, `gsea_gobp_collapsed.csv` | Verhaak UP +3.31, DN -3.76 |
| D | HOXA9 and MEIS1 in NPM1-mutated, KMT2A-rearranged and other AML (box + points) | `three_group_hoxa9_meis1.csv`, `three_group_tests.csv` | HOXA9 median 12.2 vs 12.3 vs 9.5; NPM1 vs KMT2A-r p = 0.42 |

Note on A: that PCA was computed in Session 2 on the 277 samples before QC. I'd plot the 276 main-cohort samples from it and say in the legend that BA2294R isn't shown. The alternative is to plot all 277 with BA2294R marked; see question 3.

## Figure 2. The signature is robust to cell maturity and replicates in TCGA-LAML

**Point:** adjusting for cell maturity barely changes the NPM1 signal, and both the fold changes and a locked Beat AML classifier carry over to an independent cohort (TCGA-LAML, n = 111).

| Panel | Content | Source table | Numbers it shows |
|---|---|---|---|
| A | log2FC, main design on 270 samples vs the two-score model (+ HSC-like + Monocyte-like), the 1,019 significant genes highlighted, HOXA9 / PBX3 / CD34 labelled (MEIS1 marked as partly circular) | `deseq2_celltype_adjusted.csv` | r = 0.89; average 0.85x; 49 of top 50 still significant |
| B | log2FC Beat AML vs TCGA-LAML, all shared genes, Beat AML significant genes highlighted | `deseq2_npm1.csv`, `tcga_deseq2_npm1.csv` (joined on `gene_id`) | r = 0.74 (all), 0.84 (significant); 95.4% same direction |
| C | TCGA signature score (Beat AML top 50 up minus top 50 down) by group: NPM1-mutated, KMT2A-r, other AML | `tcga_signature_scores.csv` | ROC AUC 0.992 (0.981 to 1.000) |
| D | Classifier ROC: LASSO (33 genes) and HOXA9 + CD34, nested CV in Beat AML and the locked test in TCGA | `classifier_cv_predictions.csv`, `classifier_tcga_predictions.csv`, `classifier_cv_summary.csv`, `classifier_tcga_performance.csv` | CV 0.985 / 0.950; TCGA 0.996 / 0.966; 108 vs 103 of 111 correct |

## Figure 3. A primitive-like NEK10/FTO/USP44-high subgroup inside NPM1-mutated AML

**Point:** within NPM1-mutated AML a primitive-like subgroup (37 of 100) switches on NEK10, USP44 and the IFEX9 transcript inside FTO on top of the shared HOX programme, and carries a different co-mutation pattern; it matches the published Mer et al. "primitive" and Arza-Apalategi et al. IFEX9-high groups.

| Panel | Content | Source table | Numbers it shows |
|---|---|---|---|
| A | The 100 NPM1-mutated samples ordered by core score: NEK10, FTO, USP44 (VST) as a heat strip, the k-means split marked, and tracks underneath for primitive/committed call, FLT3-ITD, TET2, IDH2, DNMT3A, NRAS, PTPN11 | `block_nek10_fto_scores.csv`, `subgroup_per_sample.csv`, `subgroup_comutations.csv` | 37 high / 63 low; 32 of 37 primitive-like; TET2 33% vs 11%, IDH2 39% vs 13%, NRAS and PTPN11 0% vs 20% and 25% |
| B | Primitive score vs committed score, coloured by subgroup | `subgroup_per_sample.csv` | primitive score vs core score Spearman 0.91 |
| C | Shared vs subgroup-specific: log2FC high vs wild-type against low vs wild-type for the 1,019 genes, coloured by the 2x rule, HOX genes and NEK10 / FTO / USP44 labelled | `deseq2_subgroups.csv`, `deseq2_subgroups_shared_counts.csv` | 521 shared, 381 mostly high, 117 mostly low |
| D | FTO locus in TCGA-LAML, one dot per patient: coverage of FTO exons 1-8, FTO exon 9, IFEX9 exon 1, and IFEX9-to-exon-9 junction reads, in FTO-high (6), FTO-low (25) and wild-type (80) | `fto_locus_per_sample.csv` | IFEX9 exon 1 median 256 vs 1.1; exons 1-8 not higher; junction reads 61 to 434 in all 6 high |

Note on D: the coverage track itself (`figures/11_fto_coverage.png`) was drawn from bigWig files in `data/` and can't be rebuilt from `results/`. The per-sample summary makes the same argument with committed numbers. If you'd rather have the track, I can place the existing PNG as a panel, but then its font and size won't match the others.

## Figure 4. Exploratory drug response and survival

**Point:** NPM1 status on its own barely changes ex vivo drug response or (once FLT3-ITD is split out) survival; the subgroup does split drug response (venetoclax vs MEK/HDAC inhibitors), while its survival difference is inconclusive.

| Panel | Content | Source table | Numbers it shows |
|---|---|---|---|
| A | Drug volcano, NPM1-mutated vs wild-type, AUC difference adjusted for monocyte-like score, coloured by drug class | `drug_npm1_lm.csv` | 0 of 123 drugs at FDR < 0.05 |
| B | Forest plot, primitive-like minus committed-like AUC difference with 95% CI: the 15 drugs at FDR < 0.05 plus sorafenib, sunitinib and quizartinib (Mer et al.) | `drug_subgroup_npm1.csv` | venetoclax -100 (FDR 8e-6); sunitinib FDR 0.32, sorafenib 0.41, quizartinib 0.99 |
| C | Kaplan-Meier, four NPM1 x FLT3-ITD groups, with numbers at risk | `survival_beataml.csv`, `survival_cox.csv` | HR 0.41 (0.25 to 0.67) for NPM1-mut / ITD-neg vs everyone else |
| D | Kaplan-Meier within NPM1-mutated, high vs low, with numbers at risk | `survival_beataml.csv` joined to `block_nek10_fto_scores.csv`, `survival_cox.csv`, `survival_power.csv` | HR 1.19 (0.64 to 2.21); 25% power for HR 1.5 |

Note on venetoclax: the per-sample venetoclax AUCs were never written to `results/` (they came straight from `data/`), so the venetoclax scatter from `figures/09_venetoclax.png` can't be redrawn without rerunning 09. Panel B shows the same result as an effect size with its CI.

## Supplementary tables (`results/paper/supplementary_tables.xlsx`)

Sheet 1, `README`: one row per table (sheet name, what it is, source file, session, rows), then one row per column (sheet, column, meaning, units). The script stops if any column lacks a description. Tables are copied from `results/` unchanged apart from rounding for readability (p and padj kept as is).

| Sheet | Content | Source in `results/` |
|---|---|---|
| S1_cohort_beataml | Cohort flow, Beat AML | `cohort_flow.csv` |
| S2_cohort_tcga | Cohort flow, TCGA-LAML | `tcga_cohort_flow.csv` |
| S3_deseq2_npm1 | Main DESeq2 result, all 19,398 genes | `deseq2_npm1.csv` |
| S4_positive_controls | Positive controls in every model (main, TCGA, two-score, subgroups) | `deseq2_positive_controls.csv`, `tcga_positive_controls.csv`, `deseq2_celltype_controls.csv`, `deseq2_subgroups_controls.csv` |
| S5_gsea | Verhaak, Hallmark, collapsed GO BP (FDR < 0.05) | `gsea_verhaak_controls.csv`, `gsea_hallmark.csv`, `gsea_gobp_collapsed.csv` |
| S6_sensitivity | FLT3-ITD and cell-maturity sensitivity summaries | `deseq2_sensitivity_flt3_summary.csv`, `deseq2_celltype_summary.csv` |
| S7_tcga_replication | Beat AML vs TCGA summary, three-group summary and tests | `tcga_vs_beataml_summary.csv`, `tcga_three_group_summary.csv`, `tcga_three_group_tests.csv` |
| S8_subgroup | Co-mutations, cell-type score tests, shared vs specific counts | `subgroup_comutations.csv`, `celltype_score_tests.csv`, `deseq2_subgroups_shared_counts.csv` |
| S9_fto_ifex9 | FTO locus key feature tests (TCGA) | `fto_key_feature_tests.csv` |
| S10_drugs | Per-drug models, NPM1 and subgroup | `drug_npm1_lm.csv`, `drug_subgroup_npm1.csv` |
| S11_survival | Cox models, Kaplan-Meier summary, power | `survival_cox.csv`, `survival_km_summary.csv`, `survival_power.csv` |
| S12_classifier | Final model genes, CV summary, CV by size, TCGA performance | `classifier_final_model.csv`, `classifier_cv_summary.csv`, `classifier_cv_by_size.csv`, `classifier_tcga_performance.csv` |

Sheets that combine several source files stack them with a `table` column saying which file each row came from, or sit side by side with a gap, depending on whether they share columns. No per-patient table goes in the workbook (same reasoning as the gene lookup in Session 11).

Possible extra sheet if you want it: TCGA DESeq2 for all 26,069 genes (`tcga_deseq2_npm1.csv`). It makes the file about twice as big.

## Package

One new package: **openxlsx** (writes .xlsx with bold headers, column widths and a frozen top row; no Java). From r2u as `r-cran-openxlsx`, added to `setup/cloud_setup.sh`, `setup/renv_snapshot.R` and `renv.lock`. Everything else (ggplot2, patchwork, ggrepel, survival, pROC) is already locked.

## Questions before I build

1. **Subgroup colours.** Both subgroups are NPM1-mutated, so I'd suggest two shades of the NPM1 orange: primitive-like (high) dark `#A63D00`, committed-like (low) light `#F4A582`, so "orange = NPM1-mutated" holds across all four figures. The alternative is the colours already used in 07, 09 and 10 (high `#D55E00`, low `#009E73` green), which match the existing figures but make orange mean two things.
2. **Beat AML vs TCGA on one panel** (Figure 2D). Session 13 used blue and orange for the two cohorts, which would clash with KMT2A-r and NPM1 here. I'd suggest black solid lines for Beat AML and black dashed for TCGA, with the two models told apart by grey vs black.
3. **PCA (Figure 1A):** plot the 276 main-cohort samples (BA2294R left out, said in the legend), or all 277 with BA2294R marked?
4. **Anything to swap?** For example the volcano for the heatmap, or Figure 4 dropping survival for the HOX-score drug result.
