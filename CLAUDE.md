# Project: NPM1-mutated vs wild-type AML (Beat AML 2.0)

## Goal
Compare gene expression between NPM1-mutated and NPM1-wild-type AML using Beat AML 2.0 raw counts, then connect the result to menin inhibitor biology. This is a learning and portfolio project. I'm a first-year HKU BSc student with basic R and Python.

## How to work with me
- Before each analysis step, explain in plain English what it does and why it's done this way.
- Stop and ask me before any scientific decision: sample filters, thresholds, the DESeq2 design formula. Commit and push before you stop to ask, because idle sessions can shut down.
- End every session by adding an entry to LEARNING.md: what we did, key numbers, what I should understand, open questions.
- Use R. Put code in numbered scripts in scripts/ (01_download.R, 02_cohort.R, ...) that run from the repo root.
- Never commit raw data. data/ is gitignored; scripts/01_download.R fetches it again.
- Commit small outputs (CSV tables in results/, 300 dpi PNG figures in figures/) so the next session can start from them.
- Set a random seed for anything random, and save sessionInfo() to results/.

## Data
Beat AML 2.0: https://github.com/biodev/beataml2.0_data (CC-BY-4.0).
Files are stored in Git LFS. Download with https://github.com/biodev/beataml2.0_data/raw/main/<filename> and check that each file is real data, not a small LFS pointer.
Cite: Bottomly et al., Cancer Cell 2022; Tyner et al., Nature 2018.

## Decisions so far (check against these every session)
- Cohort (done in Session 1, scripts/02_cohort.R): AML only, initial diagnosis, de novo only, one sample per patient (earliest), NPM1 from the clinical summary consensus call. Excluded APL and KMT2A fusions. Main cohort was 277 patients after Session 1 (100 NPM1-mutated, 177 wild-type); 276 after the Session 2 QC (see below). `main_cohort == TRUE` in results/cohort.csv.
- KMT2A fusion cases stay in cohort.csv (`group = "KMT2A_r"`) for a later three-group plot of HOXA9 and MEIS1: NPM1-mutated, KMT2A-rearranged, other AML.
- The 8 wild-type samples with KMT2A-PTD in their notes stay in the main cohort (PTD isn't tested consistently in Beat AML).
- DESeq2 design (approved after Session 2 QC): ~ wave + specimen_type + sex + NPM1_status, NPM1 wild-type as the reference level. Wave keeps three levels (Waves1+2, Waves3+4, Both). specimen_type is two levels, BM and Blood; leukapheresis counts as Blood.
- BA2294R dropped (recorded female, expression male; possible sample swap). Main cohort is now 276 patients: 100 NPM1-mutated, 176 wild-type. It stays in cohort.csv with `main_cohort = FALSE` and `group = "QC_excluded"` (step 9 in scripts/02_cohort.R).
- Sensitivity analysis: rerun with FLT3-ITD added (~ wave + specimen_type + sex + FLT3_ITD + NPM1_status) and compare top genes. The main model stays without FLT3-ITD.
- Before enrichment, remove chrY genes and XIST from the ranked list (approved after Session 3), because a sex-adjusted model still can't handle Y genes with near-zero counts in females (RPS4Y1 and KDM5D looked "down" because of one female sample, BA3206R). "chrY" means the male-specific region only; pseudoautosomal (PAR1/PAR2) genes are on X and Y in both sexes and stay. This removes 26 of 19,398 genes (25 chrY + XIST), leaving 19,372. The chrY list comes from Ensembl GRCh37 and is saved in results/chrY_genes_grch37.csv (scripts/05_enrichment.R step 1).
- Subgroup sensitivity DESeq2 (approved in Session 5, scripts/06b_deseq2_subgroups.R): an extra model only, the main model is unchanged. Three-level group from the Session 4 k-means split of the NEK10/FTO/USP44 score: NPM1 wild-type (reference), NPM1-high (primitive-like, n = 39), NPM1-low (committed-like, n = 61); design ~ wave + specimen_type + sex + group. Outputs results/deseq2_subgroups.csv and figures/06_subgroup_shared_vs_specific.png.
- The NEK10/FTO/USP44-high group matches published subgroups (Mer et al., Nat Commun 2021 "primitive"; Arza-Apalategi et al., Leukemia 2023, IFEX9-high). Describe it that way, not as a new finding. Gene-level counts can't separate FTO from IFEX9.
- Drug response by subgroup (sorafenib, sunitinib, quizartinib, any menin inhibitor) is planned for a later session, not before.
- Positive controls: HOXA9, HOXA10, HOXB cluster, MEIS1, PBX3 should be higher in NPM1-mutated; CD34 lower; VERHAAK_AML_WITH_NPM1_MUTATED_UP enriched. If they aren't, stop and look for a bug before going further.
- Work on main where possible, or open a pull request into main at the end of each session so the next session starts from the latest work.
- End every session with a "Deviations" note in LEARNING.md: anything done differently from this list, and why.
