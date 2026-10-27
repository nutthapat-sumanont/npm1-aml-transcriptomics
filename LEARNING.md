# Learning log

## Session 1 (2026-10-26): setup, download, first look at the cohort

**What we did**
- Installed R 4.6.1 and all 12 packages from r2u, which serves CRAN and Bioconductor packages as ready-built Ubuntu binaries. Nothing compiles, so the whole install took about 1 minute. `setup/cloud_setup.sh` repeats it for new cloud sessions.
- `scripts/01_download.R` downloads six Beat AML 2.0 files into `data/` and checks none is a Git LFS pointer (a ~130 byte text stub that stands in for the real file).
- `scripts/02_cohort.R` (part 1) prints every clinical column. No filtering yet.

**Key numbers**
- Files: raw counts 117 MB, normalized expression 268 MB, drug AUCs 18 MB, mutations 3.5 MB, clinical 0.5 MB, sample mapping 0.07 MB.
- Clinical summary: 942 samples from 805 patients, 95 columns. One row is one sample, so some patients appear more than once.
- 698 samples have an RNA-seq ID, but only 671 are in the counts matrix (the rest failed RNA QC).

**What I should understand**
- Git LFS: big files on GitHub are replaced by small pointer files. Download the wrong way and you get the pointer, which still "downloads fine". Always check file size.
- The clinical `NPM1` column is the clinical lab call (usually PCR). It agrees with the exome calls in 862 of 903 samples; 39 are positive clinically but not in the exome data, which is expected because the NPM1 4-base insertion is hard for exome variant callers.

**Cohort (approved plan, `scripts/02_cohort.R` part 2)**

| Step | Filter | Samples | Patients | NPM1 mut / wt |
|---|---|---|---|---|
| 0 | All samples | 942 | 805 | 240 / 698 |
| 1 | RNA-seq in counts matrix | 671 | 613 | 174 / 493 |
| 2 | AML at sampling | 654 | 596 | 174 / 476 |
| 3 | Initial diagnosis | 445 | 443 | 125 / 320 |
| 4 | NPM1 known | 445 | 443 | 125 / 320 |
| 5 | One sample per patient | 443 | 443 | 125 / 318 |

- Final cohort: 443 patients, 125 NPM1-mutated and 318 wild-type. Saved in `results/cohort.csv`, with the step counts in `results/cohort_flow.csv`.
- Duplicate rule: keep the earliest sample, and bone marrow first if two are tied. Patient 2280 kept a bone marrow sample. Patient 2403 kept a leukapheresis sample, because it was collected a day before their bone marrow sample.
- Specimen type is unbalanced: bone marrow is 66% of mutated samples but 59% of wild-type. That's a reason to put specimen type in the DESeq2 design.

**Cohort refinement (steps 6 to 8, added after review)**

| Step | Filter | Samples | NPM1 mut / wt | % mutated |
|---|---|---|---|---|
| 5 | One sample per patient | 443 | 125 / 318 | 28.2 |
| 6 | De novo AML only | 317 | 102 / 215 | 32.2 |
| 7 | Drop APL | 301 | 102 / 199 | 33.9 |
| 8 | Drop KMT2A fusions | **277** | **100 / 177** | **36.1** |

Why each filter:
- **De novo only.** NPM1-mutated AML is almost always de novo, but secondary AML (transformed from MDS/MPN) and therapy-related AML have their own expression programs. When those cases sit mostly in the wild-type group, some genes look "NPM1-related" when they're really "secondary-AML-related". Beat AML marks a sample as not de novo if the patient had *any* earlier cancer, even an untreated one. The 1 sample with unknown status was also dropped.
- **No APL.** Acute promyelocytic leukemia (PML-RARA) is a biologically separate disease with its own treatment (ATRA and arsenic). All 18 APL cases were wild-type, so they would add a strong, unrelated signal to the wild-type side.
- **No KMT2A fusions in the main comparison.** KMT2A-rearranged AML depends on menin too, and it switches on the same HOXA/MEIS1 genes as NPM1-mutated AML. In the wild-type group, these cases would hide the very signal we want to find. They aren't thrown away: `cohort.csv` keeps them with `main_cohort = FALSE` and `group = "KMT2A_r"`, so later we can plot NPM1-mutated, KMT2A-rearranged and other AML side by side.

What the KMT2A flag means:
- `KMT2A_rearranged` comes from `consensusAMLFusions`, which Beat AML built by manually reviewing clinical cytogenetics plus RNA-seq fusion calls. So it means a **fusion** (a chromosome translocation): `MLLT3-KMT2A` is t(9;11), and `KMT2A_re` is a KMT2A rearrangement whose partner isn't named. It is **not** a KMT2A partial tandem duplication (PTD).
- KMT2A-PTD only shows up in free-text lab reports (as "MLL (PTD)" in `variantSummary`). It isn't tested consistently, so I didn't filter on it. 8 wild-type samples in the main cohort mention it.
- 2 of the 24 KMT2A-fusion cases also carry an NPM1 mutation, both t(9;11). That combination is rare. They sit in the `KMT2A_r` group, not `NPM1_mut`.

**Open questions for the DESeq2 step**
- Design formula: `~ specimenType + NPM1` at minimum. Should sex, age or sequencing cohort (wave) be added? Bone marrow is 67% of mutated samples but 58% of wild-type ones in the main cohort.
- Should the 8 wild-type samples with KMT2A-PTD in their notes stay in the comparison?

## Session 2 (2026-10-26): QC of the counts (`scripts/03_qc.R`)

**What we did**
- Built a DESeq2 dataset for the 277 main-cohort samples. The script pulls count columns in cohort order and stops with an error if `colnames(counts)` and the cohort sample IDs aren't identical. It also stops if any count is negative or not a whole number.
- Stored the planned design `~ wave + specimen_type + NPM1_status` (wt as reference) but didn't fit it. All QC is "blind" to the design.
- Library sizes, VST + PCA, a sex check, and the HOXA9 / MEIS1 / CD34 sanity check. Figures in `figures/03_*.png`, tables in `results/qc_*.csv`.

**Key numbers**
- Low-count filter: keep a gene if it has at least 10 reads in at least 100 samples (100 = size of the smaller group, NPM1-mutated). 63,677 genes in, **19,406 kept** (13,694 protein coding). All positive-control genes pass.
- Library size: median 40.5 M reads, range 20.0 to 74.9 M. **No outliers** (rule: more than 3 MADs from the median on log10 scale). Waves 3+4 run a bit deeper than Waves 1+2, which size factors handle.
- PCA (top 500 variable genes): PC1 19.7%, PC2 13.0%, PC3 7.8%, PC4 6.0%. No sample more than 4 SD out on PC1 or PC2.
  - R squared with NPM1 status: PC1 0.10, **PC2 0.43**. Specimen type and wave: under 0.05 on every PC. Sex: PC2 0.25, **PC4 0.55**.
- Sex check: 276 of 277 agree. **1 mismatch: BA2294R** (patient 2499, t(8;21) AML), recorded female but XIST is off and RPS4Y1 is fully on. That looks male. 5 recorded males have partly reduced RPS4Y1, which fits loss of chromosome Y in the leukemia cells (common in AML, especially t(8;21)), so they are not flagged.
- Sanity genes, median VST (roughly log2), mutated minus wild-type: **HOXA9 +2.7, MEIS1 +2.2, CD34 -7.3**. All three go the expected way, with very small p values.

**What I should understand**
- *Why filter low counts?* Genes with almost no reads can't show a real difference and they add to the multiple-testing burden. Tying the threshold to the smaller group size means a gene that's on only in NPM1-mutated cases still survives.
- *VST* puts counts on a log2-like scale where the noise is roughly the same for low and high genes. PCA on raw counts would be dominated by a handful of very highly expressed genes.
- *Reading the PCA:* NPM1-mutated samples sit together at the top of PC2. Wave and specimen type don't form their own clusters, which is good news: no big batch effect.
- *Sex is a confounder here.* 64% of NPM1-mutated patients are female vs 41% of wild-type. NPM1-mutated AML is known to be more common in women. Without sex in the model, X and Y genes (XIST, RPS4Y1, KDM5D, ...) will show up as "NPM1 genes" even though they are just about sex.

**Open questions (need your decision before Session 3, DESeq2)**
1. **Add sex to the design?** My suggestion: `~ wave + specimen_type + sex + NPM1_status`. Sex explains a large share of PC2 and PC4 and is unbalanced between groups.
2. **BA2294R:** drop it, keep it, or keep it but use the expression-based sex if sex goes in the design? I'd lean towards dropping it, because a sex mismatch can mean a sample mix-up, and then its NPM1 call could be wrong too.
3. **Wave "Both" level** has only 6 samples (3 wt, 3 mut). It works as a factor level, but we could ask whether those should be merged with another wave.
4. **Leukapheresis:** for the design I grouped the 10 leukapheresis samples with peripheral blood (`specimen_type = Blood`), since both are blood. The PCA plot still shows all three types separately.

**Deviations**
- Specimen type in the design is two levels (BM vs Blood), with leukapheresis counted as blood. CLAUDE.md just says "specimen_type". Reason: leukapheresis is blood, and 10 samples is too few for its own level.
- Wave uses the clinical `cohort` column as is, three levels (Waves1+2, Waves3+4, Both).
- I picked the low-count rule (10 reads in 100 samples) myself because the task asked me to choose and explain one. It's easy to change at the top of section 2 in the script.
- `dds` and `vsd` are saved as .rds files in `data/` (gitignored), not committed. Re-run `scripts/03_qc.R` to rebuild them (about 2 minutes).

**Decisions after the Session 2 QC (approved before Session 3)**
- Design: `~ wave + specimen_type + sex + NPM1_status`, wild-type as reference. Wave keeps three levels; leukapheresis counts as Blood.
- BA2294R is dropped. Main cohort is now **276** (100 mutated, 176 wild-type).
- Worth remembering: **BA3206R** (recorded female, NPM1 wild-type) has XIST about 13.5 but also partial RPS4Y1, about 10 on the VST scale (1,121 normalized reads; every other female is under 20). XIST says female, so it isn't a mismatch and it stays in. Possible explanations: some male cells in the sample (for example after a past transplant or transfusion, or a mixed sample), or a sex-chromosome abnormality. It turns out to matter for the Y genes in Session 3, see below.

## Session 3 (2026-10-26): differential expression (`scripts/04_deseq2.R`)

**What we did, step by step**
1. **Dropped BA2294R in `scripts/02_cohort.R`** as a new step 9, so `cohort.csv` is still the single place that says who is in the main cohort. The row stays in the file with `main_cohort = FALSE`, `group = "QC_excluded"` and a `qc_exclusion` reason. Nothing else in `cohort.csv` changed (checked with `all.equal`).
2. **Fitted DESeq2** with `~ wave + specimen_type + sex + NPM1_status` on the 276 samples, same low-count filter as Session 2 (10 reads in 100 samples, 19,398 genes kept). DESeq2 estimates sequencing depth per sample, a noise level (dispersion) per gene, then fits one negative binomial model per gene. The NPM1 coefficient is "mutated vs wild-type, holding wave, specimen type and sex fixed".
3. **Shrunk the fold changes with apeglm.** Raw fold changes for low-count, noisy genes can be huge just by chance. apeglm pulls those towards 0 and leaves well-measured genes almost unchanged. The p-values and the Wald statistic still come from the unshrunken fit; only the log2FC column is shrunken.
4. **Saved everything** to `results/deseq2_npm1.csv`, sorted by the absolute Wald statistic (same order as the p-value, but no ties when p rounds to 0). Extra columns: biotype, the unshrunken log2FC, rank, and up / down / ns.
5. **Figures:** `figures/04_volcano.png`, `figures/04_ma_plot.png`, `figures/04_heatmap_top50.png`.
6. **Positive-control and sex-gene checks**, then a **FLT3-ITD sensitivity run**.

**Key numbers**
- 6,713 genes at padj < 0.05. With padj < 0.05 **and** |log2FC| > 1: **1,019 genes, 307 up and 712 down** in NPM1-mutated.
- Top 10 by significance: NEK10 (up), TRH, TSPAN7, APP, SDK2, DUSP27, BAALC (down), FTO (up), CD34, LINC00865 (down). Of the top 50, 11 are up and 38 down (1 just misses the fold-change cut).
- **Positive controls all pass**, every one significant in the expected direction:

| Gene | log2FC (shrunken) | padj |
|---|---|---|
| HOXA9 | +1.60 | 5.3e-08 |
| HOXA10 | +1.51 | 8.9e-09 |
| HOXB2 to HOXB9 | +1.35 to +3.10 | 1.3e-05 to 2.5e-19 |
| MEIS1 | +1.60 | 6.8e-12 |
| PBX3 | +1.91 | 1.0e-42 (rank 16) |
| CD34 | -4.00 | 2.8e-50 (rank 9) |

  HOXB1 and HOXB13 didn't pass the low-count filter (barely expressed in blood cells, as expected).
- **Sex genes are not top hits.** XIST: log2FC -0.03, padj 0.88, rank 17,885 of 19,398. No sex gene in the top 50. But RPS4Y1 (log2FC -1.19, rank 1,517) and KDM5D (-1.04, rank 1,363) do pass the 1,019-gene cut-off, see the next section.
- **FLT3-ITD sensitivity:** **42 of the top 50 stay in the top 50**, and 49 of 50 stay significant. 866 genes are significant in both models (1,019 main, 938 with FLT3-ITD). Shrunken fold changes correlate at 0.98 across all genes. Genes that leave the top 50: POU4F1, CLIP3, HPGDS, AC096677.1, WDR86, KLRG2, SCG2, CD1C.

**What I should understand**
- *Why NPM1-mutated AML looks like this.* The "up" side is the HOX program (HOXA, HOXB, MEIS1, PBX3) that NPM1-mutant protein keeps switched on, with menin and KMT2A at those genes. That's exactly what menin inhibitors are meant to switch off. The "down" side is stem and progenitor genes (CD34, BAALC, MN1, SPARC): NPM1-mutated AML is usually CD34-negative. So more genes go down than up.
- *Reading the volcano plot:* x = how big the change is, y = how sure we are. With 276 samples, small changes become very significant, which is why we also ask for |log2FC| > 1 (a 2-fold change).
- *The heatmap* clusters samples using only the top 50 genes, and NPM1-mutated samples (orange bar) fall almost entirely into one branch. FLT3-ITD and sex don't line up with the branches.
- *Why RPS4Y1 still shows up as "down".* With sex in the model, DESeq2 fits **one** NPM1 effect for both sexes. For a Y gene, females have almost no reads, so any female with some reads moves the estimate a lot on a log scale. In males alone, RPS4Y1 barely differs (log2 ratio of medians -0.08). The "down" call comes from females: the wild-type group includes BA3206R with 1,121 reads, while mutated females average about 1. So this is an artefact of one sample, not NPM1 biology. KDM5D is the same story with a smaller push.
- *FLT3-ITD sensitivity:* FLT3-ITD is in 41% of mutated and 19% of wild-type cases here, so it could leak into the NPM1 result. 42 of 50 staying, and a 0.98 correlation, says the NPM1 signal mostly stands on its own.

**Open questions**
1. RPS4Y1 and KDM5D in the significant list are an artefact (see above). Options: leave them and note it; drop Y-chromosome genes (and XIST) from the results before enrichment; or drop BA3206R. My lean is to remove chrY genes plus XIST before gene-set enrichment, since a sex-adjusted model still can't handle near-zero counts in one sex well. Your call.
2. Next step: gene-set enrichment (fgsea with the Wald statistic as the ranking), checking that VERHAAK_AML_WITH_NPM1_MUTATED_UP is enriched, then the three-group HOXA9 / MEIS1 plot with KMT2A-rearranged cases.
3. NEK10 and FTO are the top "up" genes. Worth a quick literature look: are they known NPM1 targets, or something specific to this cohort?

**Deviations**
- BA2294R is removed in `02_cohort.R` (new step 9) rather than only inside `04_deseq2.R`, so `cohort.csv` stays the one source of truth. The Session 2 QC outputs (`results/qc_*.csv`, `figures/03_*.png`) were **not** re-run and still show 277 samples, which is right: they're the evidence for dropping it. Re-running `03_qc.R` now would use 276.
- `04_deseq2.R` rebuilds the dataset from the counts instead of loading `data/dds_filtered.rds` from Session 2, because that object had 277 samples and the old design.
- Wave "Both" has **5** samples in the main cohort (3 mutated, 2 wild-type), not 6 as the Session 2 notes said. BA2294R was in Waves1+2, so this didn't change.
- Extra outputs not asked for: `results/deseq2_npm1_summary.csv`, `results/deseq2_positive_controls.csv`, `results/deseq2_sex_genes.csv`, `results/deseq2_sensitivity_flt3_top50.csv` and `_summary.csv`. The sex-gene check also looks at six other Y genes (KDM5D, DDX3Y, UTY, EIF1AY, ZFY, USP9Y) as a wider net.
- The heatmap uses `vst(blind = FALSE)`, the usual choice after fitting; Session 2 QC used `blind = TRUE`.
- "Top 50" / "10 most significant" = ranked by absolute Wald statistic (same order as p-value, without ties at p = 0). The volcano plot floors padj at 1e-300 so nothing is lost off the top.

## Session 4 (2026-10-26): gene-set enrichment (`scripts/05_enrichment.R`)

**What we did, step by step**
1. **Removed sex-chromosome genes from the ranked list.** The DESeq2 table has no chromosome column, so the script asks the Ensembl GRCh37 server (same genome version as the Beat AML counts) for every gene on chrY once and saves the list to `results/chrY_genes_grch37.csv`. Later runs read the file and don't need the internet. Only the male-specific part of Y goes; the pseudoautosomal genes (CD99, CRLF2, IL3RA, IL9R, ...) are on both X and Y and expressed in both sexes, so they stay. XIST goes too.
2. **Preranked GSEA with fgsea.** Every remaining gene is ranked by the DESeq2 Wald statistic (log2FC divided by its standard error), so the top of the list is "clearly up in NPM1-mutated" and the bottom is "clearly down". Gene sets: MSigDB Hallmark (50 sets) and GO Biological Process (7,538 sets), MSigDB version 2026.1.Hs via msigdbr. Sets with 15 to 500 of our genes are tested. Genes are matched by Ensembl ID, which matches more genes than symbols do (our symbols date from 2014). Seed 2026.
3. **Positive controls:** VERHAAK_AML_WITH_NPM1_MUTATED_UP and _DN, each with an enrichment plot. The script stops with an error if either goes the wrong way.
4. **Dot plots** of the top 15 up and top 15 down pathways by NES, one for Hallmark and one for GO BP.
5. **Three-group plot** of HOXA9 and MEIS1 (VST) for NPM1-mutated, KMT2A-rearranged and other AML.
6. **NEK10 / FTO / USP44 / WDR86 check** inside the 100 NPM1-mutated samples. Report only, nothing filtered.

**Key numbers**
- **Sex-gene filter removed 26 of 19,398 genes** (25 male-specific chrY + XIST), leaving **19,372**. Only 2 of them were significant in Session 3: RPS4Y1 and KDM5D, the artefacts we expected. 25 PAR genes kept. Full list in `results/enrichment_removed_sex_genes.csv`.
- **Positive controls pass, very strongly:**

| Gene set | Genes | NES | FDR |
|---|---|---|---|
| VERHAAK_AML_WITH_NPM1_MUTATED_UP | 174 | +3.31 | 2.4e-36 |
| VERHAAK_AML_WITH_NPM1_MUTATED_DN | 237 | -3.76 | 1.4e-96 |

- **Hallmark:** 49 sets tested, **23 at FDR < 0.05, all 23 down** in NPM1-mutated. In fact all 49 have a negative NES. Strongest: E2F targets (NES -2.48), G2M checkpoint (-2.20), apical junction, interferon alpha and gamma response, heme metabolism, mitotic spindle, allograft rejection.
- **GO BP:** 3,540 sets tested, **468 at FDR < 0.05: 43 up, 425 down.**
  - Up: embryonic skeletal system morphogenesis (NES +2.57) and anterior/posterior pattern specification (+2.04). Their leading edges are HOXA2 to HOXA11, HOXB2 to HOXB7 and PBX3: **this is the HOX program showing up in GO's language.** A group of cilium terms (cilium organization, intraciliary transport) is also up.
  - Down: MHC class II antigen presentation (HLA-DRA, HLA-DRB1, HLA-DQ/DP, CD74), antigen receptor signalling, cell-cycle terms.
- **Three-group plot** (median VST, roughly log2):

| | NPM1-mutated (n=100) | KMT2A-r (n=24) | Other AML (n=176) |
|---|---|---|---|
| HOXA9 | 12.2 | 12.3 | 9.5 |
| MEIS1 | 12.2 | 12.1 | 10.1 |

  NPM1-mutated and KMT2A-r are each about 2 to 2.8 log2 units above other AML (Wilcoxon p < 1e-8). NPM1-mutated vs KMT2A-r: no difference (p = 0.42 and 0.40).
- **NEK10 / FTO / USP44 / WDR86 within NPM1-mutated:**
  - NEK10, FTO and USP44 rise and fall together (Spearman 0.79 to 0.83). WDR86 barely follows them (0.28 to 0.32). In wild-type samples the four are hardly correlated at all (-0.30 to 0.27), so the co-regulation is NPM1-mutated specific.
  - The pattern looks like an on/off switch: about 40 samples high, about 60 low (k-means split 39/61; bimodality coefficient 0.70, where above 0.555 hints at two groups).
  - **WDR86 is the FLT3-ITD gene**: median VST 8.6 with FLT3-ITD vs 4.8 without (p = 1.6e-08). **The NEK10 / FTO / USP44 core is not explained by FLT3-ITD** (core score p = 0.08; FTO p = 0.15, USP44 p = 0.14, NEK10 p = 0.03).
  - Nothing else in `cohort.csv` explains the split: sex p = 0.90, specimen type 0.18, centre 0.80, age 0.21. Wave (p = 0.04) and ELN2017 (p = 0.05) are borderline. That's from seven tests, and ELN2017 "Intermediate" in NPM1-mutated cases mostly means FLT3-ITD with a high allelic ratio, so it's the WDR86/FLT3 link showing through again.
  - Bottom line: there does seem to be a NEK10/FTO/USP44-high subgroup of NPM1-mutated AML (around 40%), and no clinical variable we have lines up with it.

**How to read NES and FDR**
- *Enrichment score (ES).* GSEA walks down the ranked list from most-up to most-down. Each time it meets a gene in the set, a running score goes up (more for genes with a bigger Wald statistic); each time it meets a gene not in the set, it goes down a bit. ES is the furthest the walk gets from zero. In the Verhaak UP plot, the green line shoots up within the first ~1,000 genes: most of the set sits near the top of our list.
- *NES (normalized enrichment score).* ES depends on the set's size, so fgsea divides it by the average ES of random gene sets of the same size. That puts all sets on one scale. **Sign = direction**: positive means the set's genes sit towards the "higher in NPM1-mutated" end, negative towards the "lower" end. **Size = how strongly**: around 1 is what random sets do; above about 1.5 (or below -1.5) is a clear shift; Verhaak's 3.3 and -3.8 are extreme.
- *p and FDR (padj).* The p-value asks how often a random set of that size would score at least this far out. FDR is the p-value adjusted for testing many sets (Benjamini-Hochberg), done separately within Hallmark and within GO BP. FDR < 0.05 means that among all sets we call significant, we expect under 5% to be false.
- *What GSEA does not tell you.* A significant set doesn't mean every gene in it changes. Look at the **leading edge** (the genes before the peak of the walk; the `leadingEdge` column in the tables). That's how we saw the "skeletal morphogenesis" hit is really just HOX genes. GO terms also overlap heavily, so 468 hits are maybe a few dozen real themes.

**What I should understand**
- *The positive controls pass with huge margins.* Verhaak et al. built their sets from an independent cohort on microarrays 20 years ago, and our RNA-seq ranking reproduces them. That's the strongest sign so far that the pipeline is working end to end.
- *Why Hallmark is all "down".* Partly biology: NPM1-mutated AML tends to be CD34-negative and HLA-DR-low, and both stem/progenitor and MHC class II genes fall here. But it's also a property of the gene sets. Hallmark genes as a group lean down in our ranking (median Wald statistic -0.55 vs about 0 for all genes), so the random-set baseline is centred near zero while almost every Hallmark set sits below it. Hallmark has no HOX set, so the main "up" story doesn't have a home there. GO BP and the Verhaak sets carry the up side.
- *Lower E2F / G2M (cell cycle) in NPM1-mutated* means a smaller share of cycling cells in the sample, on average. It doesn't mean NPM1-mutated AML is less aggressive. Bulk RNA-seq mixes in non-leukaemic cells, and CD34+ progenitors cycle a lot.
- *The three-group plot is the menin story in one picture.* NPM1-mutated and KMT2A-rearranged AML are different mutations, but both keep menin and KMT2A sitting on HOXA9 and MEIS1, and the two groups end up at the same high level. That shared dependency is why menin inhibitors (revumenib, ziftomenib) work in both. Other AML is spread out: some cases are also HOXA9-high without either lesion, which is worth remembering.
- *The two triangles* in the KMT2A-r box are the two cases with both a KMT2A fusion and an NPM1 mutation; both sit in the high range.

**Open questions**
1. **What is the NEK10/FTO/USP44-high subgroup?** Things we haven't checked: other mutations from the WES file (DNMT3A, IDH1/2, TET2, cohesin genes), NPM1 variant allele frequency (how much of the sample is leukaemia), and drug response in the Beat AML drug AUC data, including menin inhibitors if any were tested. Do you want this as a Session 5 step?
2. **GO BP is redundant.** Should we collapse overlapping terms (fgsea's `collapsePathways`, or clusterProfiler's `simplify`) for a cleaner portfolio figure?
3. **Cilium terms up** is unexpected for blood cells. Possibly genes that share HOX-linked regulation, or IFT/BBS genes with roles outside cilia. Worth a short literature look before calling it biology.
4. Next for menin biology: test published menin-inhibitor response signatures (genes that drop after revumenib or VTP50469 treatment) against our ranking. They aren't in MSigDB, so the gene lists would have to come from papers.

**Deviations**
- *"chrY genes" = male-specific region only.* The 25 pseudoautosomal genes stay, because they exist on both X and Y and aren't sex-specific. If we'd removed them as well, 50 genes would go instead of 26. CRLF2 and IL9R (PAR genes, both significant) are the ones that would matter.
- *chrY list from the Ensembl GRCh37 REST server*, saved to `results/chrY_genes_grch37.csv`. CLAUDE.md doesn't say where to get chromosome positions; `org.Hs.eg.db` only mapped 15,524 of our 19,398 IDs.
- *Gene sets matched by Ensembl ID*, not symbol (83% of Hallmark genes matched vs 81% by symbol, and no duplicate symbols to resolve).
- *Dot plots are not pre-filtered on FDR.* Hallmark has no set with NES > 0 at all, so an FDR filter would have left the "up" half empty without saying why. Filled dots are FDR < 0.05, open dots are not, and the Hallmark plot says no set has NES > 0. Every dot shown is in fact FDR < 0.05.
- *GSEA FDR corrected within each collection* (Hallmark, GO BP, Verhaak) rather than across all of them together.
- *Three-group VST* is a fresh blind VST on the 300 samples (276 main + 24 KMT2A-r), same low-count rule as before, not the Session 3 VST (which had only the 276).
- *Extra outputs not asked for:* `results/gsea_summary.csv`, `results/three_group_summary.csv`, `results/three_group_tests.csv`, `results/block_nek10_fto_*.csv` (scores, associations, per-gene FLT3-ITD test), and a per-gene FLT3-ITD test plus a three-gene core score, added after the plot showed WDR86 behaving differently.
- *`data/` was re-downloaded* with `scripts/01_download.R` (new container); nothing in it is committed.

## Session 5 (2026-10-26): step 1, what is the NEK10/FTO/USP44-high group? (`scripts/06_subgroup.R`)

Report only. Nothing filtered, DESeq2 model untouched. The high/low split is the one from Session 4 (39 high, 61 low), read from `results/block_nek10_fto_scores.csv`.

**What we did, step by step**
1. **Co-mutations.** Pulled DNMT3A, TET2, IDH1, IDH2, FLT3-TKD (D835/I836), NRAS, KRAS, PTPN11 and cohesin genes (SMC1A, SMC3, RAD21, STAG2) from the Beat AML exome calls; FLT3-ITD comes from the clinical table, since exomes can't call ITDs. 97 of the 100 NPM1-mutated samples have exome data; the other 3 are NA, not "wild-type". Fisher's test high vs low, BH-corrected across the 10 features.
2. **The triple combination** from Arza-Apalategi et al.: counted samples with NPM1 + FLT3-ITD + (TET2 or IDH2).
3. **Primitive vs committed** (Mer et al. 2021). Primitive score = mean z-score of CDH2, GPR12, MDFI, ZNF521, PLXNB1, HHIP; committed score = mean z-score of CD163, C1QA, CD14, MARCO, MSR1; both within the 100 NPM1-mutated samples. Same blind VST as Session 4 (all 11 markers pass the low-count filter anyway).
4. **What sits inside FTO in our annotation.** Asked the Ensembl GRCh37 server which genes overlap FTO, saved to `results/fto_locus_genes_grch37.csv`.

**Key numbers**
- **Co-mutations (high vs low, % of samples with exome data):**

| | Low (n=60) | High (n=37) | Fisher p | BH padj |
|---|---|---|---|---|
| FLT3-ITD (n=61 / 39) | 34% | 51% | 0.10 | 0.15 |
| TET2 | 12% | 32% | 0.018 | 0.044 |
| IDH2 | 15% | 35% | 0.027 | 0.053 |
| IDH1 | 12% | 14% | 0.76 | 0.85 |
| DNMT3A | 57% | 27% | 0.006 | 0.021 |
| PTPN11 | 25% | 0% | 0.0004 | 0.004 |
| NRAS | 20% | 0% | 0.003 | 0.015 |
| FLT3-TKD | 12% | 0% | 0.04 | 0.07 |
| cohesin | 17% | 5% | 0.12 | 0.15 |

- **NPM1 + FLT3-ITD + (TET2 or IDH2): 15 of 97** samples with exome data. **10 of 37 high (27%) vs 5 of 60 low (8%)**, Fisher p = 0.02. Of the 15: 5 with TET2, 10 with IDH2. NPM1 + (TET2 or IDH2) alone: 25 of 37 high vs 16 of 60 low (p = 1e-4).
- **Primitive vs committed lines up almost perfectly.** All 6 primitive markers are higher in the high group (by 3 to 9 VST units, every BH padj < 1e-7) and all 5 committed markers lower (by 3 to 4.6 units, padj < 1e-4). Primitive score vs the NEK10/FTO/USP44 core score: Spearman 0.91. A two-cluster split of primitive-minus-committed gives 35 of 39 high samples "primitive-like" and 55 of 61 low "committed-like" (Fisher p = 2e-16).
- **FTO locus:** in our annotation (Ensembl 75) the genes inside FTO are FTO itself, three antisense RNAs (RP11-357N13.1/2/3), FTO-IT1 and a pseudogene. None is called IFEX9. Only FTO-IT1 has enough reads to test, and it barely tracks FTO (Spearman 0.31).
- NPM1 variant allele frequency (exome): median 0.31 high vs 0.28 low (p = 0.025). Small difference; the high group isn't explained by lower leukaemia content.

**What I should understand**
- **This group isn't new. It matches two published NPM1 subgroups.**
  - *Arza-Apalategi et al., Leukemia 2023* found that "high FTO" in about a third of NPM1-mutated AML is really a long non-coding RNA called IFEX9 that sits inside the FTO gene. They linked it to FLT3-ITD, TET2, IDH2 and an immature profile. Our group is 39% of NPM1-mutated cases, is enriched for TET2 and IDH2, leans towards FLT3-ITD (51% vs 34%, not significant on its own), and holds two-thirds of the NPM1 + FLT3-ITD + TET2/IDH2 cases. That's the same picture.
  - *Mer et al., Nature Communications 2021* split NPM1-mutated AML into "primitive" (stem cell signature) and "committed" subtypes. Our high group is their primitive subtype, almost one to one.
  - Read together: the immature profile the Leukemia paper describes and the primitive subtype of the Nature Communications paper look like the same biology, and our NEK10/FTO/USP44 genes are markers of it in Beat AML.
- **Our "FTO up" may really be IFEX9.** Gene-level counts add up every read that lands in a gene's exons. IFEX9 doesn't have its own gene ID in our annotation (Ensembl 75), so any of its reads that overlap FTO exons get counted as FTO. With gene-level counts we can't tell FTO from IFEX9. Telling them apart would need transcript-level quantification (e.g. salmon against an annotation that includes IFEX9) or the raw reads, and Beat AML only shares counts openly. So "FTO is up in NPM1-mutated AML" in our DESeq2 table should be read as "the FTO locus is up", most likely driven by this subgroup.
- **The mutation pattern makes sense next to the primitive/committed idea.** The committed (low) group carries all the RAS-pathway mutations (NRAS, PTPN11) and the FLT3-TKD cases and more DNMT3A; the primitive (high) group carries TET2 and IDH2. These are associations from small numbers (37 vs 60), so they're a hint, not a finding.
- **Why the overlap with Mer et al. is so strong:** both scores come from the same RNA-seq, so some of the agreement is just "immature samples express immature genes". What makes it convincing is that the marker genes were picked by a different group in a different cohort, and they split our samples the same way.

**Open questions**
1. Drug response by subgroup (sorafenib, sunitinib, quizartinib, any menin inhibitor in the Beat AML drug AUC file). Mer et al. found the subtypes respond differently to kinase inhibitors. *Decision: planned for a later session, not now.*

### Session 5, sensitivity DESeq2 with the two NPM1 subgroups (`scripts/06b_deseq2_subgroups.R`)

**What we did.** An extra model next to the main one (which doesn't change). NPM1-mutated is split into the Session 4 groups, so `group` has three levels: wild-type (reference, n = 176), NPM1-high (primitive-like, n = 39) and NPM1-low (committed-like, n = 61). Design `~ wave + specimen_type + sex + group`, same low-count filter, so the 19,398 genes match the main model one to one. Each subgroup is compared with wild-type, with apeglm-shrunken fold changes.

**Key numbers**
- **(a) The positive controls hold in both subgroups.** All are padj < 1e-4 in the expected direction:

| Gene | NPM1-high vs wt | NPM1-low vs wt | Main model |
|---|---|---|---|
| HOXA9 | +1.79 | +1.37 | +1.60 |
| HOXA10 | +1.59 | +1.38 | +1.51 |
| MEIS1 | +1.98 | +1.23 | +1.60 |
| PBX3 | +2.06 | +1.80 | +1.91 |
| CD34 | -5.45 | -3.55 | -4.00 |

- **(b) The scatter** (`figures/06_subgroup_shared_vs_specific.png`): the HOXA and HOXB genes sit near the diagonal, so they go up about as much in both subgroups. NEK10, USP44 and FTO sit straight up the y-axis: up only in the high group, close to zero in the low group.
- **(c) Of the 1,019 main-model genes:**
  - *By fold change* (shared if neither subgroup's change is more than 2x the other's): **513 shared, 388 mostly high, 118 mostly low**, 0 going the wrong way in both.
  - *By significance* (padj < 0.05 in the main-model direction): 800 in both, 142 only in high, 70 only in low, 7 in neither.
  - The 142 "high only" genes are mostly up genes (116 of 142), even though the high group is smaller and has less power. That's the NEK10/FTO/USP44 module.

**What I should understand**
- *The core NPM1 signature is shared.* HOX/MEIS1/PBX3 up and CD34 down hold in both subgroups, so the menin-dependent HOX program is a feature of NPM1-mutated AML as a whole, not of one subgroup. That's what you'd expect if menin inhibitors work across NPM1-mutated AML.
- *About half the main-model list is shared and over a third leans to the primitive-like group.* Some of the "NPM1 genes" in our main table are really "primitive-subgroup genes" diluted across all 100 NPM1-mutated samples. NEK10 and FTO, our top two up genes, are the clearest examples: they're not NPM1-wide at all.
- *Two views of "shared" give different counts* because significance depends on group size (61 vs 39) and fold change doesn't. The fold-change view is fairer when the groups are different sizes; the significance view is shown so you can see the cut doesn't drive the story.
- *CD34 is lower in the primitive-like group than in the committed-like one* (-5.5 vs -3.6). Odd at first sight: "primitive" sounds like it should be CD34-high. But NPM1-mutated AML is typically CD34-negative, and Mer et al.'s primitive signature is about other stem-like genes (CDH2, ZNF521, HHIP...), not CD34. Worth checking against their paper before we write anything about it.

**Deviations**
- *Two rules for "shared vs driven by one"*, a fold-change ratio (2x) and a significance rule. You asked for the count without fixing a rule; both are in `results/deseq2_subgroups_shared_counts.csv`, and the figure uses the fold-change rule.
- *Extra outputs:* `results/deseq2_subgroups_controls.csv`, `results/deseq2_subgroups_shared_counts.csv`, `results/sessionInfo_06b_deseq2_subgroups.txt`.
- *Script numbered 06b* (first saved as 07, renamed later in the session because the next session uses 07), with the figure named `06_subgroup_shared_vs_specific.png` as you asked.
- *`scripts/run_all.R` added* (runs 01 to 06b in order, each in its own R process). Not test-run from a clean checkout, as agreed; a later session does that.

### Session 5, step 1 deviations
- *Extra mutations checked* beyond TET2/IDH2/FLT3-ITD: DNMT3A, IDH1, FLT3-TKD, NRAS, KRAS, PTPN11, cohesin, plus NPM1 VAF. These were on the Session 4 open-question list.
- *Mutation calls from the exome file only*, not the clinical `variantSummary` text. The exome sees NPM1 in only 81 of 97 samples (the clinical call is better for NPM1), so some TET2/IDH2 carriers may be missed too. The denominator for all mutation counts is the 97 with exome data.
- *Primitive/committed call* uses two-cluster k-means on primitive-minus-committed (seed 2026). Mer et al. used their own clustering on many more genes, so this is an approximation; the continuous scores and per-gene tests don't depend on it.
- *Marker genes forced through the low-count filter* in the VST; in practice all 11 pass anyway, so this changed nothing.
- *FTO-locus gene list* from the Ensembl GRCh37 REST server, saved to `results/fto_locus_genes_grch37.csv`.
- *`data/` re-downloaded* with `scripts/01_download.R` (new container).

### Session 5, rest of the session: rename, collapsed GO terms, README write-up

**What we did**
1. **Renamed** `scripts/07_deseq2_subgroups.R` to `scripts/06b_deseq2_subgroups.R` and its session info file to `results/sessionInfo_06b_deseq2_subgroups.txt`, because the next session uses 07. Updated every reference (run_all.R, CLAUDE.md, README.md, this file).
2. **Collapsed redundant GO BP terms** (new step 4b in `scripts/05_enrichment.R`). `fgsea::collapsePathways` takes the significant terms from most to least significant and asks, for each one, whether it's still enriched once you account for a more significant term it overlaps with. If not, it's folded into that "parent" term. New outputs: `results/gsea_gobp_collapsed.csv` (every FDR < 0.05 term, whether it's a main term, its parent and how many terms fold into it) and `figures/05_gsea_dotplot_gobp_collapsed.png`.
3. **Rewrote README.md as a short paper-style write-up.** Data, Methods, Results (eight figures, each with a one-paragraph caption) and How to reproduce are filled in. Summary, Background, Question and hypothesis, Interpretation, Limitations, Next steps and AI use are headings with a one-line note for you to write.
4. **Collected candidate limitations and next steps** from all sessions (list below).

**Key numbers**
- GO BP terms at FDR < 0.05: **468, down to 154 after collapsing** (18 up, 136 down).
- The collapsed "up" list still leads with the HOX terms (embryonic skeletal system morphogenesis, NES +2.57; anterior/posterior pattern specification) and the cilium group (intraciliary transport, protein localization to cilium, cilium organization). The "down" side is antigen presentation (MHC class II), antigen receptor signalling, cell adhesion and cell-cycle terms.
- Re-running `05_enrichment.R` gave identical GSEA tables (no CSV in `results/` changed), so the collapse step didn't disturb earlier results.

**What I should understand**
- *Collapsing doesn't change any result, only how many terms you have to read.* The 314 dropped terms are still significant; they just don't add anything beyond a parent term. The `parent` column in `results/gsea_gobp_collapsed.csv` shows where each one went.
- *154 terms is still a lot*, mostly because the "down" side covers many separate themes (immune, adhesion, cell cycle, ion transport). The dot plot shows the top 15 each way.
- *The cilium terms survive collapsing*, so they aren't just an echo of the HOX terms. Still an open question (Session 4, question 3).

**Candidate limitations (all sessions)**
- Bulk RNA-seq mixes leukaemic and non-leukaemic cells; blast percentage and purity aren't in the model, so some "down" pathways (cell cycle, MHC class II, immune) may reflect cell mix, not gene regulation.
- One cohort only; apart from the Verhaak gene sets, nothing is validated in an independent dataset.
- Observational data: associations, not cause and effect.
- Cohort limited to de novo AML at diagnosis; results may not carry over to secondary, therapy-related or relapsed AML.
- NPM1 status is the clinical consensus call; the exome detects NPM1 in only 81 of 97 NPM1-mutated samples, so the exome can't fully cross-check it.
- KMT2A-PTD isn't tested consistently in Beat AML; 8 wild-type samples with PTD in their notes stay in the main cohort.
- Sex is unbalanced (64% female in NPM1-mutated vs 41% in wild-type); chrY genes and XIST had to be removed before enrichment, and one female sample (BA3206R) has partial Y expression.
- One sample dropped as a possible swap (BA2294R); others could exist that the sex check can't catch.
- Batch: the wave "Both" level has only 5 samples.
- Old gene annotation (Ensembl 75, GRCh37, 2014 symbols); only gene-level counts are openly available, so transcripts (FTO vs IFEX9) can't be separated.
- GSEA: gene sets overlap; Hallmark sets lean "down" as a group in our ranking; FDR is corrected within each collection, not across all of them.
- The NEK10/FTO/USP44 subgroup was found in the same data it's tested in; the k-means split and the bimodality (coefficient 0.70) are modest evidence of two clean groups.
- The primitive/committed call uses 11 marker genes, an approximation of Mer et al.'s method.
- Co-mutation tests rest on small numbers (37 vs 60 with exome data) and several features; the exome may miss subclonal mutations; FLT3-ITD comes from a clinical assay, not the exome.
- Thresholds are choices: padj < 0.05 with |log2FC| > 1, the low-count filter (10 reads in 100 samples), the 2x rule for "shared vs mostly one subgroup".
- Age, ELN risk and blast percentage aren't in the DESeq2 model.
- The menin link is inferred from HOXA9/MEIS1 expression; no menin inhibitor response data have been analysed yet.
- No survival or drug response analysis yet.
- `scripts/run_all.R` hasn't been run from a clean checkout.

**Candidate next steps (all sessions)**
- Drug response by subgroup: sorafenib, sunitinib, quizartinib and any menin inhibitor in the Beat AML drug AUC file (planned for a later session).
- Test published menin-inhibitor response signatures (genes that fall after revumenib or VTP50469) against our ranking; the gene lists have to come from papers.
- Clean run of `scripts/run_all.R` from a fresh checkout.
- Replicate the NEK10/FTO/USP44 and primitive/committed split in an independent cohort (for example TCGA LAML or Leucegene).
- Separate IFEX9 from FTO with transcript-level data (needs raw reads, or a dataset that has IFEX9 in its annotation).
- Survival by subgroup (Mer et al. reported a survival difference).
- Check why CD34 is lower in the primitive-like group than in the committed-like group, against Mer et al.
- Look into the cilium terms: leading-edge genes and a short literature search.
- Cell-type deconvolution to separate cell mix from gene regulation in the "down" pathways.
- Adjust for blast percentage or NPM1 VAF as a purity proxy.
- Identify the HOXA9-high "other AML" cases (KMT2A-PTD, NUP98 fusions, other lesions?).
- Polish figures for the portfolio.

**Deviations**
- *`05_enrichment.R` was re-run* to add the collapse step. All CSV tables came out identical. Two committed files changed as a side effect: `figures/05_three_group_hoxa9_meis1.png` (the jitter positions moved, because the new step uses the random number generator before the plot; the data points are the same) and `results/sessionInfo_05_enrichment.txt` (one extra package loaded, curl).
- *Collapse settings:* `collapsePathways` on the 468 FDR < 0.05 terms, default threshold (p < 0.05), seed 2026. The uncollapsed dot plot (`figures/05_gsea_dotplot_gobp.png`) is kept.
- *README extras not asked for:* a table of the data files used, a cohort paragraph in Data, a References section, and a line in the CC-BY-4.0 attribution saying the data are used unchanged apart from our filtering and that the Beat AML team doesn't endorse this project. The Verhaak UP caption also gives the DN result.
- *Fixed a wrong citation in a code comment:* `05_enrichment.R` said Verhaak et al. 2005 was in *Haematologica*; it's *Blood* 2005;106:3747-3754. All five references in README.md were checked against PubMed.
- *`data/` re-downloaded* with `scripts/01_download.R` at the start of this session (new container); nothing in it is committed.

## Session 6 (2026-10-27): how much of the NPM1 signal is cell maturity? (`scripts/07_cell_composition.R`)

**Headline (sensitivity model, + HSC-like + Monocyte-like):** all five core positive controls survive at 97 to 117% of their size (MEIS1 is partly circular, see below), 49 of the top 50 genes stay significant (32 stay in the top 50), and fold changes correlate at 0.89 with the main model. Cell maturity explains only a small part of the NPM1 signal. The all-six model is an upper bound, not a result.

**What we did, step by step**
1. **Added `beataml2_manuscript_vg_cts.xlsx` to `scripts/01_download.R`.** This file is a normal git file in the Beat AML repo, not a Git LFS file, so the LFS media URL gives a 404. Added `raw.githubusercontent.com` as a third fallback URL. All 7 files download and none is an LFS pointer.
2. **Read how the scores were made** (Bottomly et al. 2022, STAR Methods, "Cell-type Scores") and fetched the gene lists (van Galen et al. 2019, Table S3), saved as `results/van_galen_top30_genes.csv`.
3. **Compared the six scores** between NPM1-mutated and wild-type, and between NEK10/FTO/USP44-high and -low NPM1-mutated samples. Figure `figures/07_celltype_scores.png`.
4. **Sensitivity DESeq2 with cell-type scores** (approved scores: HSC-like + Monocyte-like; plus all six as an upper bound), all on the same 270 samples, with the main design refit on those 270 for a fair comparison. Figure `figures/07_celltype_adjusted_lfc.png`.

**What the scores are**
- van Galen et al. sequenced single cells from AML bone marrow and trained a classifier that labels each malignant cell as one of six types along the normal blood hierarchy: **HSC-like** (stem cell), **Progenitor-like**, **GMP-like** (granulocyte-monocyte progenitor), **Promono-like** (promonocyte), **Monocyte-like** and **cDC-like** (conventional dendritic cell). For each type they listed the top 30 marker genes.
- Bottomly et al. turned each 30-gene list into one number per bulk sample: centre and scale each gene, separately for bone marrow and for blood samples, then take the first principal component of the 30 genes (like a WGCNA eigengene), flipped so that higher = higher average expression of the markers. Waves 1+2 were scored first and Waves 3+4 were projected onto the same PCs.
- So a score means "how much does this sample look like cell type X, compared with other Beat AML samples of the same specimen type". It is **not** a percentage of cells, and scores sit around 0 (roughly -10 to +9).
- Signature genes worth knowing: **MEIS1 is one of the 30 HSC-like genes** (and so is KMT2A). The cDC-like list is mostly HLA class II genes (HLA-DRA, HLA-DRB1, HLA-DQ/DP). The GMP-like list contains **TRH (our rank 2 gene) and POU4F1 (rank 29)**. HOXA9, HOXA10, PBX3 and CD34 are not in any list.

**Key numbers**
- File: 560 samples x 6 scores, exactly the samples with `manuscript_rnaseq = "yes"` in the clinical file.
- **270 of 276 main-cohort samples have scores: 97 of 100 NPM1-mutated, 173 of 176 wild-type** (38 of 39 high, 59 of 61 low). The 6 without: 5 are the wave "Both" samples (Bottomly scored none of those) and BA2644R (patient 2280; Bottomly used that patient's blood sample BA2011R, we kept the bone marrow one).
- **NPM1-mutated vs wild-type** (Glass's delta = mean difference in wild-type SDs; BH across six):

| Score | Median wt | Median mut | Glass d | padj |
|---|---|---|---|---|
| HSC-like | 0.36 | -0.69 | -0.75 | 4e-07 |
| Progenitor-like | 0.35 | 1.42 | +0.06 | 0.13 |
| GMP-like | 0.19 | -0.53 | -0.24 | 0.017 |
| Promono-like | -0.23 | 1.33 | +0.54 | 7e-05 |
| Monocyte-like | 1.17 | 0.45 | -0.02 | 0.70 |
| cDC-like | 0.88 | -1.60 | -0.87 | 1e-05 |

- **NPM1-high vs NPM1-low:** high has more Progenitor-like (d +0.66, padj 0.004) and much less Monocyte-like (d -1.22, padj 8e-07) and cDC-like (d -1.84, padj 2e-10); Promono-like lower (padj 0.015). HSC-like and GMP-like don't differ.
- The six scores are correlated (Progenitor vs Monocyte -0.80; GMP vs Promono +0.62); one PC explains 52% of their variance.
- **The six scores together predict NPM1 status with AUC 0.93** (logistic regression, in-sample). One at a time: 0.51 to 0.70. HSC-like + Monocyte-like: 0.75. GMP-like + Promono-like: 0.80.

**The sensitivity DESeq2 (steps 3 and 4)**

Three fits, all on the same 270 samples and the same 19,398 genes. Wave has two levels here because no wave "Both" sample has scores. Scores scaled to SD 1.
- `main_270`: the main design, `~ wave + specimen_type + sex + NPM1_status`. Shows what changes just from losing 6 samples (almost nothing: fold changes correlate 0.999 with the main model, 1,017 vs 1,019 significant genes).
- **`two_score`, the sensitivity model:** `~ wave + specimen_type + sex + HSC_like + Monocyte_like + NPM1_status`.
- `all_six`: all six scores. **This is an upper bound, not a result.** The six scores predict NPM1 status with AUC 0.93, so this model regresses out a lot of NPM1 biology along with cell maturity. It tells us the most that "cell type" could possibly explain, nothing more.

**Key numbers, sensitivity model (HSC-like + Monocyte-like)**
- **All five core positive controls survive** (padj < 0.05, expected direction), with almost no change in size. % = adjusted log2FC as a share of the log2FC from the main design on the same 270.

| Gene | Main (276) | Main design, 270 | **Two-score model** | % kept | Partly circular? |
|---|---|---|---|---|---|
| HOXA9 | +1.60 | +1.57 | **+1.66** (padj 3e-07) | 105% | no |
| HOXA10 | +1.51 | +1.47 | **+1.46** (5e-07) | 99% | no |
| MEIS1 | +1.60 | +1.62 | **+1.89** (6e-14) | 117% | **yes, MEIS1 is an HSC-like signature gene** |
| PBX3 | +1.91 | +1.90 | **+1.84** (2e-32) | 97% | no |
| CD34 | -4.00 | -4.02 | **-3.91** (2e-45) | 97% | no |

  All 8 HOXB genes that pass the low-count filter also survive. Of all the positive controls (HOXA9, HOXA10, HOXB cluster, MEIS1, PBX3, CD34), **MEIS1 is the only one inside a van Galen signature**, so it's the only one flagged as partly circular.
- **Top 50:** 49 of the main-model top 50 are still significant (padj < 0.05 and |log2FC| > 1), and all 50 go the same direction. The one that isn't, PHKA2, wasn't significant in the main model either (it just misses the fold-change cut). **32 of 50 stay in the top 50**; the rest move down the ranking but stay significant. Biggest movers: CLIP3 (rank 35 to 1,644), FGD5, CLEC9A, LRP6.
- **Fold changes correlate at 0.89** with the main model (Pearson; Spearman 0.82), and across the 1,019 main-model genes the adjusted fold change is on average **0.85x** the main-design one (slope through zero).
- 952 genes significant (285 up, 667 down), 703 shared with the main model's 1,019. Of the 316 that drop out, 252 still have padj < 0.05 but fall just below the 2-fold cut, so most aren't lost, just a bit smaller.
- NEK10 shrinks from +5.1 to +3.8 and FTO from +2.5 to +1.9; USP44 doesn't change (+3.9 to +3.8). All stay highly significant.

**Upper bound only (all six scores), not a result**
- The five core controls still all survive, but shrink: HOXA9 to 58% of its main-design size, HOXA10 52%, PBX3 82%, **CD34 70%** (-4.02 to -2.80). MEIS1 stays at 100% (partly circular: it's in the HSC-like score). TRH and POU4F1 are partly circular here too (both GMP-like genes).
- 44 of 50 top genes still significant (lost: FTO, MAN1A1, CLIP3, PHKA2, SCG2, CD1C); 21 of 50 stay in the top 50. Fold-change correlation 0.71, average fold change 0.54x. 677 significant genes.
- Read this as: even if we let cell type take everything it possibly could, including real NPM1 biology, the HOX/MEIS1/PBX3-up and CD34-down core is still there.

**What I should understand**
- *Bottom line: cell maturity explains only a small part of the NPM1 signal.* Adjusting for where a sample sits on the stem-to-monocyte axis keeps every positive control at roughly full size and keeps 49 of the top 50 genes significant. Fold changes shrink by about 15% on average. Cell maturity is a real influence, but not what drives the result.
- *NPM1-mutated AML looks less stem-like and more promonocyte-like*, with low cDC-like (HLA class II) signal. That matches what's known: NPM1-mutated AML is often CD34-negative and HLA-DR-low with monocytic features, and Bottomly et al. themselves report NPM1 enriched for Promono-like. Good sanity check.
- *The primitive-like (high) group isn't "more HSC-like".* It is Progenitor-like and lacks the mature monocyte/cDC signal. Committed-like (low) samples are the monocytic ones. So "primitive vs committed" in Mer et al. is mostly a progenitor-vs-monocytic split inside NPM1-mutated AML.
- *Why adjusting for cell type can also remove real NPM1 biology (CD34 is the clearest case).* Adjusting for a variable only makes sense if it's a confounder: something that affects both NPM1 status and gene expression from outside. Cell maturity isn't like that. The NPM1 mutation itself pushes the leukaemia towards a CD34-negative, more mature-looking state while keeping the HOX program on. So maturity is partly a **consequence** of NPM1 (a *mediator*). If you adjust for a mediator, you subtract the part of the NPM1 effect that works through it, and call what's left "the NPM1 effect". For CD34, low expression *is* the NPM1 phenotype: in the upper-bound model it shrinks from -4.0 to -2.8, and some of that 30% loss is very likely real NPM1 biology, not a correction. On top of that, the scores are built from the same RNA-seq, so they're themselves expression outcomes; a score that contains a gene adjusts that gene for itself (MEIS1 here, TRH and POU4F1 in the all-six model).
- *Why MEIS1 went up (+1.62 to +1.89) after adjusting.* MEIS1 is one of the 30 genes behind the HSC-like score, and NPM1-mutated samples have lower HSC-like scores. The model "corrects" for that lower score and so makes the NPM1 effect on MEIS1 look bigger. Because MEIS1 helps build the score, this number is partly circular and shouldn't be read as evidence either way. HOXA9 and PBX3, which aren't in any signature, are the cleaner evidence.
- *Why the all-six model is only an upper bound.* With AUC 0.93, the six scores are almost a stand-in for NPM1 status, so most of what NPM1 does to expression can be soaked up by them. Whatever survives that model is very robust; whatever doesn't could still be real NPM1 biology.
- *Why AUC matters.* If the scores together almost tell you who is NPM1-mutated, putting them all in the model leaves little NPM1 variation for DESeq2 to work with. That's why we picked two scores (AUC 0.75) for the actual sensitivity model.

**Open questions**
1. The two-score model is a sensitivity check; the main model stays as is. Should the README's Limitations section now say that cell maturity was checked (and how), rather than listing it as untested?
2. CLIP3 drops from rank 35 to 1,644 with only two scores added, and FGD5 and CLEC9A also fall. Worth a look at what they are (CLEC9A is a dendritic-cell gene, so it probably tracks the cDC/monocyte axis rather than NPM1).
3. A proper deconvolution (CIBERSORTx-style cell proportions, including non-leukaemic cells like T cells) would be a different and complementary way to ask the same question. The scores only describe the malignant hierarchy.
4. The 6 samples without scores (5 wave "Both", BA2644R) could be scored by projecting them onto Bottomly's PCs, but that needs the PC rotations, which aren't in the public files.

**Deviations**
- *Extra download in `07_cell_composition.R`*: van Galen et al. 2019 Table S3 (from the Cell supplementary file URL) is downloaded once to `data/` (not committed); the 180 genes are saved to `results/van_galen_top30_genes.csv`. You asked only for the Bottomly file; the gene lists were needed to check which positive controls sit inside a score. GMP-like uses the "tumor-derived, combined" column because Table S3 says "see combined signature" for the per-cell-type GMP-like column. 162 of the 180 symbols match our 2014 gene symbols; none of the 18 non-matching ones is a positive control.
- *`01_download.R` got a third fallback URL* (raw.githubusercontent.com), because the cell-type file isn't stored in Git LFS and the LFS media URL gives a 404 for it.
- *Sensitivity models use 270 samples, not 276*, and wave has two levels instead of three (no wave "Both" sample has scores). That's why there is an extra `main_270` fit, as agreed.
- *Same 19,398 genes as the main model* (the gene list is taken from `results/deseq2_npm1.csv`) rather than re-running the low-count filter on 270 samples, so every gene compares one to one.
- *"Survives"* = padj < 0.05 in the expected direction (the Session 3 rule). *"Still significant"* for the top 50 = padj < 0.05 and |log2FC| > 1.
- *Extra checks not asked for*: Glass's delta effect sizes (the measure Bottomly et al. used), score-vs-score correlations, score-vs-gene correlations and R squared (`results/celltype_score_gene_*.csv`), the AUC table (`results/celltype_score_predict_npm1.csv`), a fold-change slope, and the HOXB cluster in the controls table.
- *Extra outputs*: `results/celltype_score_coverage.csv`, `celltype_score_tests.csv`, `celltype_score_correlations.csv`, `celltype_signature_overlap.csv`, `deseq2_celltype_adjusted.csv`, `deseq2_celltype_controls.csv`, `deseq2_celltype_summary.csv`, `sessionInfo_07_cell_composition.txt`, and a second figure, `figures/07_celltype_adjusted_lfc.png`.
- *`scripts/run_all.R`* now includes 07. Still not test-run from a clean checkout.
- *`data/` re-downloaded* with `scripts/01_download.R` (new container); nothing in it is committed.
