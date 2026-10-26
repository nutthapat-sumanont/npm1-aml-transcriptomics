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
