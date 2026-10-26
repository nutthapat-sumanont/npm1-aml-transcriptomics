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

**Open questions for the DESeq2 step**
- De novo only, or all AML? The wild-type group is 32% secondary or therapy-related AML, but the mutated group is only 18%.
- KMT2A-rearranged cases: 31 are wild-type and 3 are NPM1-mutated, which is unusual and worth checking. They're menin-dependent, so leaving them in the wild-type group weakens the signal we care about.
- APL: 18 cases, all wild-type. Exclude them or keep them?
- Design formula: `~ specimenType + NPM1` at minimum. Should sex, age or sequencing cohort (wave) be added?
