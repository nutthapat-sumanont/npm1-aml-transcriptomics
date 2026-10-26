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
