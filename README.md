# NPM1-mutated vs wild-type AML (Beat AML 2.0)

A learning project: differential expression between NPM1-mutated and NPM1-wild-type AML using Beat AML 2.0 raw RNA-seq counts, connected to menin inhibitor biology. The session-by-session log, with all numbers and decisions, is in [LEARNING.md](LEARNING.md).

## Pipeline

Scripts run from the repo root, in order:

| Script | What it does |
|---|---|
| `scripts/01_download.R` | Downloads the Beat AML 2.0 files into `data/` (not committed) |
| `scripts/02_cohort.R` | De novo AML at diagnosis, one sample per patient, no APL or KMT2A fusions: 276 patients (100 NPM1-mutated, 176 wild-type) |
| `scripts/03_qc.R` | Filtering, library sizes, PCA, sex check |
| `scripts/04_deseq2.R` | DESeq2, `~ wave + specimen_type + sex + NPM1_status`, plus a FLT3-ITD sensitivity model |
| `scripts/05_enrichment.R` | GSEA (Hallmark, GO BP, Verhaak NPM1 controls), HOXA9/MEIS1 three-group plot, NEK10/FTO/USP44 check |
| `scripts/06_subgroup.R` | What the NEK10/FTO/USP44-high group within NPM1-mutated AML corresponds to |
| `scripts/06b_deseq2_subgroups.R` | Sensitivity DESeq2: each NPM1 subgroup vs wild-type (extra model; the main model is unchanged) |

`scripts/run_all.R` runs all of them in order.

## Main findings so far

- The positive controls hold: HOXA/HOXB, MEIS1 and PBX3 up, CD34 down, and the Verhaak NPM1 gene sets strongly enriched (NES +3.3 and -3.8).
- NPM1-mutated and KMT2A-rearranged AML reach the same high HOXA9 and MEIS1 levels, the shared menin dependency behind menin inhibitors.

![HOXA9 and MEIS1 in NPM1-mutated, KMT2A-rearranged and other AML](figures/05_three_group_hoxa9_meis1.png)

### A known subgroup inside NPM1-mutated AML

About 40% of NPM1-mutated samples (39 of 100) have high NEK10, FTO and USP44. This group **matches two published NPM1 subgroups**; it is not a new discovery:

- **Mer et al. 2021, "primitive" subtype.** Their primitive markers (CDH2, GPR12, MDFI, ZNF521, PLXNB1, HHIP) are all higher in our group and their committed markers (CD163, C1QA, CD14, MARCO, MSR1) all lower. 35 of 39 high samples look primitive and 55 of 61 low samples look committed.
- **Arza-Apalategi et al. 2023, IFEX9-high.** They showed that "high FTO" in about a third of NPM1-mutated AML is really an intragenic lncRNA, IFEX9, linked to FLT3-ITD, TET2, IDH2 and an immature profile. Our group is enriched for TET2 (32% vs 12%) and IDH2 (35% vs 15%), leans towards FLT3-ITD (51% vs 34%, not significant), and holds 10 of the 15 samples with NPM1 + FLT3-ITD + TET2 or IDH2.

**Caveat:** gene-level counts can't separate FTO from IFEX9. IFEX9 has no gene ID of its own in the Ensembl 75 annotation used by Beat AML, so its reads that overlap FTO exons are counted as FTO. Our "FTO up" may really be IFEX9.

![Primitive and committed scores by group](figures/06_subgroup_primitive_committed.png)

**Shared vs subgroup-specific genes.** In an extra DESeq2 model that compares each subgroup with wild-type, HOXA9, HOXA10, MEIS1 and PBX3 are up and CD34 down in both subgroups, so the HOX program belongs to NPM1-mutated AML as a whole. Of the 1,019 genes from the main model, 513 change about equally in both subgroups, 388 mostly in the primitive-like group (including NEK10, FTO and USP44) and 118 mostly in the committed-like group (rule: one subgroup's fold change at least 2x the other's). Scatter plot: `figures/06_subgroup_shared_vs_specific.png`.

## Data and citations

- Beat AML 2.0 data: https://github.com/biodev/beataml2.0_data (CC-BY-4.0). Bottomly D et al. Integrative analysis of drug response and clinical outcome in acute myeloid leukemia. *Cancer Cell* 2022. Tyner JW et al. Functional genomic landscape of acute myeloid leukaemia. *Nature* 2018.
- Arza-Apalategi S, Heuts BMH, Dooijes MTM, et al. High expression of an intragenic long noncoding RNA misinterpreted as high FTO oncogene expression in NPM1 mutant acute myeloid leukemia. *Leukemia* 2023;37:901-904. https://doi.org/10.1038/s41375-023-01844-w
- Mer AS, Heath EM, Madani Tonekaboni SA, et al. Biological and therapeutic implications of a unique subtype of NPM1 mutated AML. *Nature Communications* 2021;12:1054. https://doi.org/10.1038/s41467-021-21233-0
- Verhaak RGW et al. Mutations in nucleophosmin (NPM1) in acute myeloid leukemia (AML): association with other gene abnormalities and previously established gene expression signatures and their favorable prognostic significance. *Blood* 2005.
