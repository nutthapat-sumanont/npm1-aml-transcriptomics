#!/usr/bin/env bash
# Cloud setup for the NPM1 AML project (Ubuntu 24.04 "noble").
# Installs R plus every package the analysis needs as prebuilt Ubuntu binaries
# from r2u (CRAN + Bioconductor as .deb files), so nothing compiles from source.
# Typical run time: 2 to 4 minutes.
set -euo pipefail
export DEBIAN_FRONTEND=noninteractive

SUDO=""; [ "$(id -u)" -ne 0 ] && SUDO="sudo"

# Skip everything if a previous run already finished.
if command -v Rscript >/dev/null 2>&1 && Rscript -e 'suppressMessages({library(DESeq2);library(apeglm);library(fgsea);library(clusterProfiler);library(org.Hs.eg.db);library(msigdbr);library(pROC);library(patchwork);library(survival);library(recount3);library(rtracklayer);library(renv)})' >/dev/null 2>&1; then
  echo "R and packages already installed"; exit 0
fi

$SUDO apt-get update -qq
$SUDO apt-get install -y -qq --no-install-recommends ca-certificates gnupg wget

# r2u repository (binary CRAN/Bioconductor packages) and CRAN's R build
$SUDO wget -q -O /etc/apt/trusted.gpg.d/cranapt_key.asc \
  https://eddelbuettel.github.io/r2u/assets/dirk_eddelbuettel_key.asc
echo "deb [arch=amd64] https://r2u.stat.illinois.edu/ubuntu noble main" | \
  $SUDO tee /etc/apt/sources.list.d/cranapt.list >/dev/null
$SUDO wget -q -O /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc \
  https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc
echo "deb [arch=amd64] https://cloud.r-project.org/bin/linux/ubuntu noble-cran40/" | \
  $SUDO tee /etc/apt/sources.list.d/cran_r.list >/dev/null
# Prefer r2u packages over Ubuntu's older r-cran-* builds
printf 'Package: *\nPin: release o=CRAN-Apt Project\nPin: release l=CRAN-Apt Packages\nPin-Priority: 700\n' | \
  $SUDO tee /etc/apt/preferences.d/99cranapt >/dev/null

$SUDO apt-get update -qq
$SUDO apt-get install -y -qq --no-install-recommends \
  r-base-core r-base-dev \
  r-bioc-deseq2 r-bioc-apeglm r-bioc-fgsea r-bioc-clusterprofiler \
  r-bioc-org.hs.eg.db r-cran-msigdbr r-cran-ggplot2 r-cran-ggrepel \
  r-cran-pheatmap r-cran-readxl r-cran-data.table r-cran-dplyr \
  r-cran-proc r-cran-patchwork r-cran-httr r-cran-survival \
  r-bioc-recount3 r-bioc-rtracklayer r-cran-jsonlite r-cran-tidyr r-cran-renv

Rscript -e 'pk <- c("DESeq2","apeglm","fgsea","clusterProfiler","org.Hs.eg.db","msigdbr","ggplot2","ggrepel","pheatmap","readxl","data.table","dplyr","pROC","patchwork","httr","survival","recount3","rtracklayer","jsonlite","tidyr","renv"); for (p in pk) suppressPackageStartupMessages(library(p, character.only = TRUE)); cat("All", length(pk), "packages load OK\n")'
