#!/bin/bash
###############################################################################
# capture_versions.sh - record the EXACT tool versions installed on this machine.
#
# TOOLS_VERSIONS.txt lists the versions we STANDARDIZE to. This script prints what is
# ACTUALLY installed right now, so you can verify a cluster module stack or a freshly
# created conda env matches. Run it inside whatever env you want to audit:
#
#     conda activate shakira && bash refs/capture_versions.sh > versions_actual.txt
#   or on the module stack:
#     module load BCFtools PLINK/1.90-foss-2016a R/4.2.0-foss-2022b ; bash refs/capture_versions.sh
###############################################################################
set -uo pipefail
line(){ printf '%-14s %s\n' "$1" "$2"; }
ver(){  command -v "$1" >/dev/null 2>&1 && line "$1" "$($2 2>&1 | head -1)" || line "$1" "NOT FOUND"; }

echo "# tool versions captured $(date -u +%Y-%m-%dT%H:%M:%SZ) on $(hostname)"
ver bcftools   "bcftools --version"
ver samtools   "samtools --version"
ver tabix      "tabix --version"
ver bgzip      "bgzip --version"
ver bedtools   "bedtools --version"
ver plink      "plink --version"
ver plink2     "plink2 --version"
ver admixture  "admixture --version"
ver parallel   "parallel --version"
ver R          "R --version"
ver Rscript    "Rscript --version"
ver pandoc     "pandoc --version"
ver python3    "python3 --version"
ver gatk       "gatk --version"
ver nextflow   "nextflow -v"
ver nf-core    "nf-core --version"
ver java       "java -version"
ver singularity "singularity --version"
ver snakemake  "snakemake --version"
# R package versions actually loaded by the analysis scripts:
command -v Rscript >/dev/null 2>&1 && Rscript -e '
  for (p in c("readr","dplyr","tidyr","stringr","ggplot2","tibble"))
    cat(sprintf("%-14s %s\n", p, tryCatch(as.character(packageVersion(p)), error=function(e) "NOT INSTALLED")))' 2>/dev/null
