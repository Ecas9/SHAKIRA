#!/bin/bash
###############################################################################
# run_analysis_cheaha.sh
#
# Standalone SLURM driver for the post-ADMIXTURE ancestry analysis: runs
# analyze_ancestry.R against a finished run_ancestry_cheaha.sh output to produce
# the ancestry plots + concordance vs the CCLE master reference.
#
# This is the analysis counterpart to run_ancestry_cheaha.sh — submit it after the
# ADMIXTURE job finishes. It uses a Cheaha R module (no Snakemake).
#
# WHAT IT PRODUCES (under OUT_DIR):
#   ancestry_study_proportions.csv   per cell line: AFR/EUR/EAS/SAS/AMR + dominant
#   ancestry_concordance.csv         per matched cell line: ours vs reference
#   concordance_summary.txt          dominant-super-pop agreement + per-pop r / RMSE
#   plot_ancestry_stacked.pdf        study ancestry stacked bars
#   plot_concordance_scatter.pdf     ours vs reference %, faceted by super-pop
#   plot_dominant_confusion.pdf      dominant super-pop confusion matrix
#
# PREREQUISITES on Cheaha:
#   * an R module whose library has the tidyverse packages
#     (readr, dplyr, tidyr, stringr, ggplot2, tibble) — the script checks and stops
#     with a clear message if any are missing (then install.packages(...) once, or
#     load an R module that bundles tidyverse).
#   * cheaha/analyze_ancestry.R present (shipped alongside this script).
#   * a finished ADMIXTURE run (RUN_DIR/admixture/prunedData.<K>.Q etc.) and the CCLE
#     master CSV.
#
# SUBMIT:   sbatch cheaha/run_analysis_cheaha.sh
###############################################################################
#SBATCH --chdir=/path/to/work/admixture_run/
#SBATCH --job-name=ancestry_analysis
#SBATCH --partition=express
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=./logs/%x_%j.out
#SBATCH --error=./logs/%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=youremail@example.org

set -euo pipefail

# =============================================================================
# ============================  EDIT BLOCK  ===================================
# =============================================================================
REPO_DIR="/path/to/shakira"
HELPER_DIR="$REPO_DIR/cheaha"                    # contains analyze_ancestry.R

# Finished ADMIXTURE run dir (the OUTDIR from run_ancestry_cheaha.sh). The R script
# reads RUN_DIR/admixture/prunedData.<K>.Q, RUN_DIR/pruned/prunedData.{fam,pop}.
#RUN_DIR="/path/to/work/admixture_run"
RUN_DIR="/path/to/work/amr"

# CCLE ancestry master reference CSV (cLine_ID, AFR/EUR/EAS/SAS/AMR, dominant_superpop,
# match_keys, ...). Copy it to the cluster and point here.
REF_CSV="/path/to/shakira/Validation/ccle_ancestry_master.csv"

# Where to write the plots + tables.
OUT_DIR="$RUN_DIR/analysis"

# Supervised K used in the ADMIXTURE run (must match run_ancestry_cheaha.sh K_SUPERVISED).
K_SUPERVISED=3

# Cheaha R module (set the exact version; bare "R" loads the default). Must provide tidyverse.
MOD_R="R/4.2.0-foss-2022b"
# =============================================================================
# ==========================  END EDIT BLOCK  =================================
# =============================================================================

module reset 2>/dev/null || module purge 2>/dev/null || true
module load "$MOD_R"

# fail early with a clear message if the tidyverse packages aren't in this R's library
Rscript -e 'p <- c("readr","dplyr","tidyr","stringr","ggplot2","tibble");
            m <- p[!sapply(p, requireNamespace, quietly = TRUE)];
            if (length(m)) { cat("[ERROR] missing R packages:", paste(m, collapse=", "),
              "\n  -> load an R module that bundles tidyverse, or run install.packages() once.\n");
              quit(status = 1) }'

mkdir -p "$OUT_DIR"

# CV-error-vs-K plot first: it only needs the CV logs, so it still runs for a K=3 panel run
# even if the concordance analysis below (which expects the CCLE master + all 5 super-pops)
# does not apply / errors.
Rscript "$HELPER_DIR/plot_cv_error.R" "$RUN_DIR" "$OUT_DIR"

Rscript "$HELPER_DIR/analyze_ancestry.R" "$RUN_DIR" "$REF_CSV" "$OUT_DIR" "$K_SUPERVISED"

echo "[DONE] analysis written to: $OUT_DIR"
echo "[DONE] CV error curve:      $OUT_DIR/plot_cv_error.pdf  ($OUT_DIR/cv_error.csv)"
echo "[DONE] headline metrics:    $OUT_DIR/concordance_summary.txt"
