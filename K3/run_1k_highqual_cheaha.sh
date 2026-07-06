#!/bin/bash
###############################################################################
# run_1k_highqual_cheaha.sh
#
# One-job SLURM driver for the high-quality 1KG reference flow, parameterized by K.
# Given a PRE-BUILT merged panel bed (from a pipeline builder such as K3_pipeline.sh /
# 15_mill_pipeline.sh — NOT rebuilt here), it:
#   1. optionally drops a sample list (e.g. Asian pops for K=3)            -> <panel>
#   2. builds the supervised .pop (pop_jhu.py)
#   3. supervised ADMIXTURE at K                                          -> <panel>.K.Q
#   4. flags 1KG samples admixed at K (max ancestry < THR; 1k_analyze_samples.R)
#                                                                         -> 1K_98_co_exclude.txt
#   5. plink --remove those samples                                      -> Final_no_admixed_1K_samples
#   6. rebuild .pop + supervised ADMIXTURE at K on the cleaned panel
#
# K=3 is the default; pass another value (e.g. 5) as the first arg or set K below.
# NOTE for K!=3: K is the ADMIXTURE component count; the REFERENCE COMPOSITION (which
# superpops are in the panel) is set by the upstream builder + REMOVE_LIST. For 5
# superpops, build an all-pops panel and set REMOVE_LIST="".
#
# SUBMIT:   sbatch run_1k_highqual_cheaha.sh            # K=3
#           sbatch run_1k_highqual_cheaha.sh 5          # K=5
###############################################################################
#SBATCH --chdir=/path/to/work/admixture_run/K3_var
#SBATCH --job-name=1k_highqual_K_fix
#SBATCH --partition=short
#SBATCH --output=./logs/%x_%j.log
#SBATCH --error=./logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=04:00:00
#SBATCH --mem-per-cpu=3G
#SBATCH --mail-user=youremail@example.org
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# =============================================================================
# ============================  EDIT BLOCK  ===================================
# =============================================================================
REPO_DIR="/path/to/shakira"          # checkout with K3/1k_analyze_samples.R + pop_jhu.py
WORKDIR="/path/to/work/admixture_run/K3_var"   # holds the panel bed + outputs

K="${1:-3}"                                       # ADMIXTURE K (default 3; arg overrides)
START_PREFIX="prunedData"         # pre-built merged panel bed (from a builder)
REMOVE_LIST="/path/to/data/K3_10_mill_filter_fixed_freq/noAsia.txt"                          # samples to drop first (Asian pops for K=3); "" to skip
THR="0.98"                                        # "pure ancestry" threshold for the exclude list

POP_TABLE="/path/to/shakira/resources/1K_pops.txt"   # sample->pop table for pop_jhu.py

MOD_ANACONDA="Anaconda3/2022.05"
MOD_R="R"                                         # R module with tidyverse (for 1k_analyze_samples.R)
CONDA_ENV_PLINK="plink2"                          # provides plink
CONDA_ENV_POP="pandoc"                            # provides python3 + pandas + tabulate (pop_jhu.py)
CONDA_ENV_ADMIX="admix2"                          # provides admixture
# =============================================================================
# ==========================  END EDIT BLOCK  =================================
# =============================================================================

NPROC="${SLURM_CPUS_PER_TASK:-8}"
PANEL="${START_PREFIX}_panel"

cd "$WORKDIR"
module reset 2>/dev/null || module purge 2>/dev/null || true
module load "$MOD_ANACONDA"
source "$(conda info --base)/etc/profile.d/conda.sh"
echo "[INFO] K=$K  panel=$START_PREFIX  remove=[$REMOVE_LIST]  thr=$THR  workdir=$WORKDIR"

# --- 1) optional sample removal -> the panel we analyze ----------------------
conda activate "$CONDA_ENV_PLINK"
if [ -n "$REMOVE_LIST" ]; then
    #plink --bfile "$(/path/to/work/admixture_run/pruned/$START_PREFIX)" --remove "$REMOVE_LIST" --make-bed --out "$PANEL"
    cd /path/to/work/admixture_run/pruned/
    plink --bfile $START_PREFIX --remove "$REMOVE_LIST" --make-bed --out "$PANEL"
    cd $WORKDIR
    cp /path/to/work/admixture_run/pruned/prunedData_panel.bed .
    cp /path/to/work/admixture_run/pruned/prunedData_panel.fam .
    cp /path/to/work/admixture_run/pruned/prunedData_panel.bim .
else
    for e in bed bim fam; do cp "${START_PREFIX}.${e}" "${PANEL}.${e}"; done
fi
conda deactivate

# --- 2) supervised .pop for the panel ----------------------------------------
conda activate "$CONDA_ENV_POP"
#python3 "/path/to/shakira/pop_jhu.py" "${PANEL}.fam" "${PANEL}.pop" "$POP_TABLE"

python3 /path/to/shakira/K3/pop_jhu.py --fam "${PANEL}.fam" --out "${PANEL}.pop" --sample-pop "$POP_TABLE"
conda deactivate

# --- 3) supervised ADMIXTURE at K -> ${PANEL}.${K}.Q -------------------------
conda activate "$CONDA_ENV_ADMIX"
admixture "${PANEL}.bed" --supervised -j"$NPROC" "$K"
conda deactivate

# --- 4) flag admixed 1KG samples -> 1K_98_co_exclude.txt ---------------------
module load "$MOD_R"
Rscript "$REPO_DIR/K3/1k_analyze_samples.R" "$K" "$PANEL" "$WORKDIR" "$WORKDIR" "$THR"
echo "[INFO] excluding $(wc -l < 1K_98_co_exclude.txt) admixed sample(s)"

# --- 5) remove the admixed samples -> high-quality reference -----------------
conda activate "$CONDA_ENV_PLINK"
plink --bfile "$PANEL" --remove 1K_98_co_exclude.txt --make-bed --out Final_no_admixed_1K_samples
conda deactivate

# --- 6) rebuild .pop + supervised ADMIXTURE at K on the cleaned panel --------
conda activate "$CONDA_ENV_POP"
#python3 "$/path/to/shakira/pop_jhu.py" Final_no_admixed_1K_samples.fam Final_no_admixed_1K_samples.pop "$POP_TABLE"
python3 /path/to/shakira/K3/pop_jhu.py --fam Final_no_admixed_1K_samples.fam --out Final_no_admixed_1K_samples.pop --sample-pop "$POP_TABLE"
conda deactivate
conda activate "$CONDA_ENV_ADMIX"
admixture Final_no_admixed_1K_samples.bed --supervised -j"$NPROC" "$K"
conda deactivate

echo "[DONE] high-quality reference: $WORKDIR/Final_no_admixed_1K_samples.{bed,bim,fam}"
echo "[DONE] exclude list:           $WORKDIR/1K_98_co_exclude.txt"
echo "[DONE] ancestry Q:             $WORKDIR/Final_no_admixed_1K_samples.${K}.Q"
