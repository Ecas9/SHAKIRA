#!/bin/bash
###############################################################################
# run_snakemake.sh - launch the ancestry pipeline.
#
# Snakemake is the lightweight controller (like the former nextflow driver): it stays on the
# submit/login node and submits one Slurm job per rule via profiles/slurm. The upstream variant
# calling (rnavar/sarek) runs separately via Validation/run_{rnavar,sarek}.slurm.
#
#   bash run_snakemake.sh                 # full run
#   bash run_snakemake.sh -n              # dry run (DAG preview)
#   bash run_snakemake.sh --unlock        # release a stale lock
#
# Portability knobs (env vars):
#   CONDA_ENV     conda env with snakemake + the slurm executor plugin  (default: shakira)
#   SMK_PROFILE   Snakemake profile directory                           (default: profiles/slurm)
#                 -> set SMK_PROFILE=profiles/local to run on a single machine with no scheduler.
# The `module` calls below are Lmod-specific (UAB Cheaha and similar); they are skipped
# automatically on clusters without environment modules.
###############################################################################
set -euo pipefail
cd "$(dirname "$0")"                        # repo root

CONDA_ENV="${CONDA_ENV:-shakira}"
SMK_PROFILE="${SMK_PROFILE:-profiles/slurm}"

# Optional HPC environment modules (Lmod). Harmless no-ops off Lmod clusters.
if command -v module >/dev/null 2>&1; then
    module reset 2>/dev/null || true
    module load Anaconda3 2>/dev/null || true
fi

# Activate the controller conda env (snakemake + slurm executor plugin).
if ! command -v conda >/dev/null 2>&1; then
    echo "[ERROR] conda not on PATH. Install Miniforge/Miniconda, then:" >&2
    echo "        conda env create -f environment.yml" >&2
    exit 1
fi
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV" 2>/dev/null || {
    echo "[ERROR] conda env '$CONDA_ENV' not found. Create it with:" >&2
    echo "        conda env create -f environment.yml     # env name: shakira" >&2
    echo "        (or point CONDA_ENV=<your_env> at an env that has snakemake + the slurm plugin)" >&2
    exit 1
}

exec snakemake -s workflow/Snakefile --use-conda --profile "$SMK_PROFILE" "$@"
