#!/bin/bash
###############################################################################
# init_nfcore_singularity.sh - one-time bootstrap of the variant-calling launcher:
#   Nextflow  ->  nf-core  ->  pull rnavar & sarek  ->  pre-stage Singularity images.
#
# Run ONCE per cluster/checkout (or as a Snakemake setup rule) on an internet-connected
# login/submit node, before the caller step. It pulls the pinned pipeline revisions and
# stages their Singularity images into a shared cache so compute nodes never hit the network.
#
# Env: the nextflow conda env (envs/nextflow.yaml, name: nextflow), plus a container runtime
#      (Singularity or Apptainer) on PATH. On Lmod clusters these may come from `module load`
#      (e.g. UAB Cheaha: `module load Java/19.0.2 Singularity`); off Lmod, provide them yourself.
###############################################################################
set -euo pipefail

RNAVAR_VER="1.3.0"
SAREK_VER="3.9.0"
: "${SCRATCH:=/scratch/$USER}"
export NXF_SINGULARITY_CACHEDIR="${NXF_SINGULARITY_CACHEDIR:-$SCRATCH/.singularity_cache}"
export NXF_HOME="${NXF_HOME:-$SCRATCH/.nextflow}"
mkdir -p "$NXF_SINGULARITY_CACHEDIR" "$NXF_HOME"

# --- 1) toolchain present? ---------------------------------------------------
command -v nextflow >/dev/null || { echo "[ERROR] nextflow not on PATH - conda activate the nextflow env (envs/nextflow.yaml)"; exit 1; }
if   command -v singularity >/dev/null; then SING=singularity
elif command -v apptainer   >/dev/null; then SING=apptainer
else echo "[ERROR] singularity/apptainer not found - install one, or on Lmod clusters: module load Singularity (or Apptainer)"; exit 1; fi
echo "[INFO] nextflow $(nextflow -v 2>&1 | awk 'NR==1{print $NF}')  |  $SING $($SING --version 2>&1 | awk '{print $NF}')"
echo "[INFO] NXF_SINGULARITY_CACHEDIR=$NXF_SINGULARITY_CACHEDIR"

# --- 2) pull the pipelines at pinned revisions -------------------------------
nextflow pull nf-core/rnavar -r "$RNAVAR_VER"
nextflow pull nf-core/sarek  -r "$SAREK_VER"

# --- 3) pre-stage Singularity images into the cache (recommended) ------------
# nf-core download writes every container the pipeline uses into the cache, so runs are
# fully offline afterwards. Flags follow nf-core tools 2.14; adjust if your CLI differs.
if command -v nf-core >/dev/null; then
  for spec in "rnavar:$RNAVAR_VER" "sarek:$SAREK_VER"; do
    p="${spec%%:*}"; v="${spec##*:}"
    echo "[INFO] staging $p $v Singularity images -> $NXF_SINGULARITY_CACHEDIR"
    nf-core download "$p" -r "$v" \
      --container-system singularity \
      --container-cache-utilisation amend \
      --compress none \
      --outdir "$SCRATCH/nfcore_${p}_${v}" \
      || echo "[WARN] nf-core download for $p failed - images will pull at first run via -profile singularity"
  done
else
  echo "[WARN] nf-core CLI not found; skipping pre-stage. Images pull on first run (needs internet on the compute node)."
fi

echo "[DONE] nf-core rnavar ($RNAVAR_VER) + sarek ($SAREK_VER) ready."
echo "[NEXT] submit Validation/run_rnavar.slurm  (or run_sarek.slurm)."
