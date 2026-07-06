#!/bin/bash
###############################################################################
# run_rnavar.slurm   (nf-core/rnavar 1.3.0)
#
# STAR 2-pass -> MarkDuplicates -> SplitNCigarReads -> (BQSR) -> HaplotypeCaller
# -> RNA hard filters, for the AFR/EUR/AMR CCLE set. Output VCFs feed ancestry calling.
# This sbatch job is only the lightweight Nextflow DRIVER; STAR/GATK steps are
# dispatched to Slurm by nfcore_hpc.config.
#
#-----------------------------------------------------------------------------
# REFERENCE SET (Ensembl-style, no 'chr' prefix, self-consistent):
#   FASTA  Homo_sapiens.GRCh38.dna.primary_assembly.fa     (contigs 1,2,...,X,MT)
#   GTF    Homo_sapiens.GRCh38.112.gtf
#   VCFs   *_rename.vcf.gz (+ matching .tbi)               (renamed to no-chr; see refs/prep_refs.sh)
#
# STAR INDEX NOTE: rnavar 1.3.0 uses a recent STAR and will rebuild its own index; we do NOT pass
# --star_index and add --save_reference so the rebuilt index is written under --outdir. Reuse it on
# later runs with:  --star_index <outdir>/.../star  (skips the one-time rebuild).
#-----------------------------------------------------------------------------
# BEFORE sbatch:  mkdir -p slurm_logs ; set the anchors (or export ANCESTRY_WORK / ANCESTRY_GATK) ;
#                 edit the SBATCH partition/account below.  (-resume is on, so re-submit to continue.)
###############################################################################
#SBATCH --job-name=rnavar_driver
#SBATCH --partition=long              # EDIT: your long-running partition (Cheaha: long)
#SBATCH --time=150:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=slurm_logs/%x_%j.out
#SBATCH --error=slurm_logs/%x_%j.err
##SBATCH --account=<your_account>     # EDIT/uncomment if your cluster requires --account
##SBATCH --mail-user=you@example.org  # EDIT/uncomment for email notifications
##SBATCH --mail-type=ALL
set -euo pipefail
: "${SCRATCH:=/scratch/$USER}"

# ===================== EDIT =====================
# Configure via env anchors (same as the Snakemake config) or edit the fallbacks after ':-'.
WORK="${ANCESTRY_WORK:-/path/to/ancestry_work}"           # writable working area
REFBUNDLE="${ANCESTRY_GATK:-/path/to/gatk_bundle}"        # GATK known-sites bundle (bare-contig *_rename)
FASTA="$WORK/ref/FASTA/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
GTF="$WORK/ref/GTF/Homo_sapiens.GRCh38.112.gtf"
DBSNP="$REFBUNDLE/dbsnp_146.hg38_rename.vcf.gz"
INDELS="$REFBUNDLE/Mills_and_1000G_gold_standard_rename.indels.hg38.vcf.gz"
SHEET="${SAMPLE_SHEET:-$WORK/table/ccle_anc_rnavar_AFR_EUR_AMR.csv}"
OUTDIR="$WORK/rnavar"
WORKDIR="$WORK/rnavar_work"
READLEN=101
# ================================================

PARAMS="${SLURM_SUBMIT_DIR:-$PWD}/rnavar_params.yaml"

# HPC modules (Lmod). Skipped automatically off Lmod; provide Java + the nextflow env yourself otherwise.
if command -v module >/dev/null 2>&1; then
    module reset 2>/dev/null || true
    module load Java/19.0.2 2>/dev/null || true
    module load Anaconda3 2>/dev/null || true
fi
NFCORE_ENV="${NFCORE_ENV:-nextflow}"                      # conda env with nextflow + nf-core (envs/nextflow.yaml)
if command -v conda >/dev/null 2>&1; then
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "$NFCORE_ENV" 2>/dev/null || { echo "[ERROR] conda env '$NFCORE_ENV' not found (envs/nextflow.yaml)"; exit 1; }
fi
command -v nextflow >/dev/null 2>&1 || { echo "[ERROR] nextflow not on PATH"; exit 1; }

export NXF_SINGULARITY_CACHEDIR="${NXF_SINGULARITY_CACHEDIR:-$SCRATCH/.singularity_cache}"
mkdir -p "$NXF_SINGULARITY_CACHEDIR" "$OUTDIR" "$WORKDIR"
cat > "$PARAMS" <<YAML
input:  "$SHEET"
outdir: "$OUTDIR"
fasta:  "$FASTA"
gtf:    "$GTF"
igenomes_ignore: true
aligner: star
star_twopass: true
read_length: $READLEN
save_reference: true
dbsnp:            "$DBSNP"
dbsnp_tbi:        "$DBSNP.tbi"
known_indels:     "$INDELS"
known_indels_tbi: "$INDELS.tbi"
star_max_memory_bamsort: 64424509440   # STAR --limitBAMsortRAM ~60 GiB (>48 GiB STAR needed)
star_bins_bamsort: 100                  # STAR --outBAMsortingBinsN (more bins = less peak RAM)
YAML
echo "[INFO] wrote params file: $PARAMS"; cat "$PARAMS"
nextflow run nf-core/rnavar -r 1.3.0 \
    -profile singularity \
    -c "${SLURM_SUBMIT_DIR:-$PWD}/nfcore_hpc.config" \
    -work-dir "$WORKDIR" \
    -params-file "$PARAMS" \
    -resume
