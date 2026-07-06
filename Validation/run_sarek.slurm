#!/bin/bash
###############################################################################
# run_sarek.slurm   (nf-core/sarek 3.9.0)  - ALTERNATIVE caller to rnavar
#
# BWA-MEM2 -> MarkDuplicates -> BQSR -> HaplotypeCaller (GERMLINE, per-sample), producing
# per-sample VCFs that feed the same ADMIXTURE ancestry stage as rnavar. Use this when the
# input is DNA (WGS/WES). For CCLE RNA-seq, rnavar (STAR + SplitNCigarReads) is the correct
# caller - Sarek has no RNA-aware split/hard-filter step, so treat this path as off-label for
# RNA and prefer run_rnavar.slurm unless you specifically have DNA reads.
#
# SAMPLE SHEET: build the sarek-format sheet with the SAME script as rnavar, adding
# `--samplesheet sarek`:
#   python Validation/map_srr_to_ancestry.py --ancestry Validation/ccle_ancestry_master.csv \
#       --runtable Validation/SraRunTable.csv --fastq-dir <FASTQ_DIR> \
#       --samplesheet sarek --out-prefix <TABLE_DIR>/ccle_anc
#   -> writes <prefix>_sarek_<POPS>.csv  (columns: patient,sex,status,sample,lane,fastq_1,fastq_2)
#
# DOWNSTREAM: set `caller: sarek` in config.yaml. Sarek writes per-sample VCFs under
#   <outdir>/variant_calling/haplotypecaller/<sample>/<sample>.haplotypecaller.vcf.gz
# (no ".filtered" in the name - the ancestry glob is caller-aware via config vcf_glob.sarek).
#
# BEFORE sbatch:  mkdir -p slurm_logs ; set the anchors (or export ANCESTRY_WORK / ANCESTRY_GATK) ;
#                 edit the SBATCH partition/account below.
###############################################################################
#SBATCH --job-name=sarek_driver
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
GATK_DIR="${ANCESTRY_GATK:-/path/to/gatk_bundle}"         # GATK known-sites bundle (bare-contig *_rename)
FASTA="$WORK/ref/FASTA/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
DBSNP="$GATK_DIR/dbsnp_146.hg38_rename.vcf.gz"
INDELS="$GATK_DIR/Mills_and_1000G_gold_standard_rename.indels.hg38.vcf.gz"
SHEET="${SAMPLE_SHEET:-$WORK/table/ccle_anc_sarek_AFR_EUR_AMR.csv}"
OUTDIR="$WORK/sarek"
WORKDIR="$WORK/sarek_work"
CONFIG="${SLURM_SUBMIT_DIR:-$PWD}/nfcore_hpc.config"      # production resource/infra config
# ================================================

PARAMS="${SLURM_SUBMIT_DIR:-$PWD}/sarek_params.yaml"

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
input:            "$SHEET"
outdir:           "$OUTDIR"
tools:            "haplotypecaller"
step:             "mapping"
wes:              false
aligner:          "bwa-mem2"
igenomes_ignore:  true
fasta:            "$FASTA"
fasta_fai:        "$FASTA.fai"
dbsnp:            "$DBSNP"
dbsnp_tbi:        "$DBSNP.tbi"
known_indels:     "$INDELS"
known_indels_tbi: "$INDELS.tbi"
joint_germline:   false
save_reference:   true
YAML
echo "[INFO] wrote params file: $PARAMS"; cat "$PARAMS"
nextflow run nf-core/sarek -r 3.9.0 \
    -profile singularity \
    -c "$CONFIG" \
    -work-dir "$WORKDIR" \
    -params-file "$PARAMS" \
    -resume
