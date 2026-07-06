#!/bin/bash
#SBATCH --chdir=/path/to/work/admixture_run/
#SBATCH --job-name=rnavar_admixture
#SBATCH --partition=express              # ADMIXTURE (esp. the CV sweep) dominates runtime
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=./%x_%j.out
#SBATCH --error=./%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=youremail@example.org

set -euo pipefail

# =============================================================================
# ============================  EDIT BLOCK  ===================================
# =============================================================================
# --- repo / helper locations -------------------------------------------------
# Set REPO_DIR to your checkout. We DON'T auto-derive it from the script path
# because sbatch copies the script into a private spool dir, so $BASH_SOURCE is
# unreliable. HELPER_DIR must contain pop.py (shipped here).
REPO_DIR="/path/to/shakira"
HELPER_DIR="$REPO_DIR/cheaha"

# --- INPUT: choose ONE mode --------------------------------------------------
# "discover" = find EVERY *.haplotypecaller.filtered.vcf.gz under RNAVAR_DIR and
#              WRITE a 2-col sample sheet to SAMPLE_SHEET (sample = the filename
#              prefix, i.e. name with .haplotypecaller.filtered.vcf.gz stripped).
# "sheet"    = READ an existing 2-column TSV at SAMPLE_SHEET (sample <ws> vcf_path).
INPUT_MODE="discover"
RNAVAR_DIR="/path/to/work/rnavar/variant_calling/"
SAMPLE_SHEET="/path/to/work/table/samples.tsv"   # written (discover) | read (sheet)

# --- GRCh38 1000 Genomes reference panel -------------------------------------
# A SINGLE combined panel VCF, all autosomes, bare contigs (1..22). The one you use.
# MUST be bgzipped + tabix-indexed (.tbi/.csi).
KG_GRCH38="/path/to/data/QC/1kgenome.vcf.gz"
# HIGH-QUALITY REFERENCE: list of admixed 1KG samples to EXCLUDE from the reference
# (FID<tab>sample; the sample column is used), produced by the K3 high-quality pipeline
# (run_1k_highqual_cheaha.sh -> 1k_analyze_samples.R). Set "" to use the full panel.
KG_EXCLUDE="/path/to/data/K3_15_mill_filter/1K_98_co_exclude.txt"

# --- population label input (make_pop) ---------------------------------------
# 1KG sample->population panel (sample <TAB> pop; has header). pop.py's own default.
KG_PANEL="$REPO_DIR/1K_pops.txt"
# NOTE: the pop->super-population mapping is kept INSIDE pop.py (its built-in
# 28-population map, covering all 1KG codes incl. CHD), so NO external
# populations.tsv is required — we don't pass --pop-table. (pop.py --pop-table can
# still overlay an authoritative table later if you ever produce one.)

# --- output location ---------------------------------------------------------
OUTDIR="/path/to/work/admixture_run"

# --- analysis parameters -----------------------------------------------------
K_SUPERVISED=3          # supervised ADMIXTURE K = number of reference super-populations. MUST equal
                        # the count in SUPERPOPS_OF_INTEREST (3 = AFR/EUR/AMR for this run). Use 5 only
                        # if you also restore EAS/SAS to SUPERPOPS_OF_INTEREST.
K_MIN=1                 # cross-validation sweep lower bound
K_MAX=5                 # cross-validation sweep upper bound
MAX_MISSING=0.8         # keep sites genotyped in >= this fraction of samples (80%). For RNA-seq
                        # across heterogeneous cell lines a strict 0.9 keeps only ~6k ubiquitously
                        # expressed SNPs; lower (0.8/0.5/0.2) to retain more ancestry-informative
                        # sites and fewer all-missing cells. Raise toward 0.9 for denser cohorts.
LD_WINDOW=50            # plink --indep-pairwise window
LD_STEP=10              # plink --indep-pairwise step
LD_R2=0.1               # plink --indep-pairwise r^2 threshold

# --- variant-level QC, tuned for SUPERVISED ANCESTRY (see Stage 4a / Stage 7b) ---------
# Defaults intentionally differ from generic GWAS QC; each filter is skipped when empty.
# PER-SUPERPOPULATION reference MAF (Stage 4a) — keep a 1KG position if it is COMMON
# (MAF > REF_MAF) in AT LEAST ONE superpopulation of interest (a UNION). Ported from the
# validated K3 pipeline (variants/K3/MAF.sh + 1kGenome_filter.sh) via
# cheaha/ref_maf_per_superpop.sh. This is NOT a single pooled --maf: pooled MAF is dominated
# by the largest super-pop and discards population-private ancestry-informative markers.
SUPERPOPS_OF_INTEREST="AFR EUR AMR"   # reference super-pops of interest (this run = K=3). Whole
                                      # non-of-interest super-pops (EAS+SAS, incl. CHD) are DROPPED
                                      # from the reference. Add "EAS SAS" for K=5 (+ K_SUPERVISED=5).
EXCLUDE_SUBPOPS="ACB ASW MXL"         # admixed SUB-pops to HARD-DROP even though their super-pop is of
                                      # interest: ACB/ASW = admixed AFR, MXL = admixed AMR (documented
                                      # in variants/K3/1k_filter.py; manuscript ref 31, Rudy 2012).
                                      # Set "" to keep them. Applied with SUPERPOPS_OF_INTEREST to the
                                      # reference SAMPLE set in Stage 4a.
REF_MAF=0.05                          # keep MAF > this within >=1 super-pop of interest (per-super-pop
                                      # UNION). Set "" to skip the MAF positions (sample keep still applies).
LRLD_REGIONS="$HELPER_DIR/high-LD-regions-hg38-GRCh38.txt"   # GRCh38 long-range high-LD regions
                        # (plinkQC, meyer-lab-cshl) EXCLUDED before LD pruning (Stage 7b) — pairwise
                        # --indep-pairwise does not fully break these megabase blocks (MHC, chr8
                        # inversion, ...) and they distort ADMIXTURE. The shipped file is chr-prefixed;
                        # Stage 7b strips "chr" to match our bare contigs (1..22). Set "" to skip.
MERGE_MAF=""            # pooled MAF on study+1KG. OFF: a globally-rare SNP can be a strong AIM
                        # (common in one super-pop, ~0 elsewhere) — pooled MAF discards it, and it
                        # is dominated by 1KG counts anyway. Use the per-superpop REF_MAF (Stage 4a).
MERGE_HWE=""            # pooled Hardy-Weinberg p-threshold. OFF: on a multi-ancestry panel the
                        # Wahlund effect makes AIMs fail HWE, so pooled --hwe strips signal. Legacy
                        # hwe 1e-6 was per-RACE-STRATUM (within one ancestry). (Set 1e-6 for legacy.)
MERGE_GENO=""           # per-variant call-rate on the merged panel. OFF: near-inert because the
                        # ~2000 fully-called 1KG samples dominate the denominator; study call-rate
                        # is controlled by MAX_MISSING (Stage 3). A loose 0.1 only guards merge dropout.

MIND=0.2               # drop individuals missing > this fraction of the pruned SNPs.
                        # ADMIXTURE aborts on any all-missing individual; RNA-seq cell lines
                        # with ~no coverage at the common-SNP set (incl. any whose VCF had 0
                        # PASS SNPs) are removed and logged to prunedData.irem.

# --- module / conda env names (edit if your cluster names differ) ------------
MOD_BCFTOOLS="BCFtools"
MOD_PLINK1="PLINK/1.90-foss-2016a"
MOD_ANACONDA="Anaconda3"
CONDA_ENV_PLINK2="plink2"             # your env that provides plink2
CONDA_ENV_ADMIX="admix"               # your env that provides admixture
# =============================================================================
# ==========================  END EDIT BLOCK  =================================
# =============================================================================

NPROC="${SLURM_CPUS_PER_TASK:-8}"
MISS_THR=$(awk -v m="$MAX_MISSING" 'BEGIN{printf "%.6f", 1-m}')   # F_MISSING <= 1-MAX_MISSING

echo "[INFO] build=GRCh38  input_mode=$INPUT_MODE  nproc=$NPROC  F_MISSING<=$MISS_THR"
echo "[INFO] outdir=$OUTDIR"

# ----------------------------------------------------------------------------
# Environment: load the module tool-stack once. The two conda-only tools
# (plink2, admixture) are activated per-stage further down — exactly the legacy
# module/conda split. `conda activate` needs the conda shell hook sourced first.
# ----------------------------------------------------------------------------
module reset 2>/dev/null || module purge 2>/dev/null || true
module load "$MOD_BCFTOOLS"
module load "$MOD_PLINK1"
module load "$MOD_ANACONDA"
# enable `conda activate` inside a non-interactive batch shell:
source "$(conda info --base)/etc/profile.d/conda.sh"



plink --bfile "$OUTDIR/pruned/prunedData_all" --mind "$MIND" \
      --keep-allele-order --allow-no-sex --make-bed --out "$OUTDIR/pruned/prunedData" --threads "$NPROC"
if [ -s "$OUTDIR/pruned/prunedData.irem" ]; then
    echo "[WARN] dropped $(wc -l < "$OUTDIR/pruned/prunedData.irem") individual(s) (> ${MIND} missing at the pruned set):"
    awk '{print $2}' "$OUTDIR/pruned/prunedData.irem" | paste -sd' ' -
fi

# ============================================================================
# Stage 8 — build the supervised .pop file: 1 label per .fam row IN .fam ORDER;
#   1KG individuals get their super-population, study (cell-line) individuals are
#   left blank so ADMIXTURE estimates them against the fixed reference.
# ============================================================================
python3 "$HELPER_DIR/pop.py" --fam "$OUTDIR/pruned/prunedData.fam" \
        --sample-pop "$KG_PANEL" --out "$OUTDIR/pruned/prunedData.pop"

# guard: .pop must be 1:1 with .fam AND have at least one labeled reference individual
nfam=$(wc -l < "$OUTDIR/pruned/prunedData.fam")
npop=$(wc -l < "$OUTDIR/pruned/prunedData.pop")
[ "$nfam" -eq "$npop" ] || { echo "[ERROR] .pop ($npop) != .fam ($nfam) — not 1:1"; exit 1; }
[ "$(grep -c . "$OUTDIR/pruned/prunedData.pop" || true)" -gt 0 ] || { echo "[ERROR] .pop is all blank"; exit 1; }

# ============================================================================
# Stage 9 — supervised ADMIXTURE at K, plus the cross-validation K sweep on the
#   SAME pruned bed. ADMIXTURE writes <prefix>.K.Q/.P to the CWD and reads
#   <prefix>.pop next to the bed, so we cd into the output dirs.
# ============================================================================
conda activate "$CONDA_ENV_ADMIX"

# supervised estimate (reads ../pruned/prunedData.pop)
( cd "$OUTDIR/admixture" && \
  admixture --supervised "../pruned/prunedData.bed" "$K_SUPERVISED" -j"${NPROC}" )

# cross-validation sweep -> cv/log.<K>.out  (unsupervised CV, as in the legacy script)
for K in $(seq "$K_MIN" "$K_MAX"); do
    ( cd "$OUTDIR/admixture/cv" && \
      admixture --cv "../../pruned/prunedData.bed" "$K" -j"${NPROC}" | tee "log.${K}.out" )
done
conda deactivate

echo "[DONE] supervised Q: $OUTDIR/admixture/prunedData.${K_SUPERVISED}.Q"
echo "[DONE] pick K from CV error:  grep -h CV $OUTDIR/admixture/cv/log.*.out"