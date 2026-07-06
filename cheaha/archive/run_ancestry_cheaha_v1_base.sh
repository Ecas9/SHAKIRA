#!/bin/bash
###############################################################################
# run_ancestry_cheaha.sh
#
# Standalone SLURM driver for the RNAvar -> ADMIXTURE ancestry pipeline, using
# Cheaha's PRE-EXISTING environment modules + conda envs (NOT the Snakemake
# conda-env files). Self-contained alternative to ../workflow/Snakefile; runs the
# whole pipeline as ONE batch job.
#
# REFERENCE BUILD: GRCh38 only. The study (rnavar) VCFs and the 1000G panel are
# both GRCh38 with bare contigs (1..22), so there is NO liftover and NO per-chr
# fan-out — the flow is linear. (An hg19/liftover variant lives in the Snakemake
# version via ref_build: hg19 if ever needed; it is intentionally not here.)
#
# WHAT IT DOES:
#   per-sample rnavar HaplotypeCaller VCFs (GRCh38, bare contigs)
#     -> ingest (PASS-only, autosomes 1..22, biallelic SNPs)
#     -> combine into one cohort VCF        [the "combine all filtered VCFs" step]
#     -> site-missingness QC
#     -> harmonize against the GRCh38 1000G panel (subset to study SNPs, biallelic)
#     -> plink beds -> set var-ids (+ reference-side MAF) -> bmerge -> intersect
#     -> variant QC on the merged panel (long-range-LD exclude; optional geno/maf/hwe)
#     -> LD-prune -> drop all-missing individuals (--mind) -> build .pop labels
#     -> supervised ADMIXTURE  +  cross-validation K sweep
#
# QC NOTE (why this differs from generic GWAS QC): for SUPERVISED ancestry the legacy
# pooled `--geno 0.1 --maf 0.05 --hwe 1e-6` trio is NOT re-applied to the merged panel.
# That trio came from a per-RACE-STRATUM REGARDS QC (within ONE ancestry group); on a
# pooled multi-ancestry panel it strips ancestry-informative markers (HWE/MAF) or is
# inert (geno, since ~2000 1KG samples dominate). Marker quality is instead controlled by
# MAX_MISSING (study call-rate, Stage 3), a PER-SUPERPOPULATION reference MAF (Stage 4a, the
# >=1-superpop MAF union ported from the K3 pipeline) and long-range-LD exclusion (Stage 7b).
# The legacy pooled knobs remain available but default OFF.
#
# PREREQUISITES on Cheaha (adjust names in the EDIT block if yours differ):
#   * modules:  BCFtools, PLINK/1.90-foss-2016a, Anaconda3
#   * conda envs you already have:  `plink2` (has plink2), `admix` (has admixture)
#   * python3 on PATH (for pop.py) — Anaconda3 base provides it
#   * cheaha/pop.py present (shipped alongside this script)
#   * the GRCh38 panel VCF must be bgzip + tabix-indexed (.tbi/.csi)
#
# SUBMIT:   sbatch cheaha/run_ancestry_cheaha.sh
#   (edit the SBATCH header for your account/partition, then the EDIT block)
#
# NOTE on #SBATCH --chdir: SLURM applies it (and resolves %x_%j.out/.err relative to it) BEFORE
# the script runs, so the target dir must EXIST at submit time. On the FIRST run, create it once:
#     mkdir -p /path/to/work/admixture_run
#   (the script also `mkdir -p`s OUTDIR internally, but that happens too late for --chdir.)
# NOTE on #SBATCH lines generally: SLURM parses them BEFORE the shell, so they CANNOT reference
# EDIT-block variables. Change them here or override at submit time, keeping --time within the
# partition cap (`short` caps at 12h; use --partition=medium/long for more):
#     sbatch --partition=medium --time=48:00:00 --cpus-per-task=24 cheaha/run_ancestry_cheaha.sh
###############################################################################
#SBATCH --chdir=/path/to/work/amr
#SBATCH --job-name=rnavar_admixture
#SBATCH --partition=short                # ADMIXTURE (esp. the CV sweep) dominates runtime
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=./logs/%x_%j.out
#SBATCH --error=./logs/%x_%j.err
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
# OPTIONAL extra per-sample exclusion list (FID<tab>sample; the sample column is used).
# DROPPED by default: the admixed reference samples are now removed by EXCLUDE_SUBPOPS
# (ACB/ASW/MXL, Stage 4a) and the per-superpop MAF, so the old <98%-ancestry list
# (1K_98_co_exclude.txt) is redundant. Set a path only to ALSO drop a custom sample list.
KG_EXCLUDE=""

# --- population label input (make_pop) ---------------------------------------
# 1KG sample->population panel (sample <TAB> pop; has header). pop.py's own default.
KG_PANEL="$HELPER_DIR/1K_pops.txt"   # clean 2,946-sample panel in cheaha/ (repo-root copy was a duplicated error)
# NOTE: the pop->super-population mapping is kept INSIDE pop.py (its built-in
# 28-population map, covering all 1KG codes incl. CHD), so NO external
# populations.tsv is required — we don't pass --pop-table. (pop.py --pop-table can
# still overlay an authoritative table later if you ever produce one.)

# --- output location ---------------------------------------------------------
OUTDIR="/path/to/work/amr"

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
STUDY_COVERAGE=0.80     # require >= this fraction of the STUDY samples genotyped per merged SNP
                        # (Stage 7b). Implemented as napkin-math plink --geno on the merged bed:
                        # assuming the 1KG reference is ~fully called, the equivalent per-variant
                        # missingness threshold is (1-STUDY_COVERAGE)*n_study/n_total — the bcftools
                        # F_MISSING napkin math from variants/K3/K3_fixed_freq_no_asia.sh. This is
                        # what actually controls study coverage AFTER the 1KG join; --mind is
                        # per-INDIVIDUAL and does not. Stage 3 pre-filters study coverage before the
                        # merge; this re-enforces it on the merged set. Set "" to skip.
MERGE_GENO=""           # manual override: a FIXED plink --geno on the merged panel. When set it
                        # REPLACES the STUDY_COVERAGE napkin math. Leave "" to use STUDY_COVERAGE.

MIND=0.99               # drop individuals missing > this fraction of the pruned SNPs.
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
AUTOSOMES="$(seq -s, 1 22)"                      # "1,2,...,22" — drops X/Y/MT/alt contigs
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
module load parallel/20240722-GCCcore-13.3.0
# enable `conda activate` inside a non-interactive batch shell:
source "$(conda info --base)/etc/profile.d/conda.sh"

mkdir -p "$OUTDIR"/{normalized,cohort,qc,study,ref,bed,merge,pruned,admixture/cv}
cd "$OUTDIR"

assert_indexed() { [ -f "$1.tbi" ] || [ -f "$1.csi" ] || { echo "[ERROR] no index for $1 (run: tabix -p vcf $1)"; exit 1; }; }

# ============================================================================
# Stage 0 — resolve the sample list -> $OUTDIR/cohort/samples.resolved
#            (each line: "<sample_name><TAB><vcf_path>")
# ============================================================================
RESOLVED="$OUTDIR/cohort/samples.resolved"
# discover mode: build the sample sheet from the rnavar tree with `find`. Sample
# name = each *.haplotypecaller.filtered.vcf.gz filename with that suffix stripped
# (the per-sample prefix). Written to SAMPLE_SHEET so it exists as a reproducible
# record; both modes then read the sheet.
if [ "$INPUT_MODE" = discover ]; then
    mkdir -p "$(dirname "$SAMPLE_SHEET")"
    { printf 'sample\tvcf\n'
      find "$RNAVAR_DIR" -type f -name '*.haplotypecaller.filtered.vcf.gz' | sort | while read -r p; do
          name="$(basename "$p")"; name="${name%.haplotypecaller.filtered.vcf.gz}"
          printf '%s\t%s\n' "$name" "$p"
      done
    } > "$SAMPLE_SHEET"
    echo "[INFO] discovered $(($(wc -l < "$SAMPLE_SHEET") - 1)) sample(s) -> wrote sheet: $SAMPLE_SHEET"
fi
[ -s "$SAMPLE_SHEET" ] || { echo "[ERROR] sample sheet missing/empty: $SAMPLE_SHEET"; exit 1; }
# read the sheet (skip header, drop blank / '#'-comment lines) -> "name<TAB>path" list
tail -n +2 "$SAMPLE_SHEET" | awk 'NF && $1!~/^#/ {print $1"\t"$2}' > "$RESOLVED"
NSAMP=$(wc -l < "$RESOLVED")
[ "$NSAMP" -ge 1 ] || { echo "[ERROR] no samples resolved (mode=$INPUT_MODE)"; exit 1; }
echo "[INFO] $NSAMP sample(s) to process:"; cut -f1 "$RESOLVED" | paste -sd' ' -

# ============================================================================
# Stage 1 — ingest each rnavar VCF: drop non-PASS (RNA hard filters
#   FS/QD/SnpCluster/LowQual), keep autosomes, split multiallelics, keep
#   strictly biallelic SNPs. bgzip + tabix each. -> normalized/<sample>.norm.snps.vcf.gz
#   Per-sample work is independent, so it is run through GNU `parallel` when available
#   (the OLD pipeline's per-sample fan-out idiom; `module load parallel`), falling back
#   to a serial loop otherwise. ingest_one echoes its output path on stdout, which is
#   collected into cohort/merge.list (avoids concurrent appends to one file).
# ============================================================================
export AUTOSOMES OUTDIR                              # consumed by cheaha/ingest_one.sh
if command -v parallel >/dev/null 2>&1 && [ "$NSAMP" -gt 1 ]; then
    echo "[INFO] Stage 1: parallel ingest of $NSAMP sample(s) (GNU parallel -j $NPROC)"
    parallel --colsep '\t' -j "$NPROC" --halt soon,fail=1 \
        bash "$HELPER_DIR/ingest_one.sh" {1} {2} :::: "$RESOLVED" > "$OUTDIR/cohort/merge.list"
else
    echo "[INFO] Stage 1: serial ingest of $NSAMP sample(s) (GNU parallel not found)"
    : > "$OUTDIR/cohort/merge.list"
    while IFS=$'\t' read -r name path; do
        echo "[ingest] $name" >&2
        bash "$HELPER_DIR/ingest_one.sh" "$name" "$path" >> "$OUTDIR/cohort/merge.list"
    done < "$RESOLVED"
fi
echo "[COUNT] Stage 1: $(wc -l < "$OUTDIR/cohort/merge.list") per-sample VCF(s) ingested"

# ============================================================================
# Stage 2 — combine all per-sample VCFs into one cohort VCF.
#   `-m none` keeps split SNP alleles as separate biallelic records (so two
#   samples' different ALTs at one site are NOT re-collapsed into a multiallelic
#   that QC would then drop); norm + biallelic re-select keeps it tidy.
#   A 1-sample cohort is just copied through (bcftools merge needs >=2 inputs).
# ============================================================================
COHORT="$OUTDIR/cohort/merge.vcf.gz"
if [ "$(wc -l < "$OUTDIR/cohort/merge.list")" -le 1 ]; then
    bcftools view --threads "$NPROC" -Oz -o "$COHORT" "$(cat "$OUTDIR/cohort/merge.list")"
else
    bcftools merge -m none -l "$OUTDIR/cohort/merge.list" --threads "$NPROC" -Ou \
      | bcftools norm -m -any --threads "$NPROC" -Ou \
      | bcftools view -m2 -M2 -v snps --threads "$NPROC" -Oz -o "$COHORT"
fi
bcftools index -t --threads "$NPROC" "$COHORT"
echo "[COUNT] Stage 2 cohort: $(bcftools index -n "$COHORT") variants, $(bcftools query -l "$COHORT" | wc -l) samples"

# ============================================================================
# Stage 3 — site-missingness QC: keep sites genotyped in >= MAX_MISSING of
#   samples (F_MISSING <= 1-MAX_MISSING, inclusive like vcftools --max-missing),
#   strictly biallelic SNPs, drop all-missing sites (-U). This QC'd cohort IS the
#   GRCh38 study set (no liftover). Emit its SNP positions for the 1KG subset.
# ============================================================================
STUDY="$OUTDIR/study/study_all_chroms.vcf.gz"
bcftools view -i "F_MISSING <= ${MISS_THR}" -m2 -M2 -v snps -U \
    --threads "$NPROC" -Oz -o "$STUDY" "$COHORT"
bcftools index -t --threads "$NPROC" "$STUDY"

SNPS="$OUTDIR/study/sodh_snps_final.txt"
bcftools query -f '%CHROM\t%POS\n' "$STUDY" > "$SNPS"
echo "[INFO] study SNPs: $(wc -l < "$SNPS")"

# ======================================