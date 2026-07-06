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

# ============================================================================
# Stage 4a — build the reference SAMPLE keep-set + per-superpop MAF positions, via
#   cheaha/ref_maf_per_superpop.sh. Keep-set = SUPERPOPS_OF_INTEREST super-pops, MINUS the
#   admixed EXCLUDE_SUBPOPS (ACB/ASW/MXL), MINUS KG_EXCLUDE's <98% samples — so EAS/SAS and
#   the admixed sub-pops are dropped from the reference. When REF_MAF is set it also emits the
#   >=1-superpop MAF UNION positions (per-superpop two-step, not pooled). Ported from the
#   validated K3 pipeline (variants/K3/MAF.sh + 1kGenome_filter.sh + 1k_filter.py).
#   Outputs: ref/ref_keep_samples.txt and ref/maf_pass.positions (CHROM<TAB>POS).
# ============================================================================
if [ -n "$KG_EXCLUDE" ] && [ ! -s "$KG_EXCLUDE" ]; then
    echo "[ERROR] KG_EXCLUDE set but missing/empty: $KG_EXCLUDE — generate it via the K3 high-quality pipeline, or set KG_EXCLUDE=\"\""; exit 1
fi
REF_KEEP="$OUTDIR/ref/ref_keep_samples.txt"
MAF_PASS="$OUTDIR/ref/maf_pass.positions"
SUPERPOPS_OF_INTEREST="$SUPERPOPS_OF_INTEREST" EXCLUDE_SUBPOPS="$EXCLUDE_SUBPOPS" \
  REF_MAF="$REF_MAF" KG_EXCLUDE="$KG_EXCLUDE" NPROC="$NPROC" \
  bash "$HELPER_DIR/ref_maf_per_superpop.sh" "$KG_GRCH38" "$KG_PANEL" "$REF_KEEP" "$MAF_PASS" "$OUTDIR/ref/maf_work"
[ -s "$MAF_PASS" ] || MAF_PASS=""   # empty when REF_MAF="" (no MAF positions step)

# ============================================================================
# Stage 4 — harmonize the GRCh38 1KG panel: subset to the study SNP positions (further
#   intersected with the Stage 4a MAF-passing positions when enabled), split multiallelics,
#   keep biallelic SNPs (symmetric with ingest, so a 1KG multiallelic carrying a biallelic
#   SNP allele still matches a study var-id). -R is index-backed, so the panel MUST carry
#   a .tbi/.csi.
# ============================================================================
REF="$OUTDIR/ref/sodh_in_1kg_all_chroms.vcf.gz"
assert_indexed "$KG_GRCH38"

# Reference positions = study SNPs, intersected with the Stage 4a per-superpop MAF set when on.
# (Gating only the reference is enough: a study SNP absent from the MAF-passed reference is
# dropped at the Stage 7 var-id intersection, so it never reaches ADMIXTURE.)
SNPS_USE="$SNPS"
if [ -n "$MAF_PASS" ]; then
    SNPS_USE="$OUTDIR/ref/study_snps_maf_pass.txt"
    awk 'NR==FNR{a[$1"\t"$2]=1; next} (($1"\t"$2) in a)' "$MAF_PASS" "$SNPS" > "$SNPS_USE"
    echo "[INFO] study SNPs surviving per-superpop reference MAF: $(wc -l < "$SNPS_USE") of $(wc -l < "$SNPS")"
    [ -s "$SNPS_USE" ] || { echo "[ERROR] no study SNPs passed the per-superpop reference MAF (Stage 4a)"; exit 1; }
fi

# Restrict the reference panel to the Stage 4a keep-samples (of-interest super-pops, minus the
# admixed EXCLUDE_SUBPOPS, minus the <98% KG_EXCLUDE list). --force-samples tolerates any keep
# sample not present in the panel.
KG_SFILTER="-S $REF_KEEP --force-samples"
echo "[INFO] reference restricted to $(wc -l < "$REF_KEEP") keep-sample(s)"

bcftools view -R "$SNPS_USE" $KG_SFILTER --threads "$NPROC" -Ou "$KG_GRCH38" \
  | bcftools norm -m -any --threads "$NPROC" -Ou \
  | bcftools view -m2 -M2 -v snps --threads "$NPROC" -Oz -o "$REF"
bcftools index -t --threads "$NPROC" "$REF"
echo "[COUNT] Stage 4 reference: $(bcftools index -n "$REF") variants, $(bcftools query -l "$REF" | wc -l) samples"

# ============================================================================
# Stage 5 — PLINK beds. --double-id makes FID==IID==sample so make_pop can join
#   study IIDs against the 1KG panel (and find them absent -> blank label).
#   --keep-allele-order preserves REF/ALT; --vcf-half-call m treats any '*'-split
#   half-calls as missing instead of letting plink abort.
# ============================================================================
plink --vcf "$STUDY" --double-id --vcf-half-call m --keep-allele-order --allow-no-sex \
      --make-bed --out "$OUTDIR/bed/study" --threads "$NPROC"
plink --vcf "$REF"   --double-id --vcf-half-call m --keep-allele-order --allow-no-sex \
      --make-bed --out "$OUTDIR/bed/ref"   --threads "$NPROC"

# ============================================================================
# Stage 6 — set variant IDs (CHROM_POS_REF_ALT) on BOTH datasets BEFORE merging,
#   so matched variants share alleles and plink --bmerge cannot raise allele-flip
#   .missnp conflicts. plink2 lives in its own conda env.
# ============================================================================
# NOTE: reference MAF is handled PER-SUPERPOPULATION at Stage 4a (on the VCF, where the
# super-pop sample subsets exist), NOT as a pooled plink2 --maf here — a single pooled --maf
# would discard population-private ancestry-informative markers.
conda activate "$CONDA_ENV_PLINK2"
for ds in study ref; do
    plink2 --bfile "$OUTDIR/bed/${ds}" --set-all-var-ids '@_#_$r_$a' \
           --new-id-max-allele-len 100 --make-bed --out "$OUTDIR/bed/${ds}.varids" --threads "$NPROC"
done
conda deactivate

# ============================================================================
# Stage 7 — merge study + 1KG, then keep ONLY variant IDs present in BOTH
#   datasets (true intersection). Stage 7b then QCs the merged variants, LD-prune
#   follows, and finally drop (near-)all-missing individuals (ADMIXTURE rejects
#   any individual with 100% missing genotypes).
# ============================================================================
plink --bfile "$OUTDIR/bed/study.varids" --bmerge "$OUTDIR/bed/ref.varids" \
      --keep-allele-order --allow-no-sex --make-bed --out "$OUTDIR/merge/1kg_sodh_merge" --threads "$NPROC"
echo "[COUNT] Stage 7 bmerge: $(wc -l < "$OUTDIR/merge/1kg_sodh_merge.bim") variants (study $(wc -l < "$OUTDIR/bed/study.varids.bim") + ref $(wc -l < "$OUTDIR/bed/ref.varids.bim") union), $(wc -l < "$OUTDIR/merge/1kg_sodh_merge.fam") individuals"

cut -f2 "$OUTDIR/bed/study.varids.bim" | sort -u > "$OUTDIR/merge/ids_study.txt"
cut -f2 "$OUTDIR/bed/ref.varids.bim"   | sort -u > "$OUTDIR/merge/ids_ref.txt"
comm -12 "$OUTDIR/merge/ids_study.txt" "$OUTDIR/merge/ids_ref.txt" > "$OUTDIR/merge/overlap_snps.txt"
echo "[COUNT] Stage 7 study/1KG intersection (overlapping var-ids): $(wc -l < "$OUTDIR/merge/overlap_snps.txt")"
# Hard guard: plink --extract silently tolerates an empty list, so an empty intersection
# would slip through to ADMIXTURE and fail cryptically. Fail loudly here instead.
[ -s "$OUTDIR/merge/overlap_snps.txt" ] || { echo "[ERROR] study/1KG SNP intersection is empty — check contig naming + REF_MAF"; exit 1; }

plink --bfile "$OUTDIR/merge/1kg_sodh_merge" --extract "$OUTDIR/merge/overlap_snps.txt" \
      --keep-allele-order --allow-no-sex --make-bed --out "$OUTDIR/merge/1kg_sodh_merge_qced" --threads "$NPROC"
echo "[COUNT] Stage 7 intersected (qced): $(wc -l < "$OUTDIR/merge/1kg_sodh_merge_qced.bim") variants"

# ----------------------------------------------------------------------------
# Stage 7b — variant-level QC on the merged+intersected panel, BEFORE LD pruning.
#   Each filter is skipped when its knob is empty; the ancestry-correct DEFAULT enables
#   only the long-range-LD exclude (and only if LRLD_REGIONS is set). --keep-allele-order
#   is MANDATORY: a frequency-based A1/A2 flip would NOT be caught downstream because the
#   @_#_$r_$a var-ids are set once (Stage 6) and are not recomputed here.
#     * MERGE_MAF / MERGE_HWE: pooled MAF/HWE — OFF by default. On a multi-ancestry panel
#       they preferentially drop ancestry-informative markers (Wahlund / population-private
#       alleles). Reference informativeness is handled by the PER-SUPERPOP REF_MAF (Stage 4a).
#     * STUDY_COVERAGE: napkin-math --geno enforcing >=80% of STUDY samples genotyped per SNP.
#     * LRLD_REGIONS: exclude long-range high-LD blocks (MHC, chr8 inversion, ...) that
#       pairwise pruning does not fully break and that can distort ADMIXTURE structure.
# ----------------------------------------------------------------------------
QCIN="$OUTDIR/merge/1kg_sodh_merge_qced"
VQC_FLAGS=""
# Study-coverage geno: keep SNPs genotyped in >= STUDY_COVERAGE of the STUDY samples. On the
# merged bed, plink --geno is per-variant missingness over ALL merged samples; assuming the 1KG
# reference is ~fully called, the equivalent threshold is (1-STUDY_COVERAGE)*n_study/n_total
# (the napkin math from variants/K3/K3_fixed_freq_no_asia.sh). MERGE_GENO overrides it if set.
GENO_THR="$MERGE_GENO"
if [ -z "$GENO_THR" ] && [ -n "$STUDY_COVERAGE" ]; then
    n_study=$(wc -l < "$OUTDIR/bed/study.varids.fam")
    n_total=$(wc -l < "$QCIN.fam")
    GENO_THR=$(awk -v ne="$n_study" -v nt="$n_total" -v p="$STUDY_COVERAGE" \
                   'BEGIN{ if (nt>0) printf "%.10f", (1-p)*ne/nt }')
    echo "[INFO] Stage 7b study-coverage geno: >=${STUDY_COVERAGE} of $n_study study / $n_total total samples -> plink --geno $GENO_THR"
fi
[ -n "$GENO_THR" ]   && VQC_FLAGS="$VQC_FLAGS --geno $GENO_THR"
[ -n "$MERGE_MAF" ]  && VQC_FLAGS="$VQC_FLAGS --maf $MERGE_MAF"
[ -n "$MERGE_HWE" ]  && VQC_FLAGS="$VQC_FLAGS --hwe $MERGE_HWE"
if [ -n "$LRLD_REGIONS" ]; then
    [ -s "$LRLD_REGIONS" ] || { echo "[ERROR] LRLD_REGIONS set but missing/empty: $LRLD_REGIONS"; exit 1; }
    # The plinkQC GRCh38 file is chr-prefixed; strip "chr" so the range chromosomes match our
    # bare contigs (1..22) in the .bim, else plink --exclude range silently matches nothing.
    LRLD_BARE="$OUTDIR/merge/lrld_regions.bare.txt"
    awk '{ sub(/^chr/,"",$1); print }' "$LRLD_REGIONS" > "$LRLD_BARE"
    VQC_FLAGS="$VQC_FLAGS --exclude range $LRLD_BARE"
fi
if [ -n "$VQC_FLAGS" ]; then
    echo "[INFO] Stage 7b variant QC:$VQC_FLAGS"
    plink --bfile "$QCIN" $VQC_FLAGS --keep-allele-order --allow-no-sex \
          --make-bed --out "$OUTDIR/merge/1kg_sodh_merge_vqc" --threads "$NPROC"
    QCIN="$OUTDIR/merge/1kg_sodh_merge_vqc"
    nvar=$(wc -l < "$QCIN.bim")
    echo "[COUNT] Stage 7b after variant QC: $nvar variants (was $(wc -l < "$OUTDIR/merge/1kg_sodh_merge_qced.bim"))"
    [ "$nvar" -ge 1 ] || { echo "[ERROR] Stage 7b removed all variants — loosen STUDY_COVERAGE/MERGE_MAF/MERGE_HWE/LRLD_REGIONS"; exit 1; }
else
    echo "[INFO] Stage 7b variant QC: no merged-panel filters enabled"
fi

plink --bfile "$QCIN" --indep-pairwise "$LD_WINDOW" "$LD_STEP" "$LD_R2" \
      --allow-no-sex --out "$OUTDIR/pruned/indep" --threads "$NPROC"
plink --bfile "$QCIN" --extract "$OUTDIR/pruned/indep.prune.in" \
      --keep-allele-order --allow-no-sex --make-bed --out "$OUTDIR/pruned/prunedData_all" --threads "$NPROC"
echo "[COUNT] Stage 7b after LD prune: $(wc -l < "$OUTDIR/pruned/prunedData_all.bim") variants (pruned out $(wc -l < "$OUTDIR/pruned/indep.prune.out"))"

# Drop individuals that are (near-)all-missing at the pruned SNP set. ADMIXTURE aborts
# if ANY individual has 100% missing genotypes — RNA-seq cell lines with little/no
# coverage at the common-SNP set (including any sample whose VCF had 0 PASS SNPs) are
# the cause. --mind removes only those and lists them in prunedData.irem.
plink --bfile "$OUTDIR/pruned/prunedData_all" --mind "$MIND" \
      --keep-allele-order --allow-no-sex --make-bed --out "$OUTDIR/pruned/prunedData" --threads "$NPROC"
if [ -s "$OUTDIR/pruned/prunedData.irem" ]; then
    echo "[WARN] dropped $(wc -l < "$OUTDIR/pruned/prunedData.irem") individual(s) (> ${MIND} missing at the pruned set):"
    awk '{print $2}' "$OUTDIR/pruned/prunedData.irem" | paste -sd' ' -
fi
echo "[COUNT] FINAL prunedData: $(wc -l < "$OUTDIR/pruned/prunedData.bim") variants, $(wc -l < "$OUTDIR/pruned/prunedData.fam") individuals -> ADMIXTURE"

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
