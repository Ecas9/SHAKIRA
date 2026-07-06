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
# MAX_MISSING (study call-rate, Stage 3), REF_MAF (reference-informativeness, Stage 6) and
# long-range-LD exclusion (Stage 7b). The legacy knobs remain available but default OFF.
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
# NOTE on #SBATCH lines: SLURM parses them BEFORE the shell runs, so they CANNOT
# reference the EDIT-block variables. Change them here or override at submit time:
#     sbatch --time=96:00:00 --cpus-per-task=24 cheaha/run_ancestry_cheaha.sh
###############################################################################
#SBATCH --chdir=/path/to/work/amr
#SBATCH --job-name=rnavar_admixture
#SBATCH --partition=express              # ADMIXTURE (esp. the CV sweep) dominates runtime
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=./logs/%x_%j.out
#SBATCH --error=./logs/%x_%j.err
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
#SAMPLE_SHEET="/path/to/work/table/samples.tsv"   # written (discover) | read (sheet)
SAMPLE_SHEET=""
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
OUTDIR="/path/to/work/admixture_run/K3"

# --- analysis parameters -----------------------------------------------------
K_SUPERVISED=5          # supervised ADMIXTURE K = number of 1KG super-populations
K_MIN=1                 # cross-validation sweep lower bound
K_MAX=3                 # cross-validation sweep upper bound
MAX_MISSING=0.8         # keep sites genotyped in >= this fraction of samples (80%). For RNA-seq
                        # across heterogeneous cell lines a strict 0.9 keeps only ~6k ubiquitously
                        # expressed SNPs; lower (0.8/0.5/0.2) to retain more ancestry-informative
                        # sites and fewer all-missing cells. Raise toward 0.9 for denser cohorts.
LD_WINDOW=50            # plink --indep-pairwise window
LD_STEP=10              # plink --indep-pairwise step
LD_R2=0.1               # plink --indep-pairwise r^2 threshold

# --- variant-level QC, tuned for SUPERVISED ANCESTRY (see Stage 6 / Stage 7b) ---------
# Defaults intentionally differ from generic GWAS QC; each filter is skipped when empty.
REF_MAF=0.01            # MAF applied to the 1KG REFERENCE ONLY (Stage 6). The ref bed holds
                        # only 1KG samples, so this is a true reference allele-frequency filter:
                        # it drops variants ~monomorphic in the reference, which carry NO ancestry
                        # information under a supervised model. This is the principled stand-in for
                        # the legacy *pooled* `--maf 0.05`. Set "" to disable.
LRLD_REGIONS=""         # OPTIONAL but RECOMMENDED. plink range file (CHR BP1 BP2 LABEL, in GRCh38
                        # bare-contig coords) of long-range high-LD regions (MHC, chr8 inversion,
                        # ...) to EXCLUDE before pruning — pairwise --indep-pairwise does not fully
                        # break these megabase blocks and they can distort ADMIXTURE. Left unset
                        # (step skipped) rather than ship possibly-wrong coordinates; point it at a
                        # GRCh38 high-LD-regions file to enable.
MERGE_MAF=""            # pooled MAF on study+1KG. OFF: a globally-rare SNP can be a strong AIM
                        # (common in one super-pop, ~0 elsewhere) — pooled MAF discards it, and it
                        # is dominated by 1KG counts anyway. Use REF_MAF. (Set 0.05 for legacy.)
MERGE_HWE=""            # pooled Hardy-Weinberg p-threshold. OFF: on a multi-ancestry panel the
                        # Wahlund effect makes AIMs fail HWE, so pooled --hwe strips signal. Legacy
                        # hwe 1e-6 was per-RACE-STRATUM (within one ancestry). (Set 1e-6 for legacy.)
MERGE_GENO=""           # per-variant call-rate on the merged panel. OFF: near-inert because the
                        # ~2000 fully-called 1KG samples dominate the denominator; study call-rate
                        # is controlled by MAX_MISSING (Stage 3). A loose 0.1 only guards merge dropout.

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
# ============================================================================
: > "$OUTDIR/cohort/merge.list"
while IFS=$'\t' read -r name path; do
    out="$OUTDIR/normalized/${name}.norm.snps.vcf.gz"
    echo "[ingest] $name"
    bcftools view -f PASS -t "$AUTOSOMES" --threads "$NPROC" -Ou "$path" \
      | bcftools norm -m -any --threads "$NPROC" -Ou \
      | bcftools view -m2 -M2 -v snps --threads "$NPROC" -Oz -o "$out"
    bcftools index -t --threads "$NPROC" "$out"
    [ "$(bcftools index -n "$out")" -gt 0 ] || \
        echo "[WARN] $name yielded 0 PASS biallelic SNPs — it will be all-missing (likely --mind-dropped later)"
    echo "$out" >> "$OUTDIR/cohort/merge.list"
done < "$RESOLVED"

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
# Stage 4 — harmonize the GRCh38 1KG panel: subset to the study SNP positions,
#   split multiallelics, keep biallelic SNPs (symmetric with ingest, so a 1KG
#   multiallelic that carries a biallelic SNP allele still matches a study var-id).
#   -R is index-backed, so the panel MUST carry a .tbi/.csi.
# ============================================================================
REF="$OUTDIR/ref/sodh_in_1kg_all_chroms.vcf.gz"
assert_indexed "$KG_GRCH38"

# Optionally restrict the panel to high-quality samples by EXCLUDING the admixed ones
# in KG_EXCLUDE (use the sample column). --force-samples tolerates any not in the panel.
KG_SFILTER=""
if [ -n "$KG_EXCLUDE" ]; then
    [ -s "$KG_EXCLUDE" ] || { echo "[ERROR] KG_EXCLUDE set but missing/empty: $KG_EXCLUDE — generate it via the K3 high-quality pipeline, or set KG_EXCLUDE=\"\""; exit 1; }
    awk '{print $2}' "$KG_EXCLUDE" > "$OUTDIR/ref/kg_exclude_samples.txt"
    KG_SFILTER="-S ^$OUTDIR/ref/kg_exclude_samples.txt --force-samples"
    echo "[INFO] excluding $(wc -l < "$OUTDIR/ref/kg_exclude_samples.txt") admixed 1KG sample(s) from the reference"
fi

bcftools view -R "$SNPS" $KG_SFILTER --threads "$NPROC" -Ou "$KG_GRCH38" \
  | bcftools norm -m -any --threads "$NPROC" -Ou \
  | bcftools view -m2 -M2 -v snps --threads "$NPROC" -Oz -o "$REF"
bcftools index -t --threads "$NPROC" "$REF"

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
conda activate "$CONDA_ENV_PLINK2"
for ds in study ref; do
    # REF_MAF is applied to the 1KG REFERENCE ONLY: the ref bed contains only 1KG samples,
    # so --maf here is a true reference allele-frequency filter that drops variants
    # ~monomorphic in the reference (zero ancestry information under a supervised model) —
    # the principled replacement for the legacy *pooled* --maf 0.05. NOT applied to `study`:
    # MAF across 143 heterogeneous cell lines is noisy and is not what informs supervised
    # ancestry. The downstream var-id intersection adapts to the (smaller) reference set.
    EXTRA=""
    if [ "$ds" = ref ] && [ -n "$REF_MAF" ]; then EXTRA="--maf $REF_MAF"; fi
    plink2 --bfile "$OUTDIR/bed/${ds}" --set-all-var-ids '@_#_$r_$a' \
           --new-id-max-allele-len 100 $EXTRA --make-bed --out "$OUTDIR/bed/${ds}.varids" --threads "$NPROC"
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

cut -f2 "$OUTDIR/bed/study.varids.bim" | sort -u > "$OUTDIR/merge/ids_study.txt"
cut -f2 "$OUTDIR/bed/ref.varids.bim"   | sort -u > "$OUTDIR/merge/ids_ref.txt"
comm -12 "$OUTDIR/merge/ids_study.txt" "$OUTDIR/merge/ids_ref.txt" > "$OUTDIR/merge/overlap_snps.txt"
echo "[INFO] overlapping SNPs: $(wc -l < "$OUTDIR/merge/overlap_snps.txt")"
# Hard guard: plink --extract silently tolerates an empty list, so an empty intersection
# would slip through to ADMIXTURE and fail cryptically. Fail loudly here instead.
[ -s "$OUTDIR/merge/overlap_snps.txt" ] || { echo "[ERROR] study/1KG SNP intersection is empty — check contig naming + REF_MAF"; exit 1; }

plink --bfile "$OUTDIR/merge/1kg_sodh_merge" --extract "$OUTDIR/merge/overlap_snps.txt" \
      --keep-allele-order --allow-no-sex --make-bed --out "$OUTDIR/merge/1kg_sodh_merge_qced" --threads "$NPROC"

# ----------------------------------------------------------------------------
# Stage 7b — variant-level QC on the merged+intersected panel, BEFORE LD pruning.
#   Each filter is skipped when its knob is empty; the ancestry-correct DEFAULT enables
#   only the long-range-LD exclude (and only if LRLD_REGIONS is set). --keep-allele-order
#   is MANDATORY: a frequency-based A1/A2 flip would NOT be caught downstream because the
#   @_#_$r_$a var-ids are set once (Stage 6) and are not recomputed here.
#     * MERGE_MAF / MERGE_HWE: pooled MAF/HWE — OFF by default. On a multi-ancestry panel
#       they preferentially drop ancestry-informative markers (Wahlund / population-private
#       alleles). Reference-side informativeness is handled by REF_MAF (Stage 6) instead.
#     * MERGE_GENO: near-inert (1KG sample count dominates); study call-rate is MAX_MISSING.
#     * LRLD_REGIONS: exclude long-range high-LD blocks (MHC, chr8 inversion, ...) that
#       pairwise pruning does not fully break and that can distort ADMIXTURE structure.
# ----------------------------------------------------------------------------
QCIN="$OUTDIR/merge/1kg_sodh_merge_qced"
VQC_FLAGS=""
[ -n "$MERGE_GENO" ] && VQC_FLAGS="$VQC_FLAGS --geno $MERGE_GENO"
[ -n "$MERGE_MAF" ]  && VQC_FLAGS="$VQC_FLAGS --maf $MERGE_MAF"
[ -n "$MERGE_HWE" ]  && VQC_FLAGS="$VQC_FLAGS --hwe $MERGE_HWE"
if [ -n "$LRLD_REGIONS" ]; then
    [ -s "$LRLD_REGIONS" ] || { echo "[ERROR] LRLD_REGIONS set but missing/empty: $LRLD_REGIONS"; exit 1; }
    VQC_FLAGS="$VQC_FLAGS --exclude range $LRLD_REGIONS"
fi
if [ -n "$VQC_FLAGS" ]; then
    echo "[INFO] Stage 7b variant QC:$VQC_FLAGS"
    plink --bfile "$QCIN" $VQC_FLAGS --keep-allele-order --allow-no-sex \
          --make-bed --out "$OUTDIR/merge/1kg_sodh_merge_vqc" --threads "$NPROC"
    QCIN="$OUTDIR/merge/1kg_sodh_merge_vqc"
    nvar=$(wc -l < "$QCIN.bim")
    echo "[INFO] variants after Stage 7b QC: $nvar"
    [ "$nvar" -ge 1 ] || { echo "[ERROR] Stage 7b removed all variants — loosen MERGE_GENO/MERGE_MAF/MERGE_HWE"; exit 1; }
else
    echo "[INFO] Stage 7b variant QC: no merged-panel filters enabled (ancestry-correct default)"
fi

plink --bfile "$QCIN" --indep-pairwise "$LD_WINDOW" "$LD_STEP" "$LD_R2" \
      --allow-no-sex --out "$OUTDIR/pruned/indep" --threads "$NPROC"
plink --bfile "$QCIN" --extract "$OUTDIR/pruned/indep.prune.in" \
      --keep-allele-order --allow-no-sex --make-bed --out "$OUTDIR/pruned/prunedData_all" --threads "$NPROC"

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
