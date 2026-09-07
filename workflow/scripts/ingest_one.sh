#!/bin/bash
###############################################################################
# ingest_one.sh <sample> <vcf_path>
#
# Ingest ONE per-sample caller VCF (the SHAKIRA per-sample ingest step):
# harmonise contig naming to the reference panel, drop filtered records, keep
# autosomes, split multiallelics, keep strictly biallelic SNPs, force the VCF's
# internal sample name to the pipeline's sample ID, bgzip + tabix. Echoes the
# output path on stdout; the 0-SNP warning goes to stderr.
#
# A standalone script (not an exported bash function) so GNU `parallel` can invoke it
# robustly regardless of parallel's shell detection.
#
# Env (exported by the driver):
#   AUTOSOMES          1,2,...,22 (canonical, no chr prefix)
#   OUTDIR             pipeline output root
#   INGEST_THREADS     per-bcftools thread count (default 2)
#   REF_CONTIG_STYLE   contig naming of reference.kg_grch38: "nochr" (1,2,..) or
#                      "chr" (chr1,chr2,..). Default nochr. The study VCF's own
#                      style is DETECTED and converted to this. Getting this wrong
#                      is the single most common cause of an empty intersection.
#   INGEST_FILTERS     bcftools --apply-filters list. Default "PASS,." - accepts
#                      records marked PASS and records with an unset FILTER column,
#                      which is what GATK writes when no filtering step was run.
#                      Set to "PASS" to require an explicit PASS.
###############################################################################
set -euo pipefail
name="${1:?usage: ingest_one.sh <sample> <vcf_path>}"
path="${2:?usage: ingest_one.sh <sample> <vcf_path>}"
: "${AUTOSOMES:?AUTOSOMES must be exported by the driver}"
: "${OUTDIR:?OUTDIR must be exported by the driver}"
THR="${INGEST_THREADS:-2}"
REF_STYLE="${REF_CONTIG_STYLE:-nochr}"
FILTERS="${INGEST_FILTERS:-PASS,.}"

out="$OUTDIR/normalized/${name}.norm.snps.vcf.gz"
work="$OUTDIR/normalized/.${name}.tmp"
mkdir -p "$(dirname "$out")" "$work"
trap 'rm -rf "$work"' EXIT

# --- 1. detect the study VCF's contig style from its first data record ----------------
first_chrom=$(bcftools view -H -G "$path" 2>/dev/null | head -1 | cut -f1 || true)
if [ -z "$first_chrom" ]; then
    first_chrom=$(bcftools view -h -G "$path" | grep -m1 '^##contig' | sed 's/.*ID=\([^,>]*\).*/\1/' || true)
fi
case "$first_chrom" in
    chr*) IN_STYLE=chr ;;
    *)    IN_STYLE=nochr ;;
esac

# --- 2. autosome list in the STUDY VCF's own naming, so -t actually matches ------------
if [ "$IN_STYLE" = chr ]; then
    REGIONS=$(echo "$AUTOSOMES" | tr ',' '\n' | sed 's/^/chr/' | paste -sd, -)
else
    REGIONS="$AUTOSOMES"
fi

# --- 3. contig rename map, only if the study and reference styles differ --------------
RENAME=""
if [ "$IN_STYLE" != "$REF_STYLE" ]; then
    RENAME="$work/rename_chrs.txt"
    if [ "$REF_STYLE" = nochr ]; then
        echo "$AUTOSOMES" | tr ',' '\n' | awk '{print "chr"$1"\t"$1}' > "$RENAME"
    else
        echo "$AUTOSOMES" | tr ',' '\n' | awk '{print $1"\tchr"$1}' > "$RENAME"
    fi
    echo "[ingest] $name: study contigs '$IN_STYLE' -> reference '$REF_STYLE'" >&2
fi

# --- 4. sample-name harmonisation -----------------------------------------------------
# bcftools merge keys on the VCF's INTERNAL sample name, which for nf-core output is
# often a run-level ID (e.g. patient1_B30), not the pipeline's sample ID. Force it to
# $name so the merged cohort, the QC table and the ADMIXTURE .Q rows all agree, and so
# two inputs can never collide on a shared internal name.
nsamp=$(bcftools query -l "$path" | wc -l)
RH=""
if [ "$nsamp" -eq 1 ]; then
    RH="$work/sample_name.txt"
    printf '%s\n' "$name" > "$RH"
elif [ "$nsamp" -gt 1 ]; then
    echo "[WARN] $name: input has $nsamp samples; expected 1 per-sample VCF. Names left as-is." >&2
fi

# --- 5. the ingest itself -------------------------------------------------------------
# Autosome + FILTER selection happens FIRST, in the study VCF's own contig naming, so the
# rename (when needed) only has to touch the records that survive.
if [ -n "$RENAME" ]; then
    bcftools view -f "$FILTERS" -t "$REGIONS" --threads "$THR" -Ou "$path" \
      | bcftools annotate --rename-chrs "$RENAME" --threads "$THR" -Ou \
      | bcftools norm -m -any --threads "$THR" -Ou \
      | bcftools view -m2 -M2 -v snps --threads "$THR" -Oz -o "$work/body.vcf.gz"
else
    bcftools view -f "$FILTERS" -t "$REGIONS" --threads "$THR" -Ou "$path" \
      | bcftools norm -m -any --threads "$THR" -Ou \
      | bcftools view -m2 -M2 -v snps --threads "$THR" -Oz -o "$work/body.vcf.gz"
fi

if [ -n "$RH" ]; then
    bcftools reheader -s "$RH" -o "$out" "$work/body.vcf.gz"
else
    mv "$work/body.vcf.gz" "$out"
fi
bcftools index -t --threads "$THR" -f "$out"

n=$(bcftools index -n "$out")
[ "$n" -gt 0 ] || \
    echo "[WARN] $name yielded 0 biallelic autosomal SNPs after ingest (study contigs=$IN_STYLE, filters=$FILTERS) - it will fail the input QC gate" >&2
echo "$out"
