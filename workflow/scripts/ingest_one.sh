#!/bin/bash
###############################################################################
# ingest_one.sh <sample> <vcf_path>
#
# Ingest ONE rnavar HaplotypeCaller VCF (the SHAKIRA per-sample ingest step):
# drop non-PASS, keep autosomes, split multiallelics, keep strictly biallelic SNPs,
# bgzip + tabix. Echoes the output path on stdout (the driver collects these into
# cohort/merge.list); the 0-SNP warning goes to stderr.
#
# A standalone script (not an exported bash function) so GNU `parallel` can invoke it
# robustly regardless of parallel's shell detection. Env: AUTOSOMES, OUTDIR (exported by
# the driver), bcftools on PATH; INGEST_THREADS (default 2) per-bcftools thread count.
###############################################################################
set -euo pipefail
name="${1:?usage: ingest_one.sh <sample> <vcf_path>}"
path="${2:?usage: ingest_one.sh <sample> <vcf_path>}"
: "${AUTOSOMES:?AUTOSOMES must be exported by the driver}"
: "${OUTDIR:?OUTDIR must be exported by the driver}"
THR="${INGEST_THREADS:-2}"

out="$OUTDIR/normalized/${name}.norm.snps.vcf.gz"
mkdir -p "$(dirname "$out")"
bcftools view -f PASS -t "$AUTOSOMES" --threads "$THR" -Ou "$path" \
  | bcftools norm -m -any --threads "$THR" -Ou \
  | bcftools view -m2 -M2 -v snps --threads "$THR" -Oz -o "$out"
bcftools index -t --threads "$THR" "$out"
[ "$(bcftools index -n "$out")" -gt 0 ] || \
    echo "[WARN] $name yielded 0 PASS biallelic SNPs — it will be all-missing (likely --mind-dropped later)" >&2
echo "$out"
