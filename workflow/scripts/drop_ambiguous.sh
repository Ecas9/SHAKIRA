#!/bin/bash
###############################################################################
# drop_ambiguous.sh <ids_in> <ids_out> <dropped_out>
#
# Remove STRAND-AMBIGUOUS (palindromic) SNPs - those with A/T or C/G alleles -
# from a PLINK variant-ID list.
#
# Variant IDs are set to CHROM_POS_REF_ALT by the `varids` rule (plink2
# --set-all-var-ids '@_#_$r_$a'), so both alleles can be read directly off the ID:
# the last two underscore-separated fields are REF and ALT.
#
# Why this exists
# ---------------
# A/T and C/G SNPs are their own reverse complement, so allele identity alone does
# not tell you which strand a genotype was reported on. Merging two datasets that
# may be on different strands therefore risks silently inverting such a SNP.
#
# SHAKIRA is not exposed to that failure mode: study and reference genotypes are BOTH
# called against GRCh38 and are matched on exact CHROM_POS_REF_ALT, never on allele
# frequency, and every plink step runs --keep-allele-order, so no strand or A1/A2 flip
# is ever attempted. A variant whose alleles disagree simply fails to intersect.
# Excluding palindromic SNPs outright nonetheless removes the whole class of risk at
# negligible cost in marker count, which is why it is enabled by default
# (params.exclude_ambiguous in config.yaml). Set it to false to reproduce the
# pre-revision behaviour.
###############################################################################
set -euo pipefail
IN="${1:?usage: drop_ambiguous.sh <ids_in> <ids_out> <dropped_out>}"
OUT="${2:?usage: drop_ambiguous.sh <ids_in> <ids_out> <dropped_out>}"
DROPPED="${3:?usage: drop_ambiguous.sh <ids_in> <ids_out> <dropped_out>}"

awk -v dropf="$DROPPED" '
{
  n = split($1, f, "_")
  if (n < 4) { print $1; next }                 # not a CHROM_POS_REF_ALT id - keep as-is
  a = toupper(f[n-1]); b = toupper(f[n])
  if (length(a) != 1 || length(b) != 1) { print $1; next }   # indel-shaped id - keep
  if ((a=="A" && b=="T") || (a=="T" && b=="A") ||
      (a=="C" && b=="G") || (a=="G" && b=="C")) { print $1 > dropf; next }
  print $1
}' "$IN" > "$OUT"

: > "$DROPPED".touch && rm -f "$DROPPED".touch
[ -f "$DROPPED" ] || : > "$DROPPED"
echo "[ambig] input=$(wc -l < "$IN")  kept=$(wc -l < "$OUT")  dropped_AT_CG=$(wc -l < "$DROPPED")"
[ -s "$OUT" ] || { echo "[ambig][ERROR] every intersected SNP was strand-ambiguous - refusing to continue"; exit 1; }
