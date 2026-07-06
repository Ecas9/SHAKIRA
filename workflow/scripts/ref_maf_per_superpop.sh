#!/bin/bash
###############################################################################
# ref_maf_per_superpop.sh
#
# Build the 1000G REFERENCE sample + position sets for supervised ADMIXTURE:
#
#   (1) KEEP-SAMPLES  (-> <out_keep_samples>): the reference is restricted to samples whose
#       SUPERpopulation is in SUPERPOPS_OF_INTEREST (drops the whole non-of-interest
#       super-pops, e.g. EAS+SAS for a K=3 AFR/EUR/AMR run), MINUS the admixed
#       SUBpopulations in EXCLUDE_SUBPOPS (ACB/ASW = admixed AFR, MXL = admixed AMR),
#       MINUS the <98%-ancestry samples in KG_EXCLUDE. Generalizable: set
#       SUPERPOPS_OF_INTEREST="AFR EUR AMR EAS SAS" (and EXCLUDE_SUBPOPS="") for K=5.
#       Because non-of-interest super-pops are dropped by SUPERPOP membership, this also
#       correctly drops CHD (an EAS sub-pop the legacy hardcoded variants/K3/1k_filter.py
#       `pops_to_exclude` list forgot).
#
#   (2) MAF POSITIONS (-> <out_positions>, only when REF_MAF is non-empty): keep a position
#       COMMON (MAF > REF_MAF) in AT LEAST ONE super-pop of interest (a UNION, computed on
#       the keep-samples). A SNP common in AFR but rare in EUR/AMR is a real AIM and is KEPT;
#       pooled MAF would discard it.
#
# Ported / generalized from the validated K3 pipeline:
#   * variants/K3/MAF.sh             (per-super-pop two-step MAF)
#   * variants/K3/1kGenome_filter.sh (drop Asia + admixed first; the >=1-super-pop UNION)
#   * variants/K3/1k_filter.py       (pops_to_exclude: ASW/ACB/MXL admixed + EAS/SAS;
#                                     "ASW/ACB excluded due to admixture", "MXL - admixture")
#   (manuscript ref 31, Rudy 2012: 1000G places admixed African-Americans in the AFR group.)
#
# !! bcftools correctness !!  -i/-e filtering is applied BEFORE sample removal, so a single
#   `view -S pop -i 'MAF>x'` computes MAF on ALL samples, not the subset. We therefore SUBSET
#   (recomputes AC/AN on write) then FILTER MAF in a SECOND step, per super-pop.
#
# Usage:
#   ref_maf_per_superpop.sh <panel.vcf.gz> <sample_pop.tsv> <out_keep_samples> <out_positions> [workdir]
# Env:
#   SUPERPOPS_OF_INTEREST  default "AFR EUR AMR"
#   EXCLUDE_SUBPOPS        default "ACB ASW MXL"   (admixed sub-pops to hard-drop; "" to keep)
#   REF_MAF                default 0.05            ("" -> skip the MAF positions step)
#   KG_EXCLUDE             optional FID<TAB>sample file (the <98% list); column 2 is dropped
#   NPROC                  default 8
###############################################################################
set -euo pipefail

PANEL="${1:?usage: ref_maf_per_superpop.sh <panel.vcf.gz> <sample_pop.tsv> <out_keep_samples> <out_positions> [workdir]}"
SAMPLE_POP="${2:?need 1KG sample->pop table (sample <TAB> pop, header ok; e.g. 1K_pops.txt)}"
OUT_KEEP="${3:?need output keep-samples file}"
OUT_POS="${4:?need output positions file}"
WORK="${5:-$(dirname "$OUT_POS")/maf_work}"
SUPERPOPS_OF_INTEREST="${SUPERPOPS_OF_INTEREST:-AFR EUR AMR}"
EXCLUDE_SUBPOPS="${EXCLUDE_SUBPOPS:-ACB ASW MXL}"
REF_MAF="${REF_MAF:-0.05}"
NPROC="${NPROC:-8}"
mkdir -p "$WORK"

# Samples to drop (the optional <98%-ancestry list). KG_EXCLUDE is FID<TAB>sample, so the
# sample is column 2. Empty file when KG_EXCLUDE is unset — handled by the awk below.
EXCL="$WORK/exclude_samples.txt"; : > "$EXCL"
if [ -n "${KG_EXCLUDE:-}" ] && [ -s "${KG_EXCLUDE:-}" ]; then
    awk '{print $2}' "$KG_EXCLUDE" | sort -u > "$EXCL"
fi

# Build sample<TAB>superpop for the REFERENCE KEEP set: keep iff the sub-pop is NOT in
# EXCLUDE_SUBPOPS, the super-pop IS in SUPERPOPS_OF_INTEREST, and the sample is NOT in EXCL.
# Done entirely in awk (the pop->super-pop map + both filters embedded) for robustness:
#   * strips \r so a CRLF (Windows-edited) panel does not turn "YRI" into "YRI\r" -> 0 matches;
#   * reads EXCL in BEGIN via getline, which is empty-safe (the old `NR==FNR ... EXCL -` pattern
#     silently treated stdin as the exclude list when EXCL was empty, e.g. KG_EXCLUDE="");
#   * tolerant of tab/space delimiters and an optional header row; dedups via an array.
SS="$WORK/sample_superpop.tsv"
awk -v keep="$SUPERPOPS_OF_INTEREST" -v drop="$EXCLUDE_SUBPOPS" -v exf="$EXCL" '
BEGIN{
  n=split(keep,a," "); for(i=1;i<=n;i++) KEEP[a[i]]=1
  n=split(drop,a," "); for(i=1;i<=n;i++) DROP[a[i]]=1
  n=split("CHB JPT CHS CDX KHV CHD",a," ");     for(i=1;i<=n;i++) P2S[a[i]]="EAS"
  n=split("CEU TSI GBR FIN IBS",a," ");         for(i=1;i<=n;i++) P2S[a[i]]="EUR"
  n=split("YRI LWK GWD MSL ESN ASW ACB",a," "); for(i=1;i<=n;i++) P2S[a[i]]="AFR"
  n=split("MXL PUR CLM PEL",a," ");             for(i=1;i<=n;i++) P2S[a[i]]="AMR"
  n=split("GIH PJL BEB STU ITU",a," ");         for(i=1;i<=n;i++) P2S[a[i]]="SAS"
  if (exf != "") while ((getline l < exf) > 0) { gsub(/\r/,"",l); if (l!="") EX[l]=1 }
}
{ gsub(/\r/,"") }                                # CRLF safety
NR==1 && tolower($1) ~ /sample/ { next }         # skip a header row if present
(NF>=2 && $1!="" && !($1 in EX)) {
  p=$2
  if (p in DROP) next                            # hard-drop admixed sub-pop
  sup=P2S[p]
  if (sup!="" && (sup in KEEP)) seen[$1"\t"sup]=1
}
END{ for (k in seen) print k }
' "$SAMPLE_POP" | sort -k1,1 > "$SS"

cut -f1 "$SS" | sort -u > "$OUT_KEEP"
echo "[maf] reference keep-samples: $(wc -l < "$OUT_KEEP")  (super-pops: $SUPERPOPS_OF_INTEREST; dropped sub-pops: ${EXCLUDE_SUBPOPS:-none}; <98% dropped: $(wc -l < "$EXCL"))"
if [ ! -s "$OUT_KEEP" ]; then
    echo "[maf][ERROR] reference keep-samples is empty. Diagnostics for the panel ($SAMPLE_POP):"
    echo "[maf][ERROR]   rows=$(wc -l < "$SAMPLE_POP" 2>/dev/null), first 2 rows (cat -A shows ^I=tab, ^M=CR):"
    head -2 "$SAMPLE_POP" 2>/dev/null | cat -A | sed 's/^/[maf][ERROR]     /'
    echo "[maf][ERROR]   Expected 2 columns 'sample <tab/space> population' with 1000G pop codes (YRI/CEU/...)."
    echo "[maf][ERROR]   If column 2 is not a pop code (e.g. columns swapped or wrong file), fix KG_PANEL."
    exit 1
fi

# Per-super-pop MAF (skipped when REF_MAF is empty).
if [ -z "$REF_MAF" ]; then
    : > "$OUT_POS"
    echo "[maf] REF_MAF empty -> skipping per-super-pop MAF positions (keep-samples only)"
    exit 0
fi

# clear stale per-super-pop artifacts so a re-run with FEWER super-pops can't leave confusing
# leftovers next to a correct union (the union itself only ever cats this run's super-pops).
rm -f "$WORK"/*.positions
: > "$OUT_POS.union"
for sp in $SUPERPOPS_OF_INTEREST; do
    slist="$WORK/${sp}.samples"
    awk -v sp="$sp" '$2==sp{print $1}' "$SS" > "$slist"
    n=$(wc -l < "$slist")
    echo "[maf] $sp: $n reference samples (keep MAF > $REF_MAF)"
    [ "$n" -ge 1 ] || { echo "[maf][WARN] no panel samples map to $sp — skipping"; continue; }

    # STEP 1 — subset to this super-pop's samples (recomputes INFO/AC,AN for the subset).
    bcftools view -S "$slist" --force-samples --threads "$NPROC" \
        -Oz -o "$WORK/${sp}.vcf.gz" "$PANEL"
    # STEP 2 — NOW filter MAF on the subset (separate command; see header note on bcftools order),
    #          writing this super-pop's OWN named positions file so the union is auditable and
    #          NOT overwritten between super-pops.
    bcftools view -i "MAF > $REF_MAF" --threads "$NPROC" "$WORK/${sp}.vcf.gz" \
        | bcftools query -f '%CHROM\t%POS\n' > "$WORK/${sp}.positions"
    echo "[maf] $sp: $(wc -l < "$WORK/${sp}.positions") positions with MAF > $REF_MAF -> $WORK/${sp}.positions"
    cat "$WORK/${sp}.positions" >> "$OUT_POS.union"
done

# UNION = concatenate the per-super-pop named .positions files (keep a position passing in >=1
# super-pop), dedupe, natural chr order. This is a UNION, never an intersection or an overwrite.
sort -u "$OUT_POS.union" | sort -k1,1V -k2,2n > "$OUT_POS"
rm -f "$OUT_POS.union"
echo "[maf] UNION positions (MAF > $REF_MAF in >=1 of: $SUPERPOPS_OF_INTEREST): $(wc -l < "$OUT_POS")  (per-super-pop .positions kept in $WORK)"
[ -s "$OUT_POS" ] || { echo "[maf][ERROR] no positions passed per-super-pop MAF — check sample lists / threshold / panel"; exit 1; }
