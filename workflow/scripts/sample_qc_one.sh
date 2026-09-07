#!/bin/bash
###############################################################################
# sample_qc_one.sh <sample> <normalized_vcf.gz> <out_metrics.tsv>
#
# SHAKIRA per-sample INPUT QC gate. Computes the metrics the
# sample-level gate acts on, from the ingested per-sample VCF (PASS, autosomal,
# strictly biallelic SNPs):
#
#   n_snps      number of retained SNPs         -> too few = not enough information
#                                                  to place the sample on the 1KG panel
#   median_dp   median FORMAT/DP at those SNPs  -> low coverage = unreliable genotypes
#   mean_dp     mean FORMAT/DP
#   frac_het    heterozygous fraction           -> extreme values flag contamination
#                                                  (high) or clonal/LOH artefacts (low)
#   ts_tv       transition/transversion ratio   -> RNA-seq germline SNVs sit ~2.0-2.3;
#                                                  a collapsed ratio means noise calls
#
# Metrics only; thresholds live in config.yaml (params.qc) and are applied by
# sample_qc_gate.py. Writing them unconditionally means a failing sample still has
# an auditable QC row in qc/sample_qc.tsv.
#
# Env: bcftools on PATH.
###############################################################################
set -euo pipefail
name="${1:?usage: sample_qc_one.sh <sample> <vcf> <out>}"
vcf="${2:?usage: sample_qc_one.sh <sample> <vcf> <out>}"
out="${3:?usage: sample_qc_one.sh <sample> <vcf> <out>}"
mkdir -p "$(dirname "$out")"

n_snps=$(bcftools index -n "$vcf" 2>/dev/null || echo 0)

if [ "$n_snps" -gt 0 ]; then
    read -r mean_dp median_dp < <(
        bcftools query -f '[%DP]\n' "$vcf" 2>/dev/null \
          | grep -Ex '[0-9]+' | sort -n \
          | awk '{a[NR]=$1; s+=$1}
                 END{ if(NR==0){print "NA NA"} else {printf "%.2f %d\n", s/NR, a[int((NR+1)/2)]} }'
    )
    frac_het=$(
        bcftools query -f '[%GT]\n' "$vcf" 2>/dev/null \
          | awk '{ gsub(/\|/,"/"); if($0=="./."||$0==".")next; n++;
                   split($0,g,"/"); if(g[1]!=g[2]) h++ }
                 END{ if(n==0) print "NA"; else printf "%.4f\n", h/n }'
    )
    ts_tv=$(
        bcftools query -f '%REF\t%ALT\n' "$vcf" 2>/dev/null \
          | awk 'BEGIN{ ti["AG"]=ti["GA"]=ti["CT"]=ti["TC"]=1 }
                 { k=$1$2; if(length($1)!=1||length($2)!=1) next;
                   if(k in ti) t++; else v++ }
                 END{ if(v==0) print "NA"; else printf "%.3f\n", t/v }'
    )
else
    mean_dp=NA; median_dp=NA; frac_het=NA; ts_tv=NA
fi

printf "sample\tn_snps\tmean_dp\tmedian_dp\tfrac_het\tts_tv\n"  > "$out"
printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$name" "${n_snps:-0}" "${mean_dp:-NA}" \
       "${median_dp:-NA}" "${frac_het:-NA}" "${ts_tv:-NA}"     >> "$out"
