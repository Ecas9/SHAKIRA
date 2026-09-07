#!/usr/bin/env python3
"""
sample_qc_gate.py - SHAKIRA sample-level input QC gate.

Collects the per-sample metric rows written by sample_qc_one.sh, applies the
thresholds from config.yaml (params.qc), and emits:

    <out_table>   qc/sample_qc.tsv     every sample, every metric, PASS/FAIL + reason
    <out_pass>    qc/samples_pass.txt  sample<TAB>vcf for the samples that passed
    <out_fail>    qc/samples_fail.txt  sample<TAB>reason for the samples that did not

Rationale: an RNA-seq library with very few callable SNPs, or with coverage too low
for confident genotypes, cannot be placed against the 1000 Genomes panel. Before this
gate such a sample was carried all the way to ADMIXTURE, where it was either dropped
late by `plink --mind` or - worse - returned an ancestry vector driven by a handful of
sites. Failing it up front, loudly and with a reason, is both cheaper and honest.

Thresholds (config.yaml -> params.qc):
    min_snps        minimum retained PASS biallelic autosomal SNPs      (default 10000)
    min_median_dp   minimum median FORMAT/DP at those SNPs              (default 5)
    het_range       [lo, hi] acceptable heterozygous fraction           (default [0.05, 0.60])
    min_ts_tv       minimum transition/transversion ratio, 0 = off      (default 1.50)
    max_fail_frac   abort the run if more than this fraction fail       (default 0.50)

Usage:
    sample_qc_gate.py --metrics qc/*.metrics.tsv --samples cohort/samples.resolved \
        --out-table qc/sample_qc.tsv --out-pass qc/samples_pass.txt \
        --out-fail qc/samples_fail.txt \
        --min-snps 10000 --min-median-dp 5 --het-lo 0.05 --het-hi 0.60 \
        --min-ts-tv 1.5 --max-fail-frac 0.5
"""
import argparse
import glob
import os
import sys


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--metrics", nargs="+", required=True,
                    help="per-sample metrics TSVs from sample_qc_one.sh (globs ok)")
    ap.add_argument("--samples", required=True,
                    help="sample<TAB>vcf table (cohort/samples.resolved)")
    ap.add_argument("--out-table", required=True)
    ap.add_argument("--out-pass", required=True)
    ap.add_argument("--out-fail", required=True)
    ap.add_argument("--min-snps", type=int, default=10000)
    ap.add_argument("--min-median-dp", type=float, default=5.0)
    ap.add_argument("--het-lo", type=float, default=0.05)
    ap.add_argument("--het-hi", type=float, default=0.60)
    ap.add_argument("--min-ts-tv", type=float, default=1.50)
    ap.add_argument("--max-fail-frac", type=float, default=0.50)
    a = ap.parse_args()

    paths = []
    for p in a.metrics:
        paths.extend(sorted(glob.glob(p)) if any(c in p for c in "*?[") else [p])

    vcf_of = {}
    with open(a.samples) as fh:
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) >= 2 and f[0]:
                vcf_of[f[0]] = f[1]

    rows = []
    for p in paths:
        if not os.path.exists(p):
            continue
        with open(p) as fh:
            lines = [l.rstrip("\n") for l in fh if l.strip()]
        for ln in lines[1:]:                      # skip the header
            f = ln.split("\t")
            if len(f) >= 6:
                rows.append(dict(sample=f[0], n_snps=f[1], mean_dp=f[2],
                                 median_dp=f[3], frac_het=f[4], ts_tv=f[5]))
    rows.sort(key=lambda r: r["sample"])
    if not rows:
        sys.exit("[qc][ERROR] no per-sample metrics found in: %s" % " ".join(a.metrics))

    n_pass = 0
    with open(a.out_table, "w") as t, open(a.out_pass, "w") as pf, open(a.out_fail, "w") as ff:
        t.write("sample\tn_snps\tmean_dp\tmedian_dp\tfrac_het\tts_tv\tstatus\treason\n")
        for r in rows:
            reasons = []
            n = num(r["n_snps"])
            md = num(r["median_dp"])
            fh_ = num(r["frac_het"])
            tt = num(r["ts_tv"])

            if n is None or n < a.min_snps:
                reasons.append("n_snps=%s<%d" % (r["n_snps"], a.min_snps))
            if a.min_median_dp > 0 and (md is None or md < a.min_median_dp):
                reasons.append("median_dp=%s<%g" % (r["median_dp"], a.min_median_dp))
            if fh_ is not None and not (a.het_lo <= fh_ <= a.het_hi):
                reasons.append("frac_het=%s outside [%g,%g]" % (r["frac_het"], a.het_lo, a.het_hi))
            if a.min_ts_tv > 0 and tt is not None and tt < a.min_ts_tv:
                reasons.append("ts_tv=%s<%g" % (r["ts_tv"], a.min_ts_tv))

            status = "PASS" if not reasons else "FAIL"
            t.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
                r["sample"], r["n_snps"], r["mean_dp"], r["median_dp"],
                r["frac_het"], r["ts_tv"], status, ";".join(reasons) or "-"))
            if status == "PASS":
                n_pass += 1
                pf.write("%s\t%s\n" % (r["sample"], vcf_of.get(r["sample"], "")))
            else:
                ff.write("%s\t%s\n" % (r["sample"], ";".join(reasons)))

    n_tot = len(rows)
    n_fail = n_tot - n_pass
    frac = n_fail / float(n_tot)
    sys.stderr.write("[qc] sample gate: %d/%d passed, %d failed (%.1f%%)\n"
                     % (n_pass, n_tot, n_fail, 100 * frac))
    if n_fail:
        sys.stderr.write("[qc] failing samples are listed with reasons in %s\n" % a.out_fail)
    if n_pass == 0:
        sys.exit("[qc][ERROR] every sample failed input QC - check the caller output, "
                 "the thresholds in config.yaml (params.qc), and qc/sample_qc.tsv")
    if frac > a.max_fail_frac:
        sys.exit("[qc][ERROR] %.1f%% of samples failed input QC (limit %.1f%%). This usually "
                 "means an upstream problem (wrong caller output, wrong reference, truncated "
                 "FASTQ) rather than %d individually bad libraries. Inspect qc/sample_qc.tsv, "
                 "then either fix the input or raise params.qc.max_fail_frac deliberately."
                 % (100 * frac, 100 * a.max_fail_frac, n_fail))


if __name__ == "__main__":
    main()
