#!/bin/bash
###############################################################################
# get_example_data.sh - fetch a small, public CCLE RNA-seq example dataset so the
# SHAKIRA pipeline can be run end-to-end without access to any restricted cohort.
#
# WHY THIS EXISTS
#   New users need to be able to test the pipeline. The six cell lines
#   below are public CCLE RNA-seq runs (SRA study SRP186687) whose DNA-based ancestry
#   is already published in Validation/ccle_ancestry_master.csv, so the expected
#   answer is known: SHAKIRA should return >95% AFR for the three AFR-dominant lines
#   and >95% EUR for the three EUR-dominant lines.
#
#   Full CCLE runs are 7-23 Gbase each. By default this script uses HTTP range
#   requests to pull only the FIRST few hundred MB of each gzipped FASTQ and then
#   truncates to whole 4-line records, giving a ~2M-read-per-mate toy set (~1.5 GB
#   total, a few minutes on a normal connection) that is still deep enough to call
#   tens of thousands of autosomal SNPs - comfortably above the pipeline's
#   params.qc.min_snps gate. Use --full for the complete runs.
#
# USAGE
#   bash resources/get_example_data.sh [options]
#     -o, --outdir DIR    output directory                (default: example_data)
#     -n, --reads N       reads per mate to keep          (default: 2000000)
#     -f, --full          download complete FASTQs (no truncation, md5-verified)
#     -j, --jobs N        parallel downloads              (default: 3)
#     -l, --list          print the sample table and exit
#     -h, --help          this message
#
# OUTPUT
#   <outdir>/fastq/<SAMPLE>_{1,2}.fastq.gz
#   <outdir>/samplesheet.csv     ready for nf-core/rnavar (Validation/run_rnavar.slurm)
#   <outdir>/expected_ancestry.csv   published DNA-based ancestry for these lines
#
# NEXT STEPS (printed again at the end)
#   1. nextflow run nf-core/rnavar -params-file Validation/rnavar_params.yaml \
#          --input <outdir>/samplesheet.csv --outdir $ANCESTRY_WORK/rnavar
#   2. snakemake -s workflow/Snakefile --use-conda --profile profiles/slurm
#   3. compare shakira_out/analysis/ancestry_study_proportions.csv against
#      <outdir>/expected_ancestry.csv
#
# REQUIREMENTS: bash, curl, awk, gzip (all standard). No SRA Toolkit needed - FASTQs
# come straight from the ENA mirror over HTTPS.
###############################################################################
set -uo pipefail

OUTDIR="example_data"
READS=2000000
FULL=0
JOBS=3
LIST=0

# sample_id | SRA run | expected dominant super-pop | published % | cell line
SAMPLES=(
  "Raji|SRR8615972|AFR|99.4|Raji (Burkitt lymphoma)"
  "P3HR1|SRR8616151|AFR|97.8|P3HR-1 (Burkitt lymphoma)"
  "PLCPRF5|SRR8615968|AFR|95.8|PLC/PRF/5 (hepatoma)"
  "MDAMB361|SRR8615581|EUR|98.4|MDA-MB-361 (breast)"
  "WM2664|SRR8618314|EUR|99.0|WM266-4 (melanoma)"
  "VMCUB1|SRR8618300|EUR|95.6|VM-CUB-1 (bladder)"
)

usage(){ sed -n '2,40p' "$0" | sed 's/^# \{0,1\}//'; exit 0; }
while [ $# -gt 0 ]; do
  case "$1" in
    -o|--outdir) OUTDIR="$2"; shift 2;;
    -n|--reads)  READS="$2";  shift 2;;
    -f|--full)   FULL=1;      shift;;
    -j|--jobs)   JOBS="$2";   shift 2;;
    -l|--list)   LIST=1;      shift;;
    -h|--help)   usage;;
    *) echo "[example] unknown option: $1 (try --help)" >&2; exit 2;;
  esac
done

if [ "$LIST" = 1 ]; then
  printf "%-10s %-12s %-4s %7s  %s\n" SAMPLE RUN POP PCT CELL_LINE
  for r in "${SAMPLES[@]}"; do IFS='|' read -r s a p q c <<< "$r"
    printf "%-10s %-12s %-4s %6s%%  %s\n" "$s" "$a" "$p" "$q" "$c"; done
  exit 0
fi

command -v curl >/dev/null || { echo "[example][ERROR] curl not found"; exit 1; }
mkdir -p "$OUTDIR/fastq" "$OUTDIR/.tmp"

# ---------------------------------------------------------------- ENA lookup
# One API call resolves every run to its exact FASTQ URLs, so no directory-layout
# guessing and no SRA Toolkit. Fail loudly if the API is unreachable.
ACCS=$(printf "%s\n" "${SAMPLES[@]}" | cut -d'|' -f2 | paste -sd, -)
REPORT="$OUTDIR/.tmp/ena_filereport.tsv"
echo "[example] querying ENA for $ACCS"
curl -fsSL --retry 3 --retry-delay 3 \
  "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${ACCS}&result=read_run&fields=run_accession,fastq_ftp,fastq_bytes,fastq_md5&format=tsv" \
  -o "$REPORT" || { echo "[example][ERROR] could not reach the ENA API. Check network/proxy."; exit 1; }
[ "$(wc -l < "$REPORT")" -gt 1 ] || { echo "[example][ERROR] ENA returned no rows for $ACCS"; exit 1; }

url_for(){ # <run> <mate 1|2>
  awk -v r="$1" -v m="$2" -F'\t' '$1==r{ n=split($2,u,";"); if (n>=m) print "https://"u[m] }' "$REPORT"
}
md5_for(){ awk -v r="$1" -v m="$2" -F'\t' '$1==r{ n=split($4,h,";"); if (n>=m) print h[m] }' "$REPORT"; }

# Bytes to request per mate when truncating. FASTQ.gz runs ~2.5-3.5x compression and
# ~250 bytes/record uncompressed at 101bp, so ~90 bytes of gzip per read; the 2x margin
# guarantees we can still fill READS complete records after the truncated tail is dropped.
RANGE_BYTES=$(( READS * 90 * 2 ))

fetch_one(){ # <sample> <run> <mate>
  local s="$1" a="$2" m="$3"
  local out="$OUTDIR/fastq/${s}_${m}.fastq.gz"
  local url; url=$(url_for "$a" "$m")
  [ -n "$url" ] || { echo "[example][ERROR] no ENA URL for $a mate $m"; return 1; }
  if [ -s "$out" ]; then echo "[example] $s mate$m already present - skipping"; return 0; fi

  if [ "$FULL" = 1 ]; then
    echo "[example] $s mate$m: full download"
    curl -fsSL --retry 3 --retry-delay 5 "$url" -o "$out.part" || return 1
    local want got; want=$(md5_for "$a" "$m")
    if [ -n "$want" ] && command -v md5sum >/dev/null; then
      got=$(md5sum "$out.part" | cut -d' ' -f1)
      [ "$want" = "$got" ] || { echo "[example][ERROR] md5 mismatch for $s mate$m"; rm -f "$out.part"; return 1; }
    fi
    mv "$out.part" "$out"
  else
    echo "[example] $s mate$m: first $(( RANGE_BYTES / 1024 / 1024 )) MB -> $READS reads"
    curl -fsSL --retry 3 --retry-delay 5 -r "0-$((RANGE_BYTES - 1))" "$url" -o "$OUTDIR/.tmp/${s}_${m}.raw.gz" || return 1
    # gzip stops with an unexpected-EOF error on the truncated stream; that is expected,
    # so ignore its status and keep only whole 4-line records (a half-written final
    # record would make the FASTQ invalid and the two mates unequal in length).
    gzip -cd "$OUTDIR/.tmp/${s}_${m}.raw.gz" 2>/dev/null \
      | awk -v max="$READS" '
          { buf[(NR-1)%4] = $0 }
          NR%4==0 { print buf[0]; print buf[1]; print buf[2]; print buf[3]
                    if (++k >= max) exit }' \
      | gzip -c > "$out"
    rm -f "$OUTDIR/.tmp/${s}_${m}.raw.gz"
  fi
  local n; n=$(gzip -cd "$out" 2>/dev/null | wc -l)
  echo "[example] $s mate$m done: $((n / 4)) reads, $(du -h "$out" | cut -f1)"
}

# ---------------------------------------------------------------- download
pids=()
for r in "${SAMPLES[@]}"; do
  IFS='|' read -r s a p q c <<< "$r"
  for m in 1 2; do
    fetch_one "$s" "$a" "$m" &
    pids+=($!)
    while [ "$(jobs -rp | wc -l)" -ge "$JOBS" ]; do wait -n 2>/dev/null || sleep 1; done
  done
done
fail=0
for pid in "${pids[@]}"; do wait "$pid" || fail=$((fail + 1)); done
[ "$fail" -eq 0 ] || { echo "[example][ERROR] $fail download(s) failed - re-run to resume (finished files are skipped)"; exit 1; }

# ---------------------------------------------------------------- sample sheet + truth
SS="$OUTDIR/samplesheet.csv"
echo "sample,fastq_1,fastq_2,strandedness" > "$SS"
EXP="$OUTDIR/expected_ancestry.csv"
echo "sample,sra_run,cell_line,expected_dominant_superpop,expected_dominant_pct" > "$EXP"
for r in "${SAMPLES[@]}"; do
  IFS='|' read -r s a p q c <<< "$r"
  f1="$(cd "$OUTDIR/fastq" && pwd)/${s}_1.fastq.gz"
  f2="$(cd "$OUTDIR/fastq" && pwd)/${s}_2.fastq.gz"
  [ -s "$f1" ] && [ -s "$f2" ] && echo "${s},${f1},${f2},auto" >> "$SS"
  echo "${s},${a},\"${c}\",${p},${q}" >> "$EXP"
done
rmdir "$OUTDIR/.tmp" 2>/dev/null || rm -rf "$OUTDIR/.tmp"

cat <<MSG

[example] ---------------------------------------------------------------
[example] $(( $(wc -l < "$SS") - 1 )) samples ready in $OUTDIR
[example]   samplesheet:       $SS
[example]   expected ancestry: $EXP
[example]
[example] Run the pipeline:
[example]   1) variant calling
[example]      nextflow run nf-core/rnavar -r 1.3.0 -profile singularity \\
[example]        -params-file Validation/rnavar_params.yaml \\
[example]        --input $SS --outdir \$ANCESTRY_WORK/rnavar
[example]   2) ancestry
[example]      snakemake -s workflow/Snakefile --use-conda --profile profiles/slurm
[example]   3) check
[example]      column -s, -t \$ANCESTRY_WORK/shakira_out/analysis/ancestry_study_proportions.csv
[example]      column -s, -t $EXP
[example]
[example] Expected: >95% AFR for Raji / P3HR1 / PLCPRF5, >95% EUR for
[example] MDAMB361 / WM2664 / VMCUB1, under the shipped K=3 AFR/EUR/AMR config.
[example] ---------------------------------------------------------------
MSG
