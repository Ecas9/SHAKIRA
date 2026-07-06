#!/usr/bin/env python3
"""
map_srr_to_ancestry.py

Resolve CCLE RNA-seq SRR files -> cell line -> ECLA ancestry, filter to AFR/EUR/AMR
superpopulations, and write a samplesheet for nf-core/rnavar (default) or sarek.

WHY THIS NEEDS INTERNET:
  The SRR -> cell-line link lives only in SRA/ENA run metadata (needs live
  internet). Run this on an internet-connected node (e.g. an HPC login node).

IDENTIFIER CHAIN:
  SRR (your filename) --[SRA/ENA run metadata]--> cell line name
  cell line name      --[ccle_ancestry_master.csv]--> AFR/EUR/AMR/... + COSMIC_ID
  ECLA keys on COSMIC_ID; cell line NAME is the bridge from your SRR-named files.

SAMPLE NAMES: rnavar 1.3.0 schema requires sample to match ^\\S+$ (no spaces).
  Cell line names with spaces/slashes (e.g. 'Hs 578T', 'PLC/PRF/5') are sanitized
  to [A-Za-z0-9._-] (space/slash -> '_'); collisions between DIFFERENT lines get a
  numeric suffix. The original name stays in the *_srr_ancestry.csv 'cell_line'
  column, and the sanitized id is added as 'sample' for cross-reference.

ANCESTRY TABLE (ccle_ancestry_master.csv, from Dutil et al. 2019 Table S5):
  7 ADMIXTURE clusters (Q1-Q7) collapsed to 1000G superpopulations via the
  published reference centroids:
      AFR = Q7 ; EUR = Q1+Q6 ; EAS = Q3+Q4 ; SAS = Q5 ; AMR(Amerindian) = Q2
  'dominant_superpop' = argmax; 'dominant_pct' = that share. AMR = Native-American
  component (not 1000G admixed-American); AA lines read as AFR with some EUR.

  NOTE ON AMR: because AMR here is the Amerindian/Native-American ADMIXTURE
  component, very few CCLE lines are AMR-dominant, and those that are tend to
  have lower dominant_pct (admixed). Expect small, low-confidence AMR counts;
  consider --min-pct 0 (or a low value) and inspect dominant_pct per AMR row
  before treating AMR as a balanced comparison class.

USAGE
  # A) READY-TO-RUN (SRA Run Selector metadata for PRJNA523380 -> SraRunTable.csv):
  python Validation/map_srr_to_ancestry.py \
      --ancestry Validation/ccle_ancestry_master.csv \
      --runtable Validation/SraRunTable.csv \
      --fastq-dir "$ANCESTRY_WORK/RNAseq/fasta" \
    --min-pct 0 --max-per-class 50\
      --out-prefix "$ANCESTRY_WORK/table/ccle_anc"
  # (add --min-pct 80 --max-per-class 50 for a balanced high-confidence subset)
  # (restrict classes with --superpops AFR,EUR  to reproduce the old AFR/EUR-only run)

  # B) no table on hand -> fetch ENA metadata live:
  python map_srr_to_ancestry.py --ancestry ccle_ancestry_master.csv --fetch-ena \
      --fastq-dir $SCRATCH/ancestry_validation/RNAseq/fasta --out-prefix ccle_anc

  --srr-list FILE        restrict to these SRR IDs (one per line)
  --min-pct 80           keep only confident calls (default 0 = all dominant AFR/EUR/AMR)
  --superpops AFR,EUR,AMR  superpops to include in the samplesheet (default AFR,EUR,AMR)
  --samplesheet rnavar   {rnavar (default), sarek}

OUTPUTS
  <prefix>_srr_ancestry.csv          every resolved SRR with full ancestry + QC flags
  <prefix>_rnavar_<POPS>.csv         rnavar samplesheet: columns sample,fastq_1,fastq_2
  <prefix>_unmatched.csv             SRRs whose cell line did not join (inspect these)
"""

import argparse, csv, glob, os, re, sys, urllib.request, urllib.parse
import pandas as pd

ENA_PROJECT = "PRJNA523380"
ENA_FIELDS = ["run_accession", "library_strategy", "sample_alias", "sample_title",
              "fastq_ftp", "fastq_aspera", "fastq_md5"]

def norm(x):
    if x is None or (isinstance(x, float) and pd.isna(x)):
        return None
    s = re.sub(r"[^A-Za-z0-9]", "", str(x)).upper()
    return s or None

def cellline_candidates(raw):
    """Normalized keys to try; handles CCLE 'NAME_TISSUE' by also trying pre-'_'."""
    if not raw:
        return set()
    cands = {norm(raw)}
    if "_" in str(raw):
        cands.add(norm(str(raw).split("_")[0]))
    return {c for c in cands if c}

def load_ancestry(path):
    df = pd.read_csv(path)
    key2row, ambiguous = {}, set()
    for _, r in df.iterrows():
        for k in str(r["match_keys"]).split("|"):
            if not k:
                continue
            if k in key2row and key2row[k]["cLine_ID"] != r["cLine_ID"]:
                ambiguous.add(k)
            else:
                key2row[k] = r
    for k in ambiguous:
        key2row.pop(k, None)
    return df, key2row

def fetch_ena():
    url = ("https://www.ebi.ac.uk/ena/portal/api/filereport?"
           + urllib.parse.urlencode({"accession": ENA_PROJECT, "result": "read_run",
                                     "fields": ",".join(ENA_FIELDS), "format": "tsv"}))
    sys.stderr.write(f"[INFO] fetching ENA metadata: {url}\n")
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=120) as resp:
                text = resp.read().decode()
            break
        except Exception as e:
            sys.stderr.write(f"[WARN] ENA fetch attempt {attempt+1} failed: {e}\n")
    else:
        sys.exit("[FATAL] could not reach ENA. Use --runtable with SraRunTable.csv instead.")
    out = {}
    for r in csv.DictReader(text.splitlines(), delimiter="\t"):
        if r.get("library_strategy", "").upper() not in ("RNA-SEQ", "RNASEQ", ""):
            continue
        out[r["run_accession"]] = {"cell_raw": r.get("sample_title") or r.get("sample_alias"),
                                   "fastq_ftp": r.get("fastq_ftp", ""), "fastq_md5": r.get("fastq_md5", "")}
    sys.stderr.write(f"[INFO] ENA returned {len(out)} RNA-seq runs\n")
    return out

def load_runtable(path):
    df = pd.read_csv(path)
    run_col = next((c for c in df.columns if c.lower() in ("run", "run_accession")), None)
    if not run_col:
        sys.exit("[FATAL] no Run column found in runtable")
    cell_col = next((c for c in df.columns
                     if c.lower() in ("cell_line", "cell line", "cell_line_name",
                                      "source_name", "sample name", "cell_type")), None)
    if not cell_col:
        sys.exit(f"[FATAL] no cell-line column found. Columns: {list(df.columns)}")
    sys.stderr.write(f"[INFO] runtable: run='{run_col}' cellline='{cell_col}'\n")
    return {str(r[run_col]): {"cell_raw": r[cell_col], "fastq_ftp": "", "fastq_md5": ""}
            for _, r in df.iterrows()}

def srrs_from_fastqdir(d):
    found = set()
    for p in glob.glob(os.path.join(d, "*.fastq.gz")) + glob.glob(os.path.join(d, "*.fq.gz")):
        m = re.match(r"(SRR\d+)", os.path.basename(p))
        if m:
            found.add(m.group(1))
    return found

def fastq_pair(d, srr):
    def pick(suffixes):
        for s in suffixes:
            hits = glob.glob(os.path.join(d, f"{srr}{s}"))
            if hits:
                return sorted(hits)[0]
        return None
    return (pick(["_1.fastq.gz", "_R1.fastq.gz", "_1.fq.gz"]),
            pick(["_2.fastq.gz", "_R2.fastq.gz", "_2.fq.gz"]))

# --- sample-name sanitizer (rnavar needs ^\S+$; we also drop slashes; keep stable) ---
_sample_map, _used = {}, set()
def safe_sample(cell_line):
    cl = str(cell_line)
    if cl in _sample_map:
        return _sample_map[cl]
    base = re.sub(r"[^A-Za-z0-9._-]+", "_", cl).strip("_") or "sample"
    name, i = base, 2
    while name in _used:            # only DIFFERENT cell lines reach here (same one is cached)
        name = f"{base}_{i}"; i += 1
    _used.add(name); _sample_map[cl] = name
    return name

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ancestry", required=True)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--runtable")
    src.add_argument("--fetch-ena", action="store_true")
    ap.add_argument("--fastq-dir")
    ap.add_argument("--srr-list")
    ap.add_argument("--min-pct", type=float, default=0.0)
    ap.add_argument("--max-per-class", type=int, default=0,
                    help="cap samples per superpop in the samplesheet (0=no cap); keeps highest dominant_pct")
    ap.add_argument("--superpops", default="AFR,EUR,AMR",
                    help="comma-separated superpops to include in the samplesheet "
                         "(default AFR,EUR,AMR; use AFR,EUR to reproduce the old run)")
    ap.add_argument("--samplesheet", choices=["rnavar", "sarek"], default="rnavar")
    ap.add_argument("--strandedness", default="reverse",
                    help="(ignored; rnavar 1.3.0 dropped the strandedness column. kept so old commands still parse)")
    ap.add_argument("--out-prefix", default="ccle_anc")
    a = ap.parse_args()
    _od = os.path.dirname(a.out_prefix)
    if _od:
        os.makedirs(_od, exist_ok=True)

    superpops = [s.strip().upper() for s in a.superpops.split(",") if s.strip()]
    valid_pops = {"AFR", "EUR", "EAS", "SAS", "AMR"}
    bad = [p for p in superpops if p not in valid_pops]
    if bad:
        sys.exit(f"[FATAL] unknown superpop(s) {bad}; valid: {sorted(valid_pops)}")
    if not superpops:
        sys.exit("[FATAL] --superpops resolved to an empty list")

    anc_df, key2row = load_ancestry(a.ancestry)
    srr_map = load_runtable(a.runtable) if a.runtable else fetch_ena()

    target = None
    if a.srr_list:
        target = {l.strip() for l in open(a.srr_list) if l.strip()}
    elif a.fastq_dir:
        target = srrs_from_fastqdir(a.fastq_dir)
        sys.stderr.write(f"[INFO] {len(target)} SRR IDs globbed from {a.fastq_dir}\n")
    if target:
        srr_map = {k: v for k, v in srr_map.items() if k in target}
        miss = target - set(srr_map)
        if miss:
            sys.stderr.write(f"[WARN] {len(miss)} of your SRRs had no metadata row: {sorted(miss)[:8]}...\n")

    rows, unmatched = [], []
    for srr, meta in sorted(srr_map.items()):
        hit = None
        for cand in cellline_candidates(meta["cell_raw"]):
            if cand in key2row:
                hit = key2row[cand]; break
        if hit is None:
            unmatched.append({"SRR": srr, "cell_raw": meta["cell_raw"]}); continue
        r1, r2 = fastq_pair(a.fastq_dir, srr) if a.fastq_dir else (None, None)
        rows.append({
            "SRR": srr, "sample": safe_sample(hit["cLine_ID"]),
            "cell_line": hit["cLine_ID"], "cell_raw": meta["cell_raw"],
            "COSMIC_ID": hit.get("COSMIC_ID"), "RRID": hit.get("RRID"),
            "sex": hit.get("sarek_sex", "NA"), "tissue": hit.get("Tissue"),
            "AFR": hit["AFR"], "EUR": hit["EUR"], "EAS": hit["EAS"], "SAS": hit["SAS"], "AMR": hit["AMR"],
            "dominant_superpop": hit["dominant_superpop"], "dominant_pct": hit["dominant_pct"],
            "high_confidence": "Y" if hit["dominant_pct"] >= 60 else "N",
            "fastq_1": r1 or "", "fastq_2": r2 or "", "fastq_ftp_ena": meta.get("fastq_ftp", ""),
        })

    res = pd.DataFrame(rows)
    pd.DataFrame(unmatched).to_csv(f"{a.out_prefix}_unmatched.csv", index=False)
    res.to_csv(f"{a.out_prefix}_srr_ancestry.csv", index=False)
    sys.stderr.write(f"[INFO] resolved {len(res)} SRRs ; {len(unmatched)} unmatched\n")
    if len(res):
        sys.stderr.write("[INFO] dominant superpop among resolved: "
                         + str(dict(res.dominant_superpop.value_counts())) + "\n")

    # ---- samplesheet for the requested superpops (default AFR + EUR + AMR) ----
    sub = res[(res.dominant_superpop.isin(superpops)) & (res.dominant_pct >= a.min_pct)].copy()
    if a.max_per_class and len(sub):
        sub = (sub.sort_values("dominant_pct", ascending=False)
                  .groupby("dominant_superpop", group_keys=False).head(a.max_per_class))
    ss, no_pair = [], []
    for _, r in sub.iterrows():
        if not r.fastq_1 or not r.fastq_2:
            no_pair.append(r.SRR); continue
        if a.samplesheet == "rnavar":
            ss.append({"sample": r["sample"], "fastq_1": r.fastq_1, "fastq_2": r.fastq_2})
        else:  # sarek (germline: status 0)
            ss.append({"patient": r["sample"], "sex": r.sex, "status": 0,
                       "sample": r["sample"], "lane": "L001",
                       "fastq_1": r.fastq_1, "fastq_2": r.fastq_2})

    pops_tag = "_".join(superpops)
    if a.samplesheet == "rnavar":
        cols = ["sample", "fastq_1", "fastq_2"]
        out_ss = f"{a.out_prefix}_rnavar_{pops_tag}.csv"
    else:
        cols = ["patient", "sex", "status", "sample", "lane", "fastq_1", "fastq_2"]
        out_ss = f"{a.out_prefix}_sarek_{pops_tag}.csv"
    pd.DataFrame(ss, columns=cols).to_csv(out_ss, index=False)
    per_pop = ", ".join(f"{p}={int((sub.dominant_superpop == p).sum())}" for p in superpops)
    sys.stderr.write(f"[INFO] {a.samplesheet} samplesheet -> {out_ss}: {len(ss)} rows "
                     f"({per_pop})\n")
    if no_pair:
        sys.stderr.write(f"[WARN] {len(no_pair)} {pops_tag} SRRs missing _1 or _2 in --fastq-dir, "
                         f"left out: {sorted(no_pair)[:8]}...\n")

if __name__ == "__main__":
    main()