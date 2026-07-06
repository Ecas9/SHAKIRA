#!/bin/bash
###############################################################################
# prep_refs.sh - preflight for the ancestry pipeline's reference files.
#
# Does three things, all idempotent (safe to re-run):
#   (1) REQUIRE the 1000G panel + its index  - hard-fail gate.
#   (2) CHECK + TRANSFORM the GATK known-sites VCFs (dbsnp, Mills) from the stock chr-prefixed
#       bundle into the bare-contig ("_rename") files rnavar/sarek expect, reheadered to be
#       GATK sequence-dictionary compatible with the Ensembl primary-assembly FASTA.
#   (3) VALIDATE the FASTA sidecars (.fai, .dict) and the 1K_pops panel.
#
# Contig convention: the study VCFs, the 1000G panel, and the Ensembl FASTA/GTF all use BARE
# contigs (1..22,X,Y,MT). The GATK bundle uses chr-prefixed (chr1..chrM). So known-sites must be
# renamed chr->bare (chrM->MT) and restricted to the FASTA's contigs.
#
# Envs: bcftools + samtools (envs/bcftools.yaml) and gatk4 (envs/gatk.yaml) for the .dict.
# Run:  bash refs/prep_refs.sh
###############################################################################
set -uo pipefail

# ============================ EDIT BLOCK (anchors mirror config.yaml) =========
# Configure via the SAME anchors the Snakefile uses: export the env vars, OR edit the fallbacks
# after ':-'. See the top-level README. Cheaha examples are shown in comments.
PANEL_DIR="${ANCESTRY_PANEL:-/path/to/1kg_panel}"     # e.g. /data/project/<lab>/<user>/QC
WORK_ROOT="${ANCESTRY_WORK:-/path/to/ancestry_work}"  # e.g. /scratch/$USER/ancestry_validation
GATK_DIR="${ANCESTRY_GATK:-/path/to/gatk_bundle}"     # e.g. /data/project/<lab>/references/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle

KG_GRCH38="$PANEL_DIR/1kgenome.vcf.gz"                                   # REQUIRED panel
KG_PANEL="resources/1K_pops.txt"                                        # clean 2-col panel (committed)
FASTA="$WORK_ROOT/ref/FASTA/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

# GATK known-sites: SRC = stock (chr-prefixed) bundle file; DST = bare-contig "_rename" target.
DBSNP_SRC="$GATK_DIR/dbsnp_146.hg38.vcf.gz";                             DBSNP_DST="$GATK_DIR/dbsnp_146.hg38_rename.vcf.gz"
INDEL_SRC="$GATK_DIR/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz";  INDEL_DST="$GATK_DIR/Mills_and_1000G_gold_standard_rename.indels.hg38.vcf.gz"
# =============================================================================

ok(){ printf '[OK]   %s\n' "$*"; }; fix(){ printf '[FIX]  %s\n' "$*"; }
warn(){ printf '[WARN] %s\n' "$*"; }; err(){ printf '[ERROR] %s\n' "$*" >&2; }
have(){ command -v "$1" >/dev/null 2>&1; }
for t in bcftools samtools; do have "$t" || { err "$t not on PATH - conda activate the bcftools env"; exit 1; }; done

is_indexed(){ [ -f "$1.tbi" ] || [ -f "$1.csi" ]; }
is_bgzf(){ [ "$(od -An -tx1 -N4 "$1" 2>/dev/null | tr -d ' \n')" = "1f8b0804" ]; }   # BGZF magic
contigs(){ bcftools view -h "$1" 2>/dev/null | awk -F'[=,]' '/^##contig/{print $3}'; }

fail=0

# ---- (1) REQUIRED: 1000G panel + index --------------------------------------
echo "== (1) required 1000G panel =="
if [ ! -s "$KG_GRCH38" ]; then err "1000G panel missing: $KG_GRCH38  (pull it from the server)"; fail=1
else
  ok "panel present: $KG_GRCH38"
  is_bgzf   "$KG_GRCH38" && ok "panel is bgzipped" || { err "panel is not BGZF - re-bgzip: bcftools view $KG_GRCH38 -Oz -o out.vcf.gz"; fail=1; }
  is_indexed "$KG_GRCH38" && ok "panel index present (.tbi/.csi)" || { err "panel NOT indexed - run: tabix -p vcf $KG_GRCH38"; fail=1; }
  if contigs "$KG_GRCH38" | grep -q '^chr'; then warn "panel has chr-prefixed contigs - the pipeline expects bare (1..22). Rename before use."; else ok "panel contigs are bare (1..22)"; fi
  nbi=$(bcftools view -H "$KG_GRCH38" 2>/dev/null | head -20000 | awk 'length($4)==1 && length($5)==1' | wc -l)
  [ "$nbi" -gt 0 ] && ok "panel carries biallelic SNP records (sampled $nbi/20000)" || warn "no biallelic SNPs in first 20k records - check the panel"
fi

# ---- (2) GATK known-sites: check + rename chr->bare, reheader to FASTA -------
echo "== (2) GATK known-sites (rename chr->bare, reheader) =="
[ -s "${FASTA}.fai" ] || { have samtools && { samtools faidx "$FASTA" && fix "built ${FASTA}.fai"; } || warn "no ${FASTA}.fai and samtools missing"; }
KEEP=""; [ -s "${FASTA}.fai" ] && KEEP=$(cut -f1 "${FASTA}.fai" | paste -sd, -)   # FASTA contigs, comma list

rename_known_sites(){    # <src_stock_chr_vcf> <dst_bare_rename_vcf>
  local SRC="$1" DST="$2" tag; tag="$(basename "$DST")"
  if [ -s "$DST" ] && is_indexed "$DST" && ! contigs "$DST" | grep -q '^chr'; then ok "$tag already bare + indexed - skip"; return 0; fi
  if [ ! -s "$SRC" ]; then
    if [ -s "$DST" ]; then warn "$tag exists but stock SRC missing ($SRC) - cannot re-verify, leaving as-is"; return 0
    else err "cannot build $tag: stock SRC missing ($SRC). Pull the GATK bundle file."; fail=1; return 1; fi
  fi
  is_bgzf "$SRC" || { err "$SRC is not BGZF (wrong format) - cannot process"; fail=1; return 1; }
  # build old<TAB>new rename map from the contigs actually present (chrN->N, chrX/Y->X/Y, chrM->MT)
  local MAP; MAP="$(mktemp)"
  contigs "$SRC" | awk '/^chr/{o=$0; n=$0; sub(/^chr/,"",n); if(n=="M")n="MT"; print o"\t"n}' > "$MAP"
  if [ ! -s "$MAP" ]; then ok "$tag source already bare - copying to target + indexing"; bcftools view "$SRC" -Oz -o "$DST"; bcftools index -t "$DST"; return 0; fi
  fix "$tag: renaming $(wc -l < "$MAP") contigs (chr->bare) + restricting to FASTA contigs + reheader"
  local TMP; TMP="$(mktemp -u).vcf.gz"
  if [ -n "$KEEP" ]; then
    bcftools annotate --rename-chrs "$MAP" "$SRC" -Ou | bcftools view -t "$KEEP" -Oz -o "$TMP"
    bcftools reheader -f "${FASTA}.fai" -o "$DST" "$TMP"           # ##contig lines := FASTA (GATK dict-compatible)
    rm -f "$TMP"
  else
    warn "no FASTA .fai - renaming only (cannot restrict/reheader to reference contigs)"
    bcftools annotate --rename-chrs "$MAP" "$SRC" -Oz -o "$DST"
  fi
  bcftools index -t "$DST"
  contigs "$DST" | grep -q '^chr' && { err "$tag still has chr contigs after rename"; fail=1; } || ok "$tag -> bare contigs, indexed"
  [ -n "$KEEP" ] && { comm -23 <(contigs "$DST"|sort -u) <(cut -f1 "${FASTA}.fai"|sort -u) | grep -q . \
      && warn "$tag has contigs not in the FASTA (GATK may reject)" || ok "$tag contigs are a subset of the FASTA"; }
  rm -f "$MAP"
}
rename_known_sites "$DBSNP_SRC" "$DBSNP_DST"
rename_known_sites "$INDEL_SRC" "$INDEL_DST"

# ---- (3) FASTA sidecars + panel validation ----------------------------------
echo "== (3) FASTA .dict + 1K_pops panel =="
DICT="${FASTA%.*}.dict"
if [ -s "$DICT" ]; then ok "sequence dictionary present: $DICT"
elif have gatk; then gatk CreateSequenceDictionary -R "$FASTA" -O "$DICT" >/dev/null 2>&1 && fix "built $DICT" || { err "gatk CreateSequenceDictionary failed"; fail=1; }
else warn "no $DICT and gatk not on PATH - conda activate the gatk env, then: gatk CreateSequenceDictionary -R $FASTA -O $DICT"; fi
if contigs_fa=$(cut -f1 "${FASTA}.fai" 2>/dev/null) && echo "$contigs_fa" | grep -q '^chr'; then
  warn "FASTA is chr-prefixed - it should be Ensembl bare-contig to match the panel/GTF"; else ok "FASTA contigs are bare (Ensembl)"; fi
if [ -s "$KG_PANEL" ]; then
  dup=$(tail -n +2 "$KG_PANEL" | sort | uniq -d | wc -l)
  [ "$dup" -eq 0 ] && ok "panel table clean: $(($(wc -l < "$KG_PANEL")-1)) rows, no duplicate lines" \
                   || warn "panel table has $dup duplicate line(s) - dedupe: (head -1; tail -n +2 $KG_PANEL | sort -u) > clean.txt"
else err "1K_pops panel missing: $KG_PANEL"; fail=1; fi

echo "============================================================"
[ "$fail" -eq 0 ] && echo "[DONE] references OK" || { echo "[DONE] references INCOMPLETE - resolve [ERROR]s above"; exit 1; }
