#!/bin/bash
#SBATCH --chdir=/path/to/data/1K
#SBATCH --job-name=1kgenome
#SBATCH --partition=medium
#SBATCH --output=./logs/%x_%j.log
#SBATCH --error=./logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --mail-user=youremail@example.org
#SBATCH --mail-type=ALL
# =============================================================================
# 1kGenome_filter.sh — build the high-quality 1KG reference panel:
#   drop Asian-pop samples + known-admixed 1KG samples, keep only positions that
#   are MAF-passing in >=1 superpopulation (AFR/AMR/EUR), -> MAF.filtered.vcf.gz.
# (Previously almost entirely commented-out exploration; the working path is now
#  uncommented. The abandoned `$NF=""` awk attempt and the debug echo are dropped;
#  the corrected `$3+1` awk is used.)
# =============================================================================
module load BCFtools
module load BEDTools

# --- inputs ------------------------------------------------------------------
onek=/path/to/data/QC/1kgenome.vcf.gz                       # 1KG panel
excludedOneKpops=/path/to/data/QC/plink2/excludedSamples.txt # known-admixed 1KG samples
cp ../K3/noAsia.txt .                                                          # Asian-pop sample table

# --- 1) drop Asian-population samples ----------------------------------------
awk -F'\t' '{print $2}' noAsia.txt > noAsia_samples.txt
bcftools view -S ^noAsia_samples.txt --threads 8 -Oz -o 1k_no_asia.vcf.gz "$onek"
bcftools index 1k_no_asia.vcf.gz

# --- 2) drop known-admixed 1KG samples ---------------------------------------
bcftools view -S ^"$excludedOneKpops" --threads 8 -Oz -o 1k_no_know_admix.vcf.gz 1k_no_asia.vcf.gz
bcftools index 1k_no_know_admix.vcf.gz

# --- 3) per-superpop MAF>0.05 positions (writes EUR/AFR/AMR.positions) --------
bash /path/to/shakira/K3/MAF.sh /path/to/data/1K/1k_no_know_admix.vcf.gz

EUR=EUR.bed
AFR=AFR.bed
AMR=AMR.bed
awk '{print $1"\t"$2"\t"$2}' EUR.positions > "$EUR"
awk '{print $1"\t"$2"\t"$2}' AFR.positions > "$AFR"
awk '{print $1"\t"$2"\t"$2}' AMR.positions > "$AMR"

# all positions present in the filtered 1KG, so the intersect can only keep shared sites
bcftools query -f '%CHROM\t%POS\n' /path/to/data/1K/1k_no_know_admix.vcf.gz > 1k.positions
awk '{print $1"\t"$2"\t"$2}' 1k.positions > 1k.bed

# --- 4) keep any 1KG position MAF-passing in >=1 superpop --------------------
# https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html
# (can't do a direct union, so -wa -wb against all superpop beds and dedupe.)
bedtools intersect -wa -wb \
    -a 1k.bed \
    -b "$EUR" "$AFR" "$AMR" \
    -sorted > shared.positions.bed

# chr,pos,pos+1 (the +1 avoids zero-length intervals / duplicate rows); then unique.
# sort -k1,1V keeps natural chr1..chr22 order (not chr1,chr10,chr11,...).
awk -F'\t' '{print $1 "\t" $2 "\t" ($3 + 1)}' shared.positions.bed > shared.positions_awked.bed
sort -k1,1V shared.positions_awked.bed | uniq > shared.positions_awked_unique.bed

# --- 5) restrict the high-quality 1KG panel to those positions ---------------
bcftools view -R shared.positions_awked_unique.bed -Oz -o /path/to/work/MAF.filtered.vcf.gz 1k_no_know_admix.vcf.gz
