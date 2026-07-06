#!/bin/bash
#SBATCH --chdir=/path/to/data/K3_15_mill_filter
#SBATCH --job-name=admix_no_1K_issue
#SBATCH --partition=express
#SBATCH --output=./logs/%x_%j.log
#SBATCH --error=./logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=2:00:00
#SBATCH --mem-per-cpu=3G
#SBATCH --mail-user=youremail@example.org
#SBATCH --mail-type=ALL

module load Anaconda3/2022.05
conda activate plink2

# ADMIXTURE K — default 3; pass another value as the first arg, e.g. `sbatch exclude_admixed_1K_samples.sh 5`.
K="${1:-3}"

start_samples=Final_high_quality_default_no_asia


plink --bfile $start_samples\
    --remove 1K_98_co_exclude.txt\
    --make-bed \
    --out Final_no_admixed_1K_samples

conda deactivate

conda activate pandoc

python3 /path/to/shakira/pop_jhu.py

conda deactivate

conda activate admix2
admixture Final_no_admixed_1K_samples.bed -j8 "$K"

conda deactivate
conda activate plink2

mkdir -p PCA #-p makessure we don't overwrite 

plink\
    --bfile Final_no_admixed_1K_samples\
    --pca \
    --allow-extra-chr \
    --out PCA/Final_no_admixed_1K_samples_PCA

