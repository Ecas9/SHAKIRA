#!/bin/bash
#SBATCH --chdir=/path/to/data/K3
#SBATCH --job-name=admixture
#SBATCH --partition=express
#SBATCH --output=%x_%j.log
#SBATCH --error=%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:10:00
#SBATCH --mem-per-cpu=3G

# ADMIXTURE K — default 3; pass another value as the first arg, e.g. `sbatch admix_k3.sh 5`.
K="${1:-3}"

module load Anaconda3/2022.05
conda activate plink2

plink --bfile Final_high_quality_default\
    --remove noAsia.txt \
    --make-bed \
    --out Final_high_quality_default_no_asia

#conda deactivate

conda activate pandoc

python3 /path/to/shakira/pop_jhu.py

conda deactivate

conda activate admix2
admixture Final_high_quality_default_no_asia.bed --supervised -j8 "$K"
