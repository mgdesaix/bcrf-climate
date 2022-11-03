#!/bin/bash
#set a job name
#SBATCH --job-name=related
#SBATCH --output=./err-out/related.%j.out
#SBATCH --error=./err-out/related.%j.err
################
#SBATCH --time=12:00:00
#SBATCH --qos=normal
#SBATCH --partition=shas
#################
# Note: 4.84G/core or task
#################
#SBATCH --mem=32G
#################

source ~/.bashrc

conda activate /projects/mgdesaix@colostate.edu/miniconda3/envs/bioinf

set -x
vcf=/scratch/summit/mgdesaix@colostate.edu/ROFI/snps/recalibrated/gatherVcfs/bcrf.vcf.gz
bcftools view -m 2 -M 2 -O z --type snps ${vcf} -o ./data/bcrf.snps.vcf.gz

plink2 --vcf ./data/bcrf.snps.vcf.gz --make-bed --allow-extra-chr --out ./data/bcrf

king=/projects/mgdesaix@colostate.edu/programs/king
${king} -b ./data/bcrf.bed --related --degree 2
