#!/bin/bash
#set a job name
#SBATCH --job-name=bcftools-filter-ROFI
#SBATCH --output=./err-out/bcftools-filter-ROFI.%A_%a.out
#SBATCH --error=./err-out/bcftools-filter-ROFI.%A_%a.err
################
#SBATCH --time=24:00:00
#SBATCH --qos=normal
#SBATCH --partition=shas
#################
# Note: 4.84G/core or task
#################
#SBATCH --mem=16G
#################

source ~/.bashrc

conda activate /projects/mgdesaix@colostate.edu/miniconda3/envs/bioinf

set -x

vcf=rofi.samtools.vcf.gz

# bcftools index ${vcf}
# bcftools stats ${vcf} > ${vcf}.stats

out_basic=./filtered/rofi.snps.25missing.samtools.vcf.gz
bcftools view -m 2 -M 2 -i 'F_MISSING < 0.25' ${vcf} -O z -o ${out_basic}
bcftools index ${out_basic}
bcftools stats ${out_basic} > ${out_basic}.stats

# mean depth = 6, so with 25% missing, 0.75 * 6 = 4.5
out_depth=./filtered/rofi.snps.25missing.4_5depth12.samtools.vcf.gz
bcftools view -i 'AVG(FORMAT/DP)>4.5 & AVG(FORMAT/DP)<12' ${out_basic} -O z -o ${out_depth}
bcftools index ${out_depth}
bcftools stats ${out_depth} > ${out_depth}.stats
