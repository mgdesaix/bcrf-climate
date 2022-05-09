#!/bin/bash
#set a job name
#SBATCH --job-name=bcftools-intersect-ROFI
#SBATCH --output=./err-out/bcftools-intersect-ROFI.%A_%a.out
#SBATCH --error=./err-out/bcftools-intersect-ROFI.%A_%a.err
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

vcf1=../gatk/gatherVcfs/filtered/rofi.snps.25missing.4_5depth12.exchet.gatk.vcf.gz
vcf2=../samtools/vcfs/filtered/rofi.snps.25missing.4_5depth12.samtools.vcf.gz


bcftools isec ${vcf1} ${vcf2} -p isec2/
