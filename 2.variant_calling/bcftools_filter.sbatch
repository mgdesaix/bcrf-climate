#!/bin/bash 
#
#all commands that start with SBATCH contain commands that are just used by SLURM for scheduling  
#################
#set a job name  
#SBATCH --job-name=bcftoolsFilter
#SBATCH --output=./err-out/bcftoolsFilter.%j.out
#SBATCH --error=./err-out/bcftoolsFilter.%j.err
#################
#SBATCH -t 4:00:00
#SBATCH --qos=normal
#SBATCH --partition=shas
#################
# Note: 4.84G/core or task
#SBATCH --nodes=1
#SBATCH --ntasks-per-node 1
#SBATCH --mem=4G
#SBATCH --mail-type=FAIL 


#echo commands to stdout
set -x
##################  
# filtered=AMRE.full.124.6mb.filtered.vcf
bcftools view -m 2 -M 2 --min-af 0.1 --max-af 0.9 -O z -i 'F_MISSING < 0.1 & QUAL > 30' BCRF.full.vcf.gz > BCRF.filtered.1maf.10missing.30qual.vcf.gz
bcftools index BCRF.filtered.1maf.10missing.30qual.vcf.gz
