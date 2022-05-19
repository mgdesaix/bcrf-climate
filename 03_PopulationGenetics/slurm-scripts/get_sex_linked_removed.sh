#!/bin/bash
#set a job name
#SBATCH --job-name=filter-sex
#SBATCH --output=./err-out/filter-sex.%j.out
#SBATCH --error=./err-out/filter-sex.%j.err
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

conda activate bioinf

set -x

# ex. fst 0.01

# > 0.01 fst snp index
fst_index=bcrf.nokinerror.snps.10missing.05maf.5_4depth12.het.sexfst_0.01_index.txt
chrom_pos=bcrf.nokinerror.snps.10missing.05maf.5_4depth12.het.chrom.pos.txt
filter_list=bcrf.nokinerror.snps.10missing.05maf.5_4depth12.het.sexfst_0.01_index.chrom.pos.txt
sed -nf <(sed 's/.*/&p/' ${fst_index}) ${chrom_pos} > ${filter_list}
sed -i 's/ /\t/g' ${filter_list}

rofi=bcrf.nokinerror.snps.10missing.05maf.5_4depth12.het.vcf.gz

# remove sex-linked snps identified FST
out_basic=bcrf.nokinerror.snps.10missing.05maf.5_4depth12.het.sexfst_0.01.vcf.gz

bcftools view -T ^${filter_list} ${rofi} -O z -o ${out_basic}
bcftools index ${out_basic}
bcftools stats ${out_basic} > ${out_basic}.stats