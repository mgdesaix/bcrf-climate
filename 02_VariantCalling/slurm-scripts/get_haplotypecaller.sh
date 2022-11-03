#!/bin/bash
#set a job name
#SBATCH --job-name=haplotypecaller-rofi
#SBATCH --output=./err-out/haplotypecaller-rofi.%A_%a.out
#SBATCH --error=./err-out/haplotypecaller-rofi.%A_%a.err
################
#SBATCH --time=24:00:00
#SBATCH --qos=normal
#SBATCH --partition=shas
#################
# Note: 4.84G/core or task
#################
#SBATCH --mem=24G
#SBATCH --array=1-519
#################

source ~/.bashrc

conda activate gatk4
module load jdk

intervals=$(awk -v N=$SLURM_ARRAY_TASK_ID 'NR == N {print $1}' ./intervals2mb-list)
# ex. interval2mb_001.list
intervalID=$(echo ${intervals} | cut -f2 -d_ | cut -f1 -d.)
# ex. 001
outname=./haplotypecaller/rofi.2mb.interval.${intervalID}.vcf.gz

bamlist=recalibrated-bam.list
reference=/projects/mgdesaix@colostate.edu/reference/ROFI/2022/leucosticte_australis_final_assembly.fasta

gatk --java-options "-Xmx24g" HaplotypeCaller -R ${reference} -I ${bamlist} -L ../gatk/intervals2mb/${intervals} -O ${outname}
