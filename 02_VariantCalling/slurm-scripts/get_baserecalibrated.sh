#!/bin/bash
#set a job name
#SBATCH --job-name=base-recalibrate
#SBATCH --output=./err-out/base-recalibrate.%A_%a.out
#SBATCH --error=./err-out/base-recalibrate.%A_%a.err
################
#SBATCH --time=24:00:00
#SBATCH --qos=normal
#SBATCH --partition=shas
#################
# Note: 4.84G/core or task
#################
#SBATCH --mem=16G
#SBATCH --array=1-16
#################

source ~/.bashrc

conda activate gatk4
module load jdk

bam=$(awk -v N=$SLURM_ARRAY_TASK_ID 'NR == N {print $1}' rofi-marked-clipped-list-short)
id=$(echo ${bam} | cut -d '_' -f1)
out_table=${id}.table.out
reference=/projects/mgdesaix@colostate.edu/reference/ROFI/2022/leucosticte_australis_final_assembly.fasta

gatk BaseRecalibrator -I ../../fastqs/overlap_clipped/${bam} --known-sites rofi.gatkANDsamtools.vcf.gz -O ./tables/${out_table} -R ${reference}

gatk ApplyBQSR -R ${reference} -I ../../fastqs/overlap_clipped/${bam} --bqsr-recal-file ./tables/${out_table} -O ./out_bams/${id}.recalibrated.bam

samtools index ./out_bams/${id}.recalibrated.bam
