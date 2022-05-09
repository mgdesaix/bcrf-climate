#!/bin/bash
#set a job name
#SBATCH --job-name=mpileup-rofi
#SBATCH --output=./err-out/mpileup-rofi.%A_%a.out
#SBATCH --error=./err-out/mpileup-rofi.%A_%a.err
################
#SBATCH --time=24:00:00
#SBATCH --qos=normal
#SBATCH --partition=shas
#################
# Note: 4.84G/core or task
#################
#SBATCH --mem=16G
#SBATCH --array=1-519
#################

source ~/.bashrc

conda activate /projects/mgdesaix@colostate.edu/miniconda3/envs/bioinf

set -x

intervals=$(awk -v N=$SLURM_ARRAY_TASK_ID 'NR == N {print $1}' intervals2mb-list)
# ex. interval2mb_001.list
intervalID=$(echo ${intervals} | cut -f2 -d_ | cut -f1 -d.)
# ex. 001
outname=./region-bcfs/rofi.2mb_interval_${intervalID}.samtools.bcf

bamlist=rofi-marked-clipped.list
reference=/projects/mgdesaix@colostate.edu/reference/ROFI/2022/leucosticte_australis_final_assembly.fasta

bcftools mpileup -Ou -f ${reference} -R ./intervals2mb/${intervals} --annotate FORMAT/AD,FORMAT/DP -b ${bamlist} | bcftools call -mv -Ob -o ${outname}

