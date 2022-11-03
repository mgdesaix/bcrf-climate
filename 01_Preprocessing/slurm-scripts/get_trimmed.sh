#!/bin/bash
#set a job name
#SBATCH --job-name=ROFI-trim_galore
#SBATCH --output=./err-out/trim_galore-ROFI.%A_%a.out
#SBATCH --error=./err-out/trim_galore-ROFI.%A_%a.err
################
#SBATCH --time=1:00:00
#SBATCH --qos=normal
#SBATCH --partition=shas
#################
# Note: 4.84G/core or task
#################
#SBATCH --mem=16G
#SBATCH --array=1-448
#################

source ~/.bashrc

conda activate cutadapt
module load jdk

read1=$(awk -v N=$SLURM_ARRAY_TASK_ID 'NR == N {print $1}' fastqs-read1-list)
read2=$(echo ${read1} | sed 's/1.fq.gz/2.fq.gz/')

cd raw/

trim_galore --paired --fastqc --cores 4 ${read1} ${read2} --output_dir ../trimmed/