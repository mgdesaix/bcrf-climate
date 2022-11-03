#!/bin/bash
#set a job name
#SBATCH --job-name=ROFI-merged
#SBATCH --output=./err-out/merged-ROFI.%A_%a.out
#SBATCH --error=./err-out/merged-ROFI.%A_%a.err
################
#SBATCH --time=24:00:00
#SBATCH --qos=normal
#SBATCH --partition=shas
#################
# Note: 4.84G/core or task
#################
#SBATCH --mem=16G
#SBATCH --ntasks=4
#SBATCH --array=1-147
#################

source ~/.bashrc

conda activate /projects/mgdesaix@colostate.edu/miniconda3/envs/bioinf

sample_id=$(awk -v N=$SLURM_ARRAY_TASK_ID 'NR == N {print $1}' rofi-sample-list)

merge_code="samtools merge --threads 4 ./merged/${sample_id}.merged.bam"

for i in `ls ./read_group/${sample_id}*bam`; do merge_code="${merge_code} ${i}"; done

${merge_code}

samtools index ./merged/${sample_id}.merged.bam