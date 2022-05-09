#!/bin/bash
#set a job name
#SBATCH --job-name=dup-remove
#SBATCH --output=./err-out/dup-remove.%A_%a.out
#SBATCH --error=./err-out/dup-remove.%A_%a.err
################
#SBATCH --time=24:00:00
#SBATCH --qos=normal
#SBATCH --partition=shas
#################
# Note: 4.84G/core or task
#################
#SBATCH --mem=64G
#SBATCH --array=1-147
#################

source ~/.bashrc

conda activate gatk4
module load jdk

bam=$(awk -v N=$SLURM_ARRAY_TASK_ID 'NR == N {print $1}' rofi-merged-list)
id=$(echo ${bam} | cut -d '.' -f1)
marked=$(echo ${id}_marked.bam)
clipped=$(echo ${id}_marked_clipped.bam)

gatk MarkDuplicatesSpark -I ./merged/${bam} -O ./marked/${marked} --remove-sequencing-duplicates --tmp-dir ~/scratch/tmp

bam clipOverlap --in ./marked/${marked} --out ./overlap_clipped/${clipped} --stats
