#!/bin/bash
#set a job name
#SBATCH --job-name=ROFI-bwa-map
#SBATCH --output=./err-out/bwa-map-ROFI.%A_%a.out
#SBATCH --error=./err-out/bwa-map-ROFI.%A_%a.err
################
#SBATCH --time=24:00:00
#SBATCH --qos=normal
#SBATCH --partition=shas
#################
# Note: 4.84G/core or task
#################
#SBATCH --mem=16G
#SBATCH --array=1-448
#################

source ~/.bashrc

conda activate /projects/mgdesaix@colostate.edu/miniconda3/envs/bioinf
module load jdk

picard=/projects/mgdesaix@colostate.edu/programs/picard.jar
reference=/projects/mgdesaix@colostate.edu/reference/ROFI/2022/leucosticte_australis_final_assembly.fasta

read1=$(awk -v N=$SLURM_ARRAY_TASK_ID 'NR == N {print $1}' fastqs-read1-trim-list)
# ex. 19N00143_CKDL190142002-1a-N703-N506_H723KCCX2_L7_1_val_1.fq.gz

read2=$(echo ${read1} | sed 's/1_val_1/2_val_2/')

# RGLB - read group library
rglb=$(echo ${read1} | cut -f3,4 -d_)

# RGPU - read group platform unit
rgpu=$(echo ${read1} | cut -f3 -d_)

# RGSM - read group sample
# ex. 17N04030
rgsm=$(echo ${read1} | cut -f1 -d_)

# RGPL - read group platform
rgpl=ILLUMINA

output=$(echo ${read1} | cut -f1,3,4 -d_)

bwa mem -t 4 ${reference} ./trimmed/${read1} ./trimmed/${read2} > ./mapped/${output}.sam

cd mapped/

samtools sort -o ${output}.bam ${output}.sam

rm ${output}.sam

java -jar ${picard} AddOrReplaceReadGroups I=${output}.bam RGLB=${rglb} RGPL=${rgpl} RGPU=${rgpu} RGSM=${rgsm} O=../read_group/${output}.\
RG.bam VALIDATION_STRINGENCY=SILENT
cd ../read_group
samtools index ${output}.RG.bam