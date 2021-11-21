#!/bin/bash
#set a job name
#SBATCH --job-name=map
#SBATCH --output=./err-out/map.%j.out
#SBATCH --error=./err-out/map.%j.err
#SBATCH --array=1-50
################
#SBATCH -t 24:00:00
#SBATCH --qos=normal
#SBATCH --partition=shas
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=
#################
# Note: 4.84G/core or task
#SBATCH --nodes=1
#SBATCH --ntasks-per-node 1
#################
#SBATCH --mem=4G
#################
set -x

# Load java
module load jdk/1.8.0

# Use -p flag so you don't get errors if directory already exists!
mkdir -p err-out

##Variables

# ex. 08N10200_R1.fastq.gz
fqgz=$(awk -v N=$SLURM_ARRAY_TASK_ID 'NR == N {print $1}' fastq-files.txt)
# ex. 08N10200
sample=$(echo ${fqgz} | cut -f1 -d'_')
# ex. Plate6
plate=Plate6

# RGID - read group ID
rgid=${plate}.${sample} # ex. Plate6.08N10200
# RGLB - read group library
rglb=${plate} # ex. Plate1
# RGPU - read group platform unit
rgpu=${rgid} # same as rgid
# RGSM - read group sample
rgsm=${sample} # 08N10200
# Output
output=${rgid}_RG.bam

samtools=./programs/samtools-1.9
BWA=./programs/bwa
reference=./reference/

${BWA}/bwa mem -t 2 ${reference} ./${sample}_R1.fastq.gz ./${sample}_R2.fastq.gz > ../bwa_dedup/aln_${rgid}.sam

cd ../bwa_dedup

${samtools}/samtools sort -o aln_${rgid}.bam aln_${rgid}.sam

java -jar ~/programs/picard.jar AddOrReplaceReadGroups INPUT=aln_${rgid}.bam RGID=${rgid} RGLB=${rglb} RGPL=illumina RGPU=${rgpu} RGSM=${rgsm} OUTPUT=${output} VALIDATION_STRINGENCY=SILENT

${samtools}/samtools index ${output}
