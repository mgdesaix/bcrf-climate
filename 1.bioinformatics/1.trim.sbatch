#!/bin/bash
#set a job name
#SBATCH --job-name=trim
#SBATCH --output=./err-out/trim.%j.out
#SBATCH --error=./err-out/trim.%j.err
#SBATCH --array=1-50
################
#SBATCH -t 24:00:00
#SBATCH --qos=normal
#SBATCH --partition=shas
#SBATCH --mail-type=FAIL
#################
# Note: 4.84G/core or task
#SBATCH --nodes=1
#SBATCH --ntasks-per-node 1
#################
#SBATCH --mem=4G
#################
set -x

#Run in directory where your raw .fq.gz files are and all _1.fq.gz are listed in a file fastq-files.txt

# ex. DNASERU3038_CKDL200146239-1a-N701-N507_H772MCCX2_L1_1.fq.gz
fqgz=$(awk -v N=$SLURM_ARRAY_TASK_ID 'NR == N {print $1}' fastq-files.txt)
sample=$(echo ${fqgz} | cut -f1 -d'_')

fastq=$(echo ${fqgz} | rev | cut -f2- -d'_' | rev)

# Use -p flag so you don't get errors if directory already exists!
mkdir -p err-out

###Identify program directories
FASTUNIQ=./programs/FastUniq/source # Directory for fastuniq

outdir=./bamfiles/Plate6/fastuniq-out


###Remove duplicates (potential PCR artifacts) using fastuniq

	echo ${fastq}_1.fq > pair.${sample}.txt
	echo ${fastq}_2.fq >> pair.${sample}.txt
	gunzip ${fastq}_1.fq.gz
	gunzip ${fastq}_2.fq.gz
	${FASTUNIQ}/fastuniq -i pair.${sample}.txt -t q -o ${outdir}/${sample}_R1.fastq -p ${outdir}/${sample}_R2.fastq
	gzip -f ${fastq}_1.fq
	gzip -f ${fastq}_2.fq
	gzip -f ${outdir}/${sample}_R1.fastq
	gzip -f ${outdir}/${sample}_R2.fastq
	rm pair.${sample}.txt
