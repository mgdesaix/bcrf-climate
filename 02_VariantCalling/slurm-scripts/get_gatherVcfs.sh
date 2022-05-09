#!/bin/bash
#set a job name
#SBATCH --job-name=gatherVcfs-rofi
#SBATCH --output=./err-out/gatherVcfs-rofi.%A_%a.out
#SBATCH --error=./err-out/gatherVcfs-rofi.%A_%a.err
################
#SBATCH --time=9:00:00
#SBATCH --qos=normal
#SBATCH --partition=shas
#################
# Note: 4.84G/core or task
#################
#SBATCH --mem=32G
#################

source ~/.bashrc

conda activate gatk4
module load jdk

set -x

# list file with the 519 2mb intervals
vcfs=rofi-intervals2mb-vcfs.list
# ex. 001
outname=./gatherVcfs/rofi.vcf.gz

gatk --java-options "-Xmx32g" GatherVcfs -I ${vcfs} -O ${outname}
