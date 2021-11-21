#!/bin/bash                                                                #set a job name
#SBATCH --job-name=coverage
#SBATCH --output=./err-out/coverage.%j.out
#SBATCH --error=./err-out/coverage.%j.err
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

# Run in directory with bamfiles

sample=${1}

touch cov/covSummary.genomecov.txt
touch cov/covSummary.genomecov.rm0.txt

bedtools genomecov -d -ibam ${sample}  > cov/${sample}_cov.genomecov.txt
awk '{if($3<500) {total+=$3; ++lines}} END {print FILENAME," ",total/lines}' cov/${sample}_cov.genomecov.txt >> cov/covSummary.genomecov.txt
cat cov/${sample}_cov.genomecov.txt|awk '$3!=0{print$0}'|awk '{if($3<500) {total+=$3; ++lines}} END {print FILENAME," ",total/lines}' cov/${sample}_cov.genomecov.txt >> cov/covSummary.genomecov.rm0.txt

rm cov/${sample}_cov.genomecov.txt
