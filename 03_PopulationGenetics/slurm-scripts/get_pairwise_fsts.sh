source ~/.bashrc

conda activate lea

pop1_list=$(awk -v N=$SLURM_ARRAY_TASK_ID 'NR == N {print $1}' ./data/pairwise.pop.list.txt)
pop2_list=$(awk -v N=$SLURM_ARRAY_TASK_ID 'NR == N {print $2}' ./data/pairwise.pop.list.txt)

pop1=$(echo ${pop1_list} | cut -f1 -d.)
pop2=$(echo ${pop2_list} | cut -f1 -d.)

vcf=../variants/bcrf.nokinerror.QC.01maf.ld02.sexFst001.vcf.gz

vcftools --gzvcf ${vcf} --weir-fst-pop ./data/${pop1_list} --weir-fst-pop ./data/${pop2_list} --out ./output/${pop1}.${pop2}
