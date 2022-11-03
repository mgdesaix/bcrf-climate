source ~/.bashrc

conda activate gradientForest

set -x


mtry=$(awk -v N=$SLURM_ARRAY_TASK_ID 'NR == N {print $1}' tuning-mtry-ntree-params)
ntree=$(awk -v N=$SLURM_ARRAY_TASK_ID 'NR == N {print $2}' tuning-mtry-ntree-params)


Rscript run-gradient-forest-tuning.r ${mtry} ${ntree}