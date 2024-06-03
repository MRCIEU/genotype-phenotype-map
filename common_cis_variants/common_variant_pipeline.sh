#!/bin/bash

USER=$(whoami)
export $(xargs <.env)
LOGS=$DATA_DIR/pipeline_metadata/logs
mkdir -p $LOGS

#STEP 1:
./run_slurm_step.sh 'Rscript 1_gwas_subset_list.R' &> $LOGS/1_gwas_subset_list.log || exit 1
./run_slurm_step.sh './2_extract_regions_from_vcf.sh' &> $LOGS/2_extract_regions_from_vcf.log || exit 1

for chr in $(seq 1 22); do
  ./run_slurm_step.sh "Rscirpt 4_colocalize_per_ld_block.R --chr $chr"
done

exit 0
#etc...

#If we need to run it in a slurm job in the future...
#STEP1_ID=$(sbatch --parsable \
#                   --account=$ACCOUNT_ID --mem=1G --time=04:00:00 \
#                   --nodes=1 --ntasks-per-node=1 --cpus-per-task=1 \
#                   --output=$LOGS/1_gwas_subset_list.log \
#                   --export=COMMAND="Rscript 1_gwas_subset_list.R" \
#                   run_slurm_step.sh
#)
