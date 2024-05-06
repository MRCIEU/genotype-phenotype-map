#!/bin/bash

USER=$(whoami)
ACCOUNT_ID=$(sacctmgr show user withassoc format=account where user=$USER | grep smed)
SLURM_LOGS=/user/work/$USER/slurm_logs
mkdir -p $SLURM_LOGS

#STEP 1:
STEP1_ID=$(sbatch --parsable \
                   --account=$ACCOUNT_ID --mem=1G --time=04:00:00 \
                   --nodes=1 --ntasks-per-node=1 --cpus-per-task=1 \
                   --output=$SLURM_LOGS/1_gwas_subset_list.log \
                   --export=FILE=1_gwas_subset_list.R \
                   run_slurm_step.sh
)

#STEP 2:
STEP2_ID=$(sbatch --parsable --dependency=afterok:${STEP1_ID} --account=$ACCOUNT_ID
                  --mem=1G \
                  --time=04:00:00 \
                  --nodes=1 \
                  --ntasks-per-node=1 \
                  --cpus-per-task=1 \
                  --output=$SLURM_LOGS/2_next_step.log \
                  --export=FILE=2_next_step.R \
                  run_slurm_step.sh
)


#STEP 3:
STEP3_ID=$(sbatch --parsable --dependency=afterok:${STEP2_ID} \
                  --account=$ACCOUNT_ID --mem=1G --time=04:00:00 \
                  --nodes=1 --ntasks-per-node=1 --cpus-per-task=1 \
                  --output=$SLURM_LOGS/3_next_step.log \
                  --export=FILE=3_next_step.R \
                  run_slurm_step.sh
)

#etc...