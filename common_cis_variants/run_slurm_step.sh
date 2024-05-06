#!/bin/bash
#SBATCH --partition=mrcieu
#SBATCH --job-name=genotype_phenotype_pipeline

module add apps/singularity/3.8.3
USER=$(whoami)
SIF=/user/work/$USER/test/data/genotype_phenotype_latest.sif

if [ -z "$COMMAND" ]
then
  echo "Error: 'COMMAND' must be set"
  exit 1
fi

singularity run -B /mnt/storage/private/mrcieu/data/ -B /user/home/$USER/ -B /user/work/$USER/ $SIF $COMMAND
