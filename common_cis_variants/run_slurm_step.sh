#!/bin/bash
#SBATCH --partition=mrcieu
#SBATCH --job-name=genotype_phenotype_pipeline

module add apps/singularity/3.8.3

if [ -z "$FILE" ]
then
  echo "Error: 'FILE' must be set"
  exit 1
fi

singularity run -B /mnt/storage/private/mrcieu/data/ \
                -B $SCRATCH_MNT_DIR:$PROJECT_LAUNCH_DIR/scratch \
                docker://andrewrrelmore/genotype_phenotype:latest $FILE