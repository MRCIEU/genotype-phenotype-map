#!/bin/bash
set -e
EXTRA_ARG=$1

if [ -f .env ]
then
  export $(cat .env | xargs)
fi

export TIMESTAMP=$(date +%Y_%m_%d-%H_%M)
snakemake_log=$DATA_DIR/pipeline_metadata/logs/snakemake_log_$TIMESTAMP.log
mkdir -p $(dirname $snakemake_log)

echo "Start time $(date)"
Rscript identify_studies_to_process.R &> $snakemake_log
echo "-----" &>> $snakemake_log
snakemake --profile ./ $EXTRA_ARG &>> $snakemake_log
echo "End time $(date)"


#IMAGE=docker://andrewrrelmore/genotype_phenotype:latest
#apptainer run -B /local-scratch \
#              -B $(pwd)/$PIPELINE:/home/$PIPELINE \
#              -B /home/$(whoami) \
#              -B /projects \
#              --env-file .env --env TIMESTAMP=$TIMESTAMP\
#              --pwd /home/ \
#              $IMAGE $COMMAND
