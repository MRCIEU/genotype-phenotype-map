#!/bin/bash
set -e
EXTRA_ARG=$1

if [ -f .env ]
then
  export $(cat .env | xargs)
fi

export TIMESTAMP=$(date +%Y_%m_%d-%H_%M)
snakemake_log=$DATA_DIR/pipeline_metadata/logs/snakemake.log
mkdir -p $(dirname $snakemake_log)

export IMAGE=docker://andrewrrelmore/genotype_phenotype:latest
export APPTAINER_VARS="-B /local-scratch -B /projects  -B /home/$(whoami)  -B $(pwd):/home/pipeline --env TIMESTAMP=$TIMESTAMP --pwd /home/pipeline"

echo "Start time $(date)"
apptainer run $APPTAINER_VARS $IMAGE Rscript pipeline_steps/identify_studies_to_process.R &> $snakemake_log
echo "-----" &>> $snakemake_log

NUM_STUDIES=$(wc -l < ${DATA_DIR}/pipeline_metadata/studies_to_process.tsv)
if [[ $NUM_STUDIES == 0 ]]; then
  echo 'Nothing to process, exiting.'
  exit 0
elif [[ $NUM_STUDIES -gt 200000 ]]; then
  echo '!!! WARNING: too many studies to ingest at one time.  This will drastically slow down snakemake, consider splitting studies into smaller chunks'
fi

apptainer run $APPTAINER_VARS $IMAGE snakemake --profile ./ $EXTRA_SNAKEMAKE_ARG &>> $snakemake_log

rm $DATA_DIR/pipeline_metadata/studies_to_process.tsv
echo "End time $(date)"
