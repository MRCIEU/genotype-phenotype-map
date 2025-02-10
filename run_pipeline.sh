#!/bin/bash
set -e
EXTRA_SNAKEMAKE_ARG=$1

if [ -f .env ]
then
  export $(cat .env | xargs)
fi

export TIMESTAMP=$(date +%Y_%m_%d-%H_%M)
snakemake_log=$DATA_DIR/pipeline_metadata/logs/snakemake.log
mkdir -p $(dirname $snakemake_log)

export IMAGE=docker://mrcieu/genotype-phenotype-map:1.0.0
export APPTAINER_VARS="--nv -B /local-scratch -B /projects  -B /home/$(whoami)  -B $(pwd):/home/pipeline --env TIMESTAMP=$TIMESTAMP --pwd /home/pipeline"

echo "Start time $(date)"
apptainer run $APPTAINER_VARS $IMAGE Rscript pipeline_steps/identify_studies_to_process.R &> $snakemake_log
echo "-----" &>> $snakemake_log

NUM_STUDIES=$(wc -l < ${DATA_DIR}/pipeline_metadata/studies_to_process.tsv)
if [[ $NUM_STUDIES == 0 ]]; then
  echo 'Nothing to process, exiting.'
  exit 0
elif [[ $NUM_STUDIES -gt 200000 ]]; then
  echo 'ERROR: too many studies to ingest at one time.  This will drastically slow down snakemake, consider splitting studies into smaller chunks'
  exit 0
fi

if [[ $EXTRA_SNAKEMAKE_ARG =~ "clean" ]]; then
  echo "Cleaning up leftover intermediate files from previous run"
  set +e
  rm $DATA_DIR/ld_blocks/*/*/*/*_complete
  EXTRA_SNAKEMAKE_ARG=
  set -e
fi

apptainer run $APPTAINER_VARS $IMAGE snakemake --profile ./ $EXTRA_SNAKEMAKE_ARG &>> $snakemake_log

rm $DATA_DIR/pipeline_metadata/studies_to_process.tsv
echo "End time $(date)"
