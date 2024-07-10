#!/bin/bash
set -e
echo 'TEST RUN'
EXTRA_SNAKEMAKE_ARG=$1

export DATA_DIR=/local-scratch/projects/genotype-phenotype-map/test/data/
export RESULTS_DIR=/local-scratch/projects/genotype-phenotype-map/test/results/
export IMAGE=docker://andrewrrelmore/genotype_phenotype:latest
export TEST_RUN=true

set +e
rm -r $DATA_DIR/study/*
rm -r $DATA_DIR/ld_blocks/*/*
rm -r $DATA_DIR/results/studies_processed.tsv
rm -r $DATA_DIR/pipeline_metadata/studies_to_process.tsv
set -e

export TIMESTAMP=test
snakemake_log=$DATA_DIR/pipeline_metadata/logs/snakemake_log_$TIMESTAMP.log
mkdir -p $(dirname $snakemake_log)

echo "Start time $(date)"
snakemake --profile ./ $EXTRA_SNAKEMAKE_ARG &>> $snakemake_log
echo "End time $(date)"

APPTAINER_VARS="-B glocal-scratch -B /projects  -B /home/$(whoami)  -B $(pwd):/home/pipeline  --env-file .env --env TIMESTAMP=$TIMESTAMP --pwd /home/pipeline"

apptainer run $APPTAINER_VARS $IMAGE Rscript identify_studies_to_process.R &> $snakemake_log
echo "-----" &>> $snakemake_log

if [[ $(wc -l $DATA_DIR/pipeline_metadata/studies_to_process.tsv) > 2 ]]; then
  echo 'Nothing to process, exiting.'
  exit 0
fi

apptainer run $APPTAINER_VARS $IMAGE snakemake --profile ./ $EXTRA_SNAKEMAKE_ARG &>> $snakemake_log
