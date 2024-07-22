#!/bin/bash
set -e
echo 'TEST RUN'
EXTRA_SNAKEMAKE_ARG=$1

export TEST_RUN=test
export TIMESTAMP=test

export DATA_DIR=/local-scratch/projects/genotype-phenotype-map/$TEST_RUN/data/
export RESULTS_DIR=/local-scratch/projects/genotype-phenotype-map/$TEST_RUN/results/

set +e
rm -r $DATA_DIR/study/*
rm -r $DATA_DIR/ld_blocks/*/*
rm -r $DATA_DIR/pipeline_metadata/studies_to_process.tsv
rm -r $RESULTS_DIR/studies_processed.tsv
rm -r $RESULTS_DIR/$TIMESTAMP/*
rm -r $RESULTS_DIR/ld_blocks/*/*
set -e

snakemake_log=$DATA_DIR/pipeline_metadata/logs/snakemake.log
mkdir -p $(dirname $snakemake_log)

export IMAGE=docker://andrewrrelmore/genotype_phenotype:latest
export APPTAINER_VARS="-B /local-scratch -B /projects -B /home/$(whoami) -B $(pwd):/home/pipeline --env TIMESTAMP=$TIMESTAMP --pwd /home/pipeline"

echo "Start time $(date)"
apptainer run $APPTAINER_VARS $IMAGE Rscript pipeline_steps/identify_studies_to_process.R &> $snakemake_log
echo "-----" &>> $snakemake_log

if [[ $(wc -l < ${DATA_DIR}/pipeline_metadata/studies_to_process.tsv) == 0 ]]; then
  echo 'Nothing to process, exiting.'
  exit 0
fi

apptainer run $APPTAINER_VARS $IMAGE snakemake --profile ./ $EXTRA_SNAKEMAKE_ARG &>> $snakemake_log
rm $DATA_DIR/pipeline_metadata/studies_to_process.tsv
echo "End time $(date)"
