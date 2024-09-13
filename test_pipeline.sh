#!/bin/bash
set -e
echo 'TEST RUN'
EXTRA_ARG=$1

#If you want to add more data to this test pipeline
#smr --beqtl-summary /some/besd_file  --extract-probe probe.txt  --query 1 --make-besd --out subset_of_besd_file

export TEST_RUN=test
export TIMESTAMP=test

export DATA_DIR=/local-scratch/projects/genotype-phenotype-map/$TEST_RUN/data/
export RESULTS_DIR=/local-scratch/projects/genotype-phenotype-map/$TEST_RUN/results/

if [[ $EXTRA_ARG =~ "delete" ]]; then
  set +e
  rm -r $DATA_DIR/pipeline_metadata/studies_to_process.tsv
  rm -r $DATA_DIR/study/*
  rm -r $DATA_DIR/ld_blocks/*/*
  rm -r $DATA_DIR/pipeline_metadata/updated_ld_blocks_to_colocalise.tsv
  rm -r $RESULTS_DIR/studies_processed.tsv
  rm -r $RESULTS_DIR/$TIMESTAMP/*
  rm -r $RESULTS_DIR/ld_blocks/*/*
  EXTRA_ARG=""
  set -e
fi

snakemake_log=$DATA_DIR/pipeline_metadata/logs/snakemake.log
mkdir -p $(dirname $snakemake_log)

export IMAGE=docker://andrewrrelmore/genotype_phenotype:latest
export APPTAINER_VARS="--nv -B /local-scratch -B /projects -B /home/$(whoami) -B $(pwd):/home/pipeline --env TIMESTAMP=$TIMESTAMP --pwd /home/pipeline"

echo "Start time $(date)"
apptainer run $APPTAINER_VARS $IMAGE Rscript pipeline_steps/identify_studies_to_process.R &> $snakemake_log
echo "-----" &>> $snakemake_log

if [[ $(wc -l < ${DATA_DIR}/pipeline_metadata/studies_to_process.tsv) == 0 ]]; then
  echo 'Nothing to process, exiting.'
  exit 0
fi

apptainer run $APPTAINER_VARS $IMAGE snakemake --profile ./ $EXTRA_ARG &>> $snakemake_log
rm $DATA_DIR/pipeline_metadata/studies_to_process.tsv
echo "End time $(date)"
