#!/bin/bash
set -e
echo 'TEST RUN'
EXTRA_ARG=$1

export TEST_RUN=test
export TIMESTAMP=test

export DATA_DIR=/local-scratch/projects/genotype-phenotype-map/$TEST_RUN/data/
export RESULTS_DIR=/local-scratch/projects/genotype-phenotype-map/$TEST_RUN/results/
export BACKUP_DIR=/local-scratch/projects/genotype-phenotype-map/$TEST_RUN/backup/

snakemake_log=$DATA_DIR/pipeline_metadata/logs/snakemake.log
mkdir -p $(dirname $snakemake_log)


if [[ $EXTRA_ARG =~ "delete" ]]; then
  set +e
  rm -rf $DATA_DIR/study/gwas_upload/*
  rm -rf $DATA_DIR/ld_blocks/gwas_upload/*
  EXTRA_ARG=""
  set -e
fi

export IMAGE=docker://mrcieu/genotype-phenotype-map:1.0.0
export APPTAINER_VARS="-B /local-scratch -B /projects -B /home/$(whoami) -B $(pwd):/home/pipeline --env TIMESTAMP=$TIMESTAMP --pwd /home/pipeline"

echo "Start time $(date)"
echo "Testing hg38 tsv"
apptainer run $APPTAINER_VARS $IMAGE Rscript worker/pipeline_worker.R --custom_message_file ../tests/data/hg38_tsv_redis_message.json

# echo "Testing hg38 vcf"
# apptainer run $APPTAINER_VARS $IMAGE Rscript worker/pipeline_worker.R --custom_message_file ../tests/data/hg38_vcf_redis_message.json

echo "End time $(date)"
