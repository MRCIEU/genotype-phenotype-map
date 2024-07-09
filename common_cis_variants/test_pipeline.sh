#!/bin/bash
set -e
EXTRA_ARG=$1

export DATA_DIR=/local-scratch/projects/genotype-phenotype-map/test/data/
export RESULTS_DIR=/local-scratch/projects/genotype-phenotype-map/test/results/

export TEST_RUN=true
echo 'TEST RUN'
set +e
rm -r $DATA_DIR/study/*
rm -r $DATA_DIR/ld_blocks/*/*
rm -r $DATA_DIR/results/studies_processed.tsv
rm -r $DATA_DIR/data/pipeline_metadata/studies_to_process.tsv
set -e

export TIMESTAMP=test
snakemake_log=$DATA_DIR/pipeline_metadata/logs/snakemake_log_$TIMESTAMP.log
mkdir -p $(dirname $snakemake_log)

echo "Start time $(date)"
Rscript identify_studies_to_process.R &> $snakemake_log
echo "-----" &>> $snakemake_log
snakemake --profile ./ $EXTRA_ARG &>> $snakemake_log
echo "End time $(date)"


#IMAGE=docker://andrewrrelmore/genotype_phenotype:latest
#PIPELINE=common_cis_variants
#apptainer run -B /local-scratch \
#              -B $(pwd)/$PIPELINE:/home/$PIPELINE \
#              -B /home/$(whoami) \
#              -B /projects \
#              --env-file $PIPELINE/.env --env TIMESTAMP=$TIMESTAMP\
#              --pwd /home/$PIPELINE/ \
#              $IMAGE $COMMAND
