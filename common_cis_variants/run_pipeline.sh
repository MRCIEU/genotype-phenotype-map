#!/bin/bash
set -e
EXTRA_ARG=$1

if [ -f .env ]
then
  export $(cat .env | xargs)
fi

#rm -r /local-scratch/projects/genotype-phenotype-map/prototype/data/study/*
#rm -r /local-scratch/projects/genotype-phenotype-map/prototype/data/study/* || echo "deleting stuff"
#rm -r /local-scratch/projects/genotype-phenotype-map/prototype/data/ld_blocks/* || echo "deleting stuff"
#rm -r /local-scratch/projects/genotype-phenotype-map/prototype/data/study/*/imput*
#rm -r /local-scratch/projects/genotype-phenotype-map/prototype/data/study/*/fine*
#rm -r /local-scratch/projects/genotype-phenotype-map/prototype/data/ld_blocks/*/*/*/fine*

export TIMESTAMP=$(date +%Y_%m_%d-%H_%M)
snakemake_log=$DATA_DIR/pipeline_metadata/logs/snakemake_log_$TIMESTAMP.log
mkdir -p $(dirname $snakemake_log)

echo "Start time $(date)"
#Rscript skip_steps.R
Rscript identify_studies_to_process.R &> $snakemake_log
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
