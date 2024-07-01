#!/bin/bash
set -e

EXTRA_ARG=$1
#rm -r /local-scratch/projects/genotype-phenotype-map/data/study/ukb-b-10003/*
#rm -r /local-scratch/projects/genotype-phenotype-map/prototype/data/study/*
#rm -r /local-scratch/projects/genotype-phenotype-map/prototype/data/ld_blocks/*
#rm -r /local-scratch/projects/genotype-phenotype-map/prototype/data/study/*/finemapped
#rm -r /local-scratch/projects/genotype-phenotype-map/prototype/data/ld_blocks/*/*/*/fine*

Rscript identify_studies_to_process.R
#Rscript skip_steps.R
export TIMESTAMP=$(date +%Y_%m_%d-%H_%M)

IMAGE=docker://andrewrrelmore/genotype_phenotype:latest
PIPELINE=common_cis_variants

snakemake --profile ./ $EXTRA_ARG &> /tmp/snakemake.log

#apptainer run -B /local-scratch \
#              -B $(pwd)/$PIPELINE:/home/$PIPELINE \
#              -B /home/$(whoami) \
#              -B /projects \
#              --env-file $PIPELINE/.env --env TIMESTAMP=$TIMESTAMP\
#              --pwd /home/$PIPELINE/ \
#              $IMAGE $COMMAND
