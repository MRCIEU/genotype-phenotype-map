#!/bin/bash
set -e

IMAGE=docker://andrewrrelmore/genotype_phenotype:latest

if [ -z "$1" ] || [ -z "$PIPELINE" ]
then
  echo "Error: 'COMMAND' must be passed, and PIPELINE must be set"
  exit 1
fi
COMMAND=$1

apptainer run -B /local-scratch \
              -B $(pwd)/$PIPELINE:/home/$PIPELINE \
              -B /home/$(whoami) \
              -B /projects \
              --pwd /home/$PIPELINE/ --env-file common_cis_variants/.env \
              $IMAGE $COMMAND
