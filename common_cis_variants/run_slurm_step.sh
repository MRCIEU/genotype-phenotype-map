#!/bin/bash
set -e

USER=$(whoami)
IMAGE=andrewrrelmore/genotype_phenotype:latest

if [ -z "$1" ]
then
  echo "Error: 'COMMAND' must be set"
  exit 1
fi
COMMAND=$1

sudo docker run -v /home/$USER:/home/$USER \
           -v /projects:/projects \
           --env-file .env -w /home/common_cis_variants/ \
           $IMAGE $COMMAND
