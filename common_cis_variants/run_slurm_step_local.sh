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

docker run -v /Users/wt23152/Documents/Projects/genotype-phenotype-map/common_cis_variants/scratch:/user/work/$USER/scratch \
           -v /Users/wt23152/Documents/Projects/genotype-phenotype-map:/user/home/$USER \
           -v /Users/wt23152/Documents/Projects/scratch:/Users/wt23152/Documents/Projects/scratch \
           --env-file .env -w /home/common_cis_variants/ \
           $IMAGE $COMMAND
