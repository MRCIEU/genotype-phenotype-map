#!/bin/bash

PROJECT_DIR=/local-scratch/projects/genotype-phenotype-map/

PROJECT_ID=021
RDFS_DIR=/projects/MRC_IEU/research/projects/ieu3/p1/021/working

#TODO: deal with syncing data directory later
#DATA_DIR=${PROJECT_DIR}/data/
RESULTS_DIR=${PROJECT_DIR}/results/

#echo "Syncing $DATA_DIR-> $RDFS_DIR/data/"
#rsync -Lavzh $DATA_DIR $RDFS_DIR/data/ --exclude ".DS_Store" --exclude "*.ipynb_checkpoints" --exclude "*.virtual_documents*"

echo "Syncing $RESULTS_DIR -> $RDFS_DIR/results/"
rsync -Lavzh $RESULTS_DIR $RDFS_DIR/results/ --exclude ".DS_Store" --exclude "*.ipynb_checkpoints" --exclude "*.virtual_documents*"

exit 0
