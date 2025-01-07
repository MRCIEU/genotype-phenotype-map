#!/bin/bash

PROJECT_DIR=/local-scratch/projects/genotype-phenotype-map/
PROJECT_ID=021

#TODO: this RDFS dir has problems with it.  Lots of issues around permissions?  Backing up to scratch dir for now
# RDFS_DIR=/projects/MRC_IEU/research/projects/ieu3/p1/021/working
BACKUP_DIR=$PROJECT_DIR/backup

#TODO: deal with syncing data directory later
DATA_DIR=${PROJECT_DIR}/data/
RESULTS_DIR=${PROJECT_DIR}/results/

echo "Syncing $DATA_DIR-> $BACKUP_DIR/data/"
rsync -Lavzh $DATA_DIR/ld_blocks $BACKUP_DIR/data/ --exclude ".DS_Store" --exclude "*.ipynb_checkpoints" --exclude "*.virtual_documents*"
rsync -Lavzh $DATA_DIR/study $BACKUP_DIR/data/ --exclude ".DS_Store" --exclude "*.ipynb_checkpoints" --exclude "*.virtual_documents*"

echo "Syncing $RESULTS_DIR -> $BACKUP_DIR/results/"
rsync -Lavzh $RESULTS_DIR $BACKUP_DIR/results/ --exclude ".DS_Store" --exclude "*.ipynb_checkpoints" --exclude "*.virtual_documents*"

exit 0
