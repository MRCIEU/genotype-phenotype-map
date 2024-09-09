#!/bin/bash

LIST=$1
NUM_PARALLEL=75

parallel --jobs $NUM_PARALLEL \
  'DIR_NAME=/local-scratch/data/ukb-seq/downloads/halldorexwas/decode_data; 
  BASENAME=$(basename {} .txt.gz);
  LOG_FILE=$DIR_NAME/$BASENAME.log;
  echo "{}";
  ./01_pullfromdecode.sh {} &>> $LOG_FILE' \
  ::: $(cat ${LIST})

