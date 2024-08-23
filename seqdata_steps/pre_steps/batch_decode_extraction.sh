#!/bin/bash

LIST=$1
NUM_PARALLEL=30

DIR_NAME=/local-scratch/data/ukb-seq/downloads/halldorexwas/decode_data

count=1
for FILE in $(cat ${LIST}); do
  BASENAME=$(basename ${FILE} .txt.gz)
  LOG_FILE=$DIR_NAME/$BASENAME.log
  echo "processing $FILE"
  ./01_pullfromdecode.sh $FILE &> $LOG_FILE &

  count=$(( count + 1 ))
  [ $(( $count % $NUM_PARALLEL )) -eq 0 ]  && wait
done
wait
