#!/bin/bash

LIST=$1
NUM_PARALLEL=250

count=1
for FILE in $(cat ${LIST}); do
  DIR_NAME=$(dirname $FILE)
  BASENAME=$(basename ${FILE} .txt.gz)
  LOG_FILE=$DIR_NAME/$BASENAME.log
  ./01_pullfromdecode.sh $FILE &> $LOG_FILE &

  count=$(( count + 1 ))
  [ $(( $count % $NUM_PARALLEL )) -eq 0 ]  && wait
done
wait
