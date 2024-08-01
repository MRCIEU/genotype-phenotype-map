#!/bin/bash

LIST=$1
NUM_PARALLEL=250

count=1
for FILE in $(cat ${LIST}); do
  ./01_pullfromdecode.sh $FILE  &

  count=$(( count + 1 ))
  [ $(( $count % $NUM_PARALLEL )) -eq 0 ]  && wait
done
wait
