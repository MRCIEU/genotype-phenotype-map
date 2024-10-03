#!/bin/bash
set -e

LD_BLOCKS_DIR=$DATA_DIR/ld_block_matrices/allele_flip
THOUSAND_GENOMES=$DATA_DIR/1000genomes/allele_flip
NUM_PARALLEL=100

count=1
{
  read
  while IFS=$'\t' read -r CHR START_BP END_BP ANCESTRY; do
    PLINK_OUTPUT=$LD_BLOCKS_DIR/${ANCESTRY}/${CHR}/${START_BP}_${END_BP}
    mkdir -p $(dirname $PLINK_OUTPUT)
    RANGE_FILE=${PLINK_OUTPUT}_range_file.tmp
    echo "$CHR $START_BP $END_BP allele_flip/${ANCESTRY}/${CHR}/${START_BP}_${END_BP}" > $RANGE_FILE

    #plink1.9 --bfile $THOUSAND_GENOMES/$ANCESTRY --chr $CHR --extract range $RANGE_FILE --r square spaces --out $PLINK_OUTPUT --keep-allele-order &
    #plink1.9 --bfile $THOUSAND_GENOMES/$ANCESTRY --chr $CHR --extract range $RANGE_FILE  --out $PLINK_OUTPUT  --keep-allele-order --freq &
    #plink1.9 --bfile $THOUSAND_GENOMES/$ANCESTRY --chr $CHR --extract range $RANGE_FILE  --out $PLINK_OUTPUT  --keep-allele-order --make-just-bim &
    plink1.9 --bfile $THOUSAND_GENOMES/$ANCESTRY --chr $CHR --extract range $RANGE_FILE  --out $PLINK_OUTPUT  --keep-allele-order --make-bed &
    count=$(( count + 1 ))
    [ $(( $count % $NUM_PARALLEL )) -eq 0 ] && wait
  done
} < ../../pipeline_steps/data/ld_blocks.tsv
