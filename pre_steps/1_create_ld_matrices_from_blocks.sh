#!/bin/bash
set -e

LD_BLOCKS_DIR=$DATA_DIR/ld_block_matrices
THOUSAND_GENOMES=$DATA_DIR/1000genomes

{
  read
  while IFS=$'\t' read -r CHR START_BP END_BP ANCESTRY; do
    PLINK_OUTPUT=$LD_BLOCKS_DIR/${ANCESTRY}/${CHR}/${START_BP}_${END_BP}
    mkdir -p $(dirname $PLINK_OUTPUT)
    RANGE_FILE=$LD_BLOCKS_DIR/range_file.tmp
    echo "$CHR $START_BP $END_BP ${ANCESTRY}/${CHR}/${START_BP}_${END_BP}" > $RANGE_FILE

    plink1.9 --bfile $THOUSAND_GENOMES/$ANCESTRY --chr $CHR --extract range $RANGE_FILE --r square spaces --out $PLINK_OUTPUT --keep-allele-order
    plink1.9 --bfile $THOUSAND_GENOMES/$ANCESTRY --chr $CHR --extract range $RANGE_FILE  --out $PLINK_OUTPUT  --keep-allele-order --freq
    plink1.9 --bfile $THOUSAND_GENOMES/$ANCESTRY --chr $CHR --extract range $RANGE_FILE  --out $PLINK_OUTPUT  --keep-allele-order --make-just-bim
  done
} < data/ld_regions.tsv
