#!/bin/bash
set -e

LD_BLOCKS_DIR=$DATA_DIR/ld_block_matrices

{
  read
  while IFS=$'\t' read -r CHR START_BP END_BP ANCESTRY; do
    mkdir -p $LD_BLOCKS_DIR/$ANCESTRY
    PLINK_OUTPUT=$LD_BLOCKS_DIR/${ANCESTRY}/${CHR}_${START_BP}_${END_BP}
    RANGE_FILE=$LD_BLOCKS_DIR/range_file.tmp
    echo "$CHR $START_BP $END_BP ${ANCESTRY}_${CHR}_${START_BP}_${END_BP}" > $RANGE_FILE

    plink1.9 --bfile $THOUSAND_GENOMES/$ANCESTRY --chr $CHR --extract range $RANGE_FILE --r square spaces --out $PLINK_OUTPUT --keep-allele-order
    plink1.9 --bfile $THOUSAND_GENOMES/$ANCESTRY --chr $CHR --extract range $RANGE_FILE  --out $PLINK_OUTPUT  --keep-allele-order --freq
    plink1.9 --bfile $THOUSAND_GENOMES/$ANCESTRY --chr $CHR --extract range $RANGE_FILE  --out $PLINK_OUTPUT  --keep-allele-order --make-just-bim
  done
} < data/ld_regions.tsv
