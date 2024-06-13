#!/bin/bash
set -e

LD_BLOCKS_DIR=$DATA_DIR/ld_block_matrices

while IFS= read -r LD_REGION; do
  CHR=$(echo $LD_REGION | awk '{print $1}')
  START_BP=$(echo $LD_REGION | awk '{print $2}')
  END_BP=$(echo $LD_REGION | awk '{print $3}')
  ANCESTRY=$(echo $LD_REGION | awk '{print $4}')

  mkdir -p $LD_BLOCKS_DIR/$ANCESTRY
  PLINK_OUTPUT=$LD_BLOCKS_DIR/${ANCESTRY}/${CHR}_${START_BP}_${END_BP}
  RANGE_FILE=$LD_BLOCKS_DIR/range_file.tmp
  echo "$CHR $START_BP $END_BP ${ANCESTRY}_${CHR}_${START_BP}_${END_BP}" > $RANGE_FILE

  plink1.9 --bfile $THOUSAND_GENOMES/$ANCESTRY --chr $CHR --extract range $RANGE_FILE --r square spaces --out $PLINK_OUTPUT --write-snplist --keep-allele-order

done < ../common_cis_variants/data/ld_regions.tsv
