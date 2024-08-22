#!/bin/bash
set -e

THOUSAND_GENOMES=$DATA_DIR/1000genomes/
THOUSAND_GENOMES_FLIPPED=$DATA_DIR/1000genomes/flipped/
ANCESTRIES='AFR AMR EUR EAS SAS'

for ancestry in $ANCESTRIES; do
  PLINK_OUTPUT=$THOUSAND_GENOMES_FLIPPED/${ancestry}
  SNP_LIST=$THOUSAND_GENOMES_FLIPPED/${ancestry}.bim_toflip
  plink1.9 --bfile $THOUSAND_GENOMES/$ancestry --flip $SNP_LIST --out $PLINK_OUTPUT --keep-allele-order --make-bed
done
