#!/bin/bash

# Rare variant information extracted from genebass hail matrix table for variants with MAF<0.01 and MAC>5
# See: genotype-phenotype-map/pre_steps/vep_annotation/pull_rarevarids_fromgenebass.py

source ~/config.sh

varlist=${DATA_DIR}/variant_annotation/genebass_rarevariant_info_hg38.tsv
outfile=${DATA_DIR}/variant_annotation/genebass_rarevariant_info_hg38_forVEP.txt

# Format as VEP input:
# chr start end REF/ALT strand snpname

awk -F'\t' 'NR > 1 {
    split($1, pos, ":");
    split($3, id, "_");
    
    chr = pos[1];
    gsub("chr","",chr);

    if(chr >=1 && chr <= 22){
        print chr, pos[2], pos[2], id[2], "+"
    }    
}' $varlist > $outfile

# Alphabetise IDs
Rscript update_varid.R $outfile

# Correct positions of indels
Rscript update_indelpositions.R $outfile