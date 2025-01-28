#!/bin/bash

seqstudies_dir="/local-scratch/projects/genotype-phenotype-map/data/study_sequencing"
chr=7
bp_min=151346135
bp_max=153300525

find ${seqstudies_dir} -type f -name "EUR_${chr}_*.tsv.gz" | while read file; do
    num=$(echo "$file" | grep -oE "EUR_${chr}_[0-9]+\.tsv\.gz" | sed -E "s/EUR_${chr}_([0-9]+)\.tsv\.gz/\1/")
    if [ "$num" -ge ${bp_min} ] && [ "$num" -le ${bp_max} ]; then
        echo "$file"
    fi
done > data/example_studiesinblock_${chr}_${bp_min}-${bp_max}.txt