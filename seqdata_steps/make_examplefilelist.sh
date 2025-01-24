#!/bin/bash

seqstudies_dir="/local-scratch/projects/genotype-phenotype-map/data/study_sequencing"
chr=1
bp_min=95684841
bp_max=97259500

find ${seqstudies_dir} -type f -name "EUR_${chr}_*.tsv.gz" | while read file; do
    num=$(echo "$file" | grep -oE 'EUR_1_[0-9]+\.tsv\.gz' | sed -E 's/EUR_1_([0-9]+)\.tsv\.gz/\1/')
    if [ "$num" -ge ${bp_min} ] && [ "$num" -le ${bp_max} ]; then
        echo "$file"
    fi
done > example_studiesinblock.txt