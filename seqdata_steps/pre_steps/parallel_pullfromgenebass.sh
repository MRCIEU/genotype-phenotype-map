#!/bin/bash
set -e

keys='/local-scratch/data/ukb-seq/downloads/genebass/study_keys.txt'

while read -r k1 k2 k3 k4 k5; do
    echo "$k1 $k2 $k3 $k4 $k5"
    python3 02_pullfromgenebass.py --trait_type ${k1} --phenocode ${k2} --pheno_sex ${k3} --coding ${k4} --modifier ${k5}
    rm *.log
done < ${keys}
