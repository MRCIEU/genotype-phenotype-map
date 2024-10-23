#!/bin/bash

keys='/local-scratch/data/ukb-seq/downloads/genebass/study_keys.txt'

while read -r k1 k2 k3 k4 k5; do
    python3 02_pullfromgenebass.py ${k1} ${k2} ${k3} ${k4} ${k5}
    rm *.log
done < ${keys}
