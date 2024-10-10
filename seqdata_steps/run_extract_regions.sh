#!/bin/bash

# To run e.g:
# ./run_extract_regions.sh bm

if [ $# -eq 0 ]; then
	echo "Please provide input study data format e.g (bm, gb, az or hd)"
	exit 1
fi

input=$1

echo "Extracting top hits and ld regions from ${input} summary statistics"

# Launch docker container
#apptainer run -B /local-scratch -B /projects  -B /home/${whoami} docker://andrewrrelmore/genotype_phenotype:latest /bin/bash

# Metadata files
backman="/local-scratch/data/ukb-seq/downloads/backmanexwas/ukb-wes-bm-studies.tsv"
azphewas="/local-scratch/data/ukb-seq/downloads/azexwas/ukb-wes-az-studies.tsv"
genebass=""
halldorsson=""

if [ ${input} == "bm" ]; then
	
	infile=${backman}
	(echo Reading study metadata from: ${infile}	
	
	for file in $(awk 'NR != 1 {print $5}' ${infile} | xargs -n 1 basename); do
    	Rscript extract_regions_from_backman.R --extracted_study_file ${file}
	done) 2>&1 | tee $(dirname ${backman})/studies_extract.log

elif [ ${input} == "az" ]; then
	
	infile=${azphewas}
	(echo Reading study metadata from: ${infile}	
	
	for file in $(awk 'NR != 1 {print $5}' ${infile} | xargs -n 1 basename); do
    	Rscript extract_regions_from_azphewas.R --extracted_study_file ${file}
	done) 2>&1 | tee $(dirname ${azphewas})/studies_extract.log

elif [ ${input} == "gb" ]; then
    infile=${genebass}
elif [ ${input} == "hd" ]; then
    infile=${halldorsson}
else
	echo "Input data format not recognised"
	exit 1
fi

exit 0
