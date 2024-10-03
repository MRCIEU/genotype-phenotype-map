#/bin/bash/!

# Launch docker container
apptainer run -B /local-scratch -B /projects  -B /home/$(whoami) docker://andrewrrelmore/genotype_phenotype:latest /bin/bash

# Metadata files
backman="/local-scratch/data/ukb-seq/downloads/backmanexwas/ukb-wes-bm-studies.tsv"

for file in $(awk 'NR != 1 {print $5}' ${backman} | xargs -n 1 basename | head -n 6); do
	Rscript extract_regions_from_backman.R --extracted_study_file ${file}
done
