#rm -r /local-scratch/projects/genotype-phenotype-map/data/study/ukb-b-10003/*
#rm -r /local-scratch/projects/genotype-phenotype-map/prototype/data/study/*
#rm -r /local-scratch/projects/genotype-phenotype-map/prototype/data/ld_blocks/*

EXTRA_ARG=$1

Rscript identify_studies_to_process.R
export TIMESTAMP=$(date +%Y_%m_%d-%H_%M)
#snakemake -s common_variant.smk -c8 &> /tmp/snakemake.log

IMAGE=docker://andrewrrelmore/genotype_phenotype:latest
PIPELINE=common_cis_variants

snakemake -s common_variant.smk --profile ./ $EXTRA_ARG &> /tmp/snakemake.log

#apptainer run -B /local-scratch \
#              -B $(pwd)/$PIPELINE:/home/$PIPELINE \
#              -B /home/$(whoami) \
#              -B /projects \
#              --env-file $PIPELINE/.env --env TIMESTAMP=$TIMESTAMP\
#              --pwd /home/$PIPELINE/ \
#              $IMAGE $COMMAND
