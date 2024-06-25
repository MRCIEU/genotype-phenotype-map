rm /local-scratch/projects/genotype-phenotype-map/data/study/ukb-b-10003/extract*
Rscript identify_studies_to_process.R
export TIMESTAMP=$(date +%Y_%m_%d-%H_%M)
snakemake -s common_variant.smk -c200 &> /tmp/snakemake.log
