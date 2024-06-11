#!/bin/bash

USER=$(whoami)
export $(xargs <.env)
LOGS=$DATA_DIR/pipeline_metadata/logs
mkdir -p $LOGS
export PIPELINE=common_cis_variants

./run_pipeline_step.sh 'Rscript 00_format_ld_matrices.R' &> $LOGS/00_format_ld_matrices.log || exit 1
exit 0

#STEP 1:
#./run_pipeline_step.sh 'Rscript 1_gwas_subset_list.R' &> $LOGS/1_gwas_subset_list.log || exit 1
#./run_pipeline_step.sh './2_extract_regions_from_vcf.sh' &> $LOGS/2_extract_regions_from_vcf.log || exit 1
./run_pipeline_step.sh 'Rscript 4_organise_extracted_regions_into_ld_regions.R' &> $LOGS/4_organise_extracted_regions_into_ld_regions.log || exit 1
echo $?
if [[ $? -ne 0 ]]; then
  exit 1
fi
exit 0

for chr in $(seq 1 22); do
  ./run_pipeline_step.sh "Rscript 5_colocalise_studies_per_ld_block.R --chr $chr" &> $LOGS/5_colocalise_studies_per_ld_block_$chr.log || exit 1
done

exit 0
#etc...

#If we need to run it in a slurm job in the future...
#STEP1_ID=$(sbatch --parsable \
#                   --account=$ACCOUNT_ID --mem=1G --time=04:00:00 \
#                   --nodes=1 --ntasks-per-node=1 --cpus-per-task=1 \
#                   --output=$LOGS/1_gwas_subset_list.log \
#                   --export=COMMAND="Rscript 1_gwas_subset_list.R" \
#                   run_pipeline_step.sh
#)
