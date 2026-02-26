#!/bin/bash

#SBATCH --partition=mrcieu
#SBATCH --job-name=vep_rare
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=15:00:00
#SBATCH --account=SMED001801
#SBATCH --output="veprare_%A.log"

# Date: 18-02-25
# Author: A.L.Hanson
# Purpose: Generate variant annotation file for all rare variants (SNPs and indels) at MAF<0.01 and MAC > 5 in Genebass

# Run on BC4 with local config

source ~/config.sh

module load apps/vep/112
module load languages/R/4.4.1

ANNOT_DIR=${DATA_DIR}/variant_annotation

echo "Running VEP with check for reference mismatch..."
echo "Start:" `date`

${VEP}/vep -i ${ANNOT_DIR}/genebass_rarevariant_info_hg38_forVEP_posupdate.txt \
--cache \
--dir_cache ${CACHE} \
--assembly GRCh38 \
--af_1kg \
--sift b \
--polyphen b \
--regulatory \
--show_ref_allele \
--symbol \
--protein \
--canonical \
--mane \
--biotype \
--check_ref \
--skipped_variants_file ${ANNOT_DIR}/vep_rare_skipped.txt \
--pick \
--no_stats \
--force_overwrite \
-v \
-o ${ANNOT_DIR}/vep_rare_variantannotations_hg38.txt

echo "End:" `date`