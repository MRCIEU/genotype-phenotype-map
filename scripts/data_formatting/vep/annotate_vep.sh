#!/bin/bash

#SBATCH --partition=mrcieu
#SBATCH --job-name=vep
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=10:00:00
#SBATCH --account=SMED001801
#SBATCH --output="vep_%A.log"

# Date: 18-11-24
# Author: A.L.Hanson
# Purpose: Generate variant annotation file for all variants (SNPs and indels) used in the hg38 EUR LD reference panel

# Run on BC4 with local config

source ~/config.sh

module load apps/vep/112
module load languages/R/4.4.1

LD_DIR=${DATA_DIR}/ld_reference_panel_hg38/EUR
ANNOT_DIR=${DATA_DIR}/variant_annotation

# Extract variants used in LD reference

echo "Extracting required variant IDs..."
awk '{print $2}' ${LD_DIR}/full.bim > ${ANNOT_DIR}/hg38_variantset.txt
grep -wf ${ANNOT_DIR}/hg38_variantset.txt ${ANNOT_DIR}/variant_info_hg38_forVEP.txt > ${ANNOT_DIR}/variant_info_hg38_forVEP_runset.txt
 
echo "Correcting position for insertions and deletions..."
# Correct position for insertion and deletions
Rscript update_indelpositions.R ${ANNOT_DIR}/variant_info_hg38_forVEP_runset.txt

# echo "Splitting by chromosome..."
# Split by chromosome
# awk -v DIR=${ANNOT_DIR} '{print > DIR"/tmp_chr"$1"_input.txt"}' ${ANNOT_DIR}/variant_info_hg38_forVEP_runset_posupdate.txt

echo "Running VEP with check for reference mismatch..."
echo "Start:" `date`

${VEP}/vep -i ${ANNOT_DIR}/variant_info_hg38_forVEP_runset_posupdate.txt \
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
--skipped_variants_file ${ANNOT_DIR}/vep_skipped.txt \
--pick \
--no_stats \
--force_overwrite \
-v \
-o ${ANNOT_DIR}/vep_variantannotations_hg38.txt

echo "End:" `date`
