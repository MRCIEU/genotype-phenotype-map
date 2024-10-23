# Date: 23-09-2024
# Author: A.Hanson
# Access Genebass matrix table format summary statistics for UKB WES with hail and export top hits per study

import pandas as pd
import hail as hl
import os
import sys

trait_type = sys.argv[1]
phenocode = sys.argv[2]
pheno_sex = sys.argv[3]
coding = sys.argv[4]
modifier = sys.argv[5]

os.environ['PYSPARK_SUBMIT_ARGS'] = '--driver-memory 10g --executor-memory 10g pyspark-shell'

DATA_DIR = "/local-scratch/data/genebass/"
SEQDATA_DIR = "/local-scratch/data/ukb-seq/"

max_pval = 0.1
max_maf = 0.01
min_mac = 10

# Read in hail matix table
mt = hl.read_matrix_table(os.path.join(DATA_DIR, "variant_results.mt"))

# Write out study metadata
#col_df = mt.cols().to_pandas()
#col_df['file_name'] = "genebass_ukbwes_" + "p" + col_df['phenocode'] + "_" + col_df['coding'] + "_rare_filtered.tsv"
#col_df.to_csv(os.path.join(SEQDATA_DIR, "downloads/genebass/genebass_wes_studymetadata.tsv"), sep = '\t', index = False)

# E.g to extract column key for first study
# col_key = mt.cols().take(1)[0]

print("Filtering with values:")
print("trait_type:", trait_type)
print("phenocode:", phenocode)
print("pheno_sex:", pheno_sex)
print("coding:", coding)
print("modifier:", modifier)

# Filter the matrix table by the column key
mt_study = mt.filter_cols(
    (mt.trait_type == trait_type) &
    (mt.phenocode == phenocode) &
    (mt.pheno_sex == pheno_sex) &
    (mt.coding == coding) &
    (mt.modifier == modifier))

print("Check dimension:", mt_study.count())

desc =  mt_study.cols().select('description').collect()
print("pheno:", desc[0].description)

outname_rare = "".join(["genebass_ukbwes_","p",col_key['phenocode'], "_", col_key['coding'], "_rare_filtered.tsv"])
#outname_common = "".join(["genebass_ukbwes_","p",col_key['phenocode'], "_", col_key['coding'], "_common_filtered.tsv"])

# Return entries (summary statistics)
study_stats = mt_study.entries()

# Filter by p-value, MAF and MAC
print("Filtering...")

# Filter rare (p<0.1, MAF<0.01, MAC>10)
study_stats_rare = study_stats.filter(
    (study_stats.Pvalue <= max_pval) & 
    (((study_stats.AF <= max_maf) & (study_stats.AC >= min_mac)) |
        ((1 - study_stats.AF <= max_maf) & (study_stats.call_stats.AN - study_stats.call_stats.AC >= min_mac))))

study_stats_rare = study_stats_rare.drop('saige_version','description_more','category')

# Write out rare summary statistics
study_stats_rare.export(os.path.join(SEQDATA_DIR, "downloads/genebass/rare", outname_rare))

# Filter common (MAF>0.01) - no p-value filter to match common variant GWAS input
## study_stats_common = study_stats.filter(
#     (study_stats.AF > max_maf) & 
#     (study_stats.AF < 1 - max_maf))

# study_stats_common = study_stats_common.drop('saige_version','description_more','category')

# # Write out common summary statistics
# study_stats_common.export(os.path.join(SEQDATA_DIR, "downloads/genebass/common", outname_common))
