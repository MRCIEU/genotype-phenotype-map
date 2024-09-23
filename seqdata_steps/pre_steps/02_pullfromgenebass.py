# Date: 23-09-2024
# Author: A.Hanson
# Access Genebass matrix table format summary statistics for UKB WES with hail and export top hits per study

import pandas as pd
import hail as hl
import os

DATA_DIR = "/local-scratch/data/genebass/"
SEQDATA_DIR = "/local-scratch/data/ukb-seq/"

max_pval = 0.1
max_maf = 0.01
min_mac = 10

mt = hl.read_matrix_table(os.path.join(DATA_DIR, "variant_results.mt"))

# Write out study metadata
col_df = mt.cols().to_pandas()
col_df.to_csv(os.path.join(SEQDATA_DIR, "downloads/genebass/genebass_wes_studymetadata.tsv"), sep = '\t', index = False)

#n_studies = mt.cols().count()
n_studies = 2

for study in range(1, n_studies + 1):

    print("Study number:", study)

    # Extract column key
    col_key = mt.cols().take(study)[study - 1]

    # Filter the matrix table by the column key
    mt_study = mt.filter_cols(
        (mt.trait_type == col_key['trait_type']) &
        (mt.phenocode == col_key['phenocode']) &
        (mt.pheno_sex == col_key['pheno_sex']) &
        (mt.coding == col_key['coding']) &
        (mt.modifier == col_key['modifier']))

    print("pheno:", col_key['description'])

    outname_rare = "".join(["genebass_ukbwes_","p",col_key['phenocode'], "_", col_key['coding'], "_rare_filtered.tsv"])
    outcome_common = "".join(["genebass_ukbwes_","p",col_key['phenocode'], "_", col_key['coding'], "_common_filtered.tsv"])

    # Return entries (summary statistics)
    study_stats = mt_study.entries()

    print("Filtering...")

    # Filter rare (p<0.1, MAF<0.01, MAC>10)
    study_stats_rare = study_stats.filter(
        (study_stats.Pvalue <= max_pval) & 
        (((study_stats.AF <= max_maf) & (study_stats.AC >= min_mac)) |
            ((1 - study_stats.AF <= max_maf) & (study_stats.call_stats.AN - study_stats.call_stats.AC >= min_mac))))

    # Write out rare summary statistics
    rare_df = study_stats_rare.to_pandas()
    rare_df.to_csv(os.path.join(SEQDATA_DIR, "downloads/genebass/rare", outname_rare), sep = '\t', index = False)

    # Filter common (p<0.1, MAF>0.01)
    study_stats_common = study_stats.filter(
        (study_stats.Pvalue <= max_pval) &
        (study_stats.AF > max_maf) & 
        (study_stats.AF < 1 - max_maf))

    # Write out common summary statistics
    commonm_df = study_stats_common.to_pandas()
    common_df.to_csv(os.path.join(SEQDATA_DIR, "downloads/genebass/common", outname_common), sep = '\t', index = False)
