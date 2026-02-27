# Date: 18-02-2025
# Author: A.Hanson
# Access Genebass matrix table format summary statistics and extract rare variant IDs to pass to VEP for variant annotation

import click
import pandas as pd
import hail as hl
import os
import sys

# Initialize Hail with a specified number of cores
hl.init(spark_conf={'spark.master': 'local[10]'})

os.environ['PYSPARK_SUBMIT_ARGS'] = '--driver-memory 10g --executor-memory 10g pyspark-shell'
DATA_DIR = "/local-scratch/data/hg38/genebass/"
OUT_DIR = "/local-scratch/data/ukb-seq/downloads/genebass/"

max_maf = 0.01
min_mac = 5

# Read in hail matix table
mt = hl.read_matrix_table(os.path.join(DATA_DIR, "variant_results.mt"))

# Filter variants
mt_snps = mt.filter_rows(((mt.call_stats.AC >= min_mac) & (mt.call_stats.AF <= max_maf))|
                         ((mt.call_stats.AN - mt.call_stats.AC >= min_mac) & (1 - mt.call_stats.AF <= max_maf)))

outname_snps = "genebass_rarevariant_info_hg38.tsv"

# Write out SNP information
mt_snps.rows().export(os.path.join(OUT_DIR, outname_snps))