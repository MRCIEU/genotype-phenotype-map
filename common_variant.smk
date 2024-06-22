import datetime
import os
import pandas as pd
import subprocess

DATA_DIR = os.getenv('DATA_DIR')
RESULTS_DIR = os.getenv('RESULTS_DIR')
PIPELINE_METADATA = DATA_DIR + 'pipeline_metadata/'
STUDY_DIR = DATA_DIR + 'study/'

onstart:
    print("##### Genotype-Phenotype Map Pipeline #####")

def clean_files_pipeline_metadata_for_processing():
    return


### INPUT DATA
subprocess.call(['Rscript', 'calculate_gwases_to_process.R'])
study_data = pd.read_csv(PIPELINE_METADATA + 'gwases_to_process.tsv', sep='\t')
studies = []

study_file_pattern = STUDY_DIR + '{study}/extraction_metadata.json'


### OUTPUT DATA
time = datetime.datetime.now()
final_report = f'{RESULTS_DIR}/report_{time:%Y_%m-%d-%H_%M}.tsv'

rule all:
    input: expand(study_file_pattern, prefix=[s for s in studies])


rule extract_regions_from_full_studies:
    input: list()
    output: list()
    shell:
        """
        Rscript extract_regions_from_opengwas.sh STUDY_DIR EXTRACTION_DIR ANCESTRY SAMPLE_SIZE P_VALUE DATA_TYPE
        """

rule organise_extracted_regions_into_ld_regions:
    input: list()
    output: list()
    shell:
        """
        Rscript 3_organise_extracted_regions_into_ld_regions.R
        """

rule impute_by_ld_region:
    input: list()
    output: list()
    shell:
        """
        export LD_PRELOAD=
        python3 4_impute_study.py
        """

rule finemap:
    input: list()
    output: list()
    shell:
        """
        Rscript 5_finemap_extracted_regions.R
        """

rule colocalise_per_ld_region:
    input: list()
    output: list()
    shell:
        """
        Rscript 6_colocalise_studies_per_ld_block.R
        """

rule mr_on_coloc_results:
    input: list()
    output: list()
    shell:
        """
        Rscript 7_perform_mr_analysis.R
        """

rule collate_info_on_all_studies:
    input: list()
    output: final_report
    shell:
        """
        Rscript 8_report_on_results_and_studies.R
        """

onsuccess:
    print('yay')

onerror:
    print('nooooo')
