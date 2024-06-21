import pandas as pd
import subprocess

onstart:
    print("##### Genotype-Phenotype Map Pipeline #####")

def clean_files_pipeline_metadata_for_processing():
    return

subprocess.call (['Rscript', 'calculate_gwases_to_process.R'])
study_data = pd.read_csv(PIPELINE_METADATA + 'gwases_to_process.tsv', sep='\t')
studies = []

study_file_pattern = STUDY_DIR + '{study}/extraction_metadata.json'

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

onsuccess:
    print('yay')

onerror:
    print('nooooo')
