import datetime
import os
import pandas as pd
from pathlib import Path
import subprocess
import re

#TODO: temporary files output: temporary()
# https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#protected-and-temporary-files

DATA_DIR = os.getenv('DATA_DIR')
RESULTS_DIR = os.getenv('RESULTS_DIR')
TIMESTAMP = os.getenv('TIMESTAMP')
PIPELINE_METADATA = DATA_DIR + 'pipeline_metadata/'
STUDY_DIR = DATA_DIR + 'study/'
LD_BLOCK_DATA_DIR = DATA_DIR + 'ld_blocks/'

onstart:
    print("##### Genotype-Phenotype Map Pipeline #####")

def clean_files_pipeline_metadata_for_processing():
    return

### INPUT DATA
process_file = PIPELINE_METADATA + 'studies_to_process.tsv'
studies_to_process = pd.read_csv(process_file, sep='\t')
ld_regions = pd.read_csv('data/ld_regions.tsv', sep='\t')

extracted_studies = [s["extracted_location"] for i,s in studies_to_process.iterrows()]
extracted_study_pattern = '{study_location}/extracted_snps.tsv'

imputation_pattern = LD_BLOCK_DATA_DIR + '{ld_block}/imputation_complete'
finemapping_pattern = LD_BLOCK_DATA_DIR + '{ld_block}/finemapping_complete'
ld_blocks = [f'{ld["pop"]}/{ld.chr}/{ld.start}_{ld.stop}' for i,ld in ld_regions.iterrows()]
ld_blocks = ['EUR/6/31571218_32682663']

### OUTPUT DATA
ld_blocks_to_process = f'{PIPELINE_METADATA}updated_ld_blocks_to_colocalise.tsv'
final_report = f'{RESULTS_DIR}report_{TIMESTAMP}.tsv'

rule all:
    input: expand(imputation_pattern, ld_block=ld_blocks),
        expand(finemapping_pattern, ld_block=ld_blocks),
        expand(extracted_study_pattern, study_location=extracted_studies),
        ld_blocks_to_process

rule extract_regions_from_studies:
    params: lambda wildcards: list(filter(bool, wildcards.study_location.split("/")))[-1]
    output: extracted_study_pattern
    run:
        print('extracting')
        study = studies_to_process[studies_to_process.study_name == params].values.flatten().tolist()
        command = ['./extract_regions_from_opengwas.sh'] + study[2:]
        command = [str(c) for c in command]
        subprocess.run(command)

rule organise_extracted_regions_into_ld_regions:
    input: expand(extracted_study_pattern, study_location=extracted_studies)
    output: temporary(ld_blocks_to_process)
    shell:
        """
        Rscript organise_extracted_regions_into_ld_regions.R --output_file {output}
        """

rule impute_by_ld_region:
    input: ld_blocks_to_process
    output: temporary(imputation_pattern)
    params:
        ld_dir=lambda wildcards, output: os.path.dirname(output[0])
    run:
        ld_blocks = pd.read_csv(ld_blocks_to_process, sep='\t')
        ld_block = ld_blocks[ld_blocks.data_dir == params.ld_dir]

        if len(ld_block) == 0:
            command = f"mkdir -p $(dirname {output}) && touch {output}"
        else:
            print(ld_block.iloc[0]['region_prefix'])
            command = f"export LD_PRELOAD= && python3 impute_region.py \
                --ld_region_prefix {ld_block.iloc[0]['region_prefix']} \
                --ld_block_dir {ld_block.iloc[0]['data_dir']}"
        subprocess.run(command, shell=True)


rule finemap_by_ld_region:
    input: ld_blocks_to_process, expand(imputation_pattern, ld_block=ld_blocks)
    output: temporary(finemapping_pattern)
    params:
        ld_dir=lambda wildcards, output: os.path.dirname(output[0])
    run:
        ld_blocks = pd.read_csv(ld_blocks_to_process, sep='\t')
        ld_block = ld_blocks[ld_blocks.data_dir == params.ld_dir]

        if len(ld_block) == 0:
            command = f"mkdir -p $(dirname {output}) && touch {output}"
        else:
            command = f"Rscript finemap_extracted_regions.R \
                --ld_region_prefix {ld_block.iloc[0]['region_prefix']} \
                --ld_block_dir {ld_block.iloc[0]['data_dir']}"
        subprocess.run(command, shell=True)

# rule colocalise_per_ld_region:
#     input: ld_blocks_to_process, expand(finemapping_pattern, ld_block=ld_blocks)
#     output: temporary(finemapping_pattern)
#     shell:
#         """
#         Rscript colocalise_studies_per_ld_block.R
#         """
#
# rule mr_on_coloc_results:
#     input: list()
#     output: list()
#     shell:
#         """
#         Rscript 7_perform_mr_analysis.R
#         """

#rule collate_info_on_all_studies:
#    input: extracted_study_pattern
#    output: final_report
#    shell:
#        """
#        Rscript 8_report_on_results_and_studies.R
#        """

onsuccess:
    print('yay')

onerror:
    print('nooooo')
