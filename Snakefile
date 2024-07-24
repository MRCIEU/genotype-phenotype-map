import os
import numpy as np
import pandas as pd
import subprocess

os.chdir('pipeline_steps')
TEST_RUN = os.getenv('TEST_RUN')
DATA_DIR = os.getenv('DATA_DIR')
RESULTS_DIR = os.getenv('RESULTS_DIR')
TIMESTAMP = os.getenv('TIMESTAMP')

PIPELINE_METADATA = DATA_DIR + 'pipeline_metadata/'
STUDY_DIR = DATA_DIR + 'study/'
LD_BLOCK_DATA_DIR = DATA_DIR + 'ld_blocks/'
LD_BLOCK_MATRICES_DIR = DATA_DIR + 'ld_block_matrices/'
LD_BLOCK_RESULTS_DIR = RESULTS_DIR + 'ld_blocks/'

### INPUT DATA FILES
studies_to_process_file = PIPELINE_METADATA + 'studies_to_process.tsv'
studies_to_process = pd.read_csv(studies_to_process_file , sep='\t')
ld_regions = pd.read_csv('data/ld_regions.tsv', sep='\t')

relevant_ancestries = np.isin(ld_regions['ancestry'], studies_to_process['ancestry'].unique())
ld_regions = ld_regions[relevant_ancestries]
ld_blocks = [f'{ld.ancestry}/{ld.chr}/{ld.start}_{ld.stop}' for i, ld in ld_regions.iterrows()]

complex_ld_blocks = ['EUR/6/19207487_21684064',
                     'EUR/8/116096495_119685456',
                     'EUR/6/31571218_32682663',
                     'EUR/6/30798168_31571217',
                     ]


simple_ld_blocks = [block for block in ld_blocks if block not in complex_ld_blocks]
complex_ld_blocks = ['EUR/6/19207487_21684064',
                     'EUR/8/116096495_119685456',
                     ]
if TEST_RUN == 'test':
    complex_ld_blocks = []

extracted_studies = [s["extracted_location"] for i,s in studies_to_process.iterrows()]
extracted_study_pattern = '{study_location}extracted_snps.tsv'

imputation_pattern = LD_BLOCK_DATA_DIR + '{simple_ld_block}/imputation_complete'
complex_imputation_pattern = LD_BLOCK_DATA_DIR + '{complex_ld_block}/complex_imputation_complete'
finemapping_pattern = LD_BLOCK_DATA_DIR + '{simple_ld_block}/finemapping_complete'
complex_finemapping_pattern = LD_BLOCK_DATA_DIR + '{complex_ld_block}/complex_finemapping_complete'

ld_block_matrices = [f'{LD_BLOCK_MATRICES_DIR}{ld.ancestry}/{ld.chr}_{ld.start}_{ld.stop}' for i, ld in ld_regions.iterrows()]
ld_blocks_lookup = [f'{LD_BLOCK_DATA_DIR}{ld.ancestry}/{ld.chr}/{ld.start}_{ld.stop}' for i, ld in ld_regions.iterrows()]
ld_block_matrices = dict(zip(ld_blocks_lookup, ld_block_matrices))


### OUTPUT DATA FILES
coloc_pattern = LD_BLOCK_RESULTS_DIR + '{simple_ld_block}/coloc_complete'
complex_coloc_pattern = LD_BLOCK_RESULTS_DIR + '{complex_ld_block}/complex_coloc_complete'

studies_processed_file = RESULTS_DIR + 'studies_processed.tsv'
ld_blocks_to_process = f'{PIPELINE_METADATA}updated_ld_blocks_to_colocalise.tsv'
raw_coloc_results = f'{RESULTS_DIR}{TIMESTAMP}/raw_coloc_results.tsv'
coloc_results = f'{RESULTS_DIR}{TIMESTAMP}/coloc_results.tsv'
all_study_regions = f'{RESULTS_DIR}{TIMESTAMP}/all_study_regions.tsv'
mr_results = f'{RESULTS_DIR}{TIMESTAMP}/mr_results.tsv'
results_metadata = f'{RESULTS_DIR}{TIMESTAMP}/results_metadata.tsv'

rule all:
    input: expand(extracted_study_pattern, study_location=extracted_studies),
        ld_blocks_to_process,
        expand(imputation_pattern, simple_ld_block=simple_ld_blocks),
        expand(complex_imputation_pattern, complex_ld_block=complex_ld_blocks),
        expand(finemapping_pattern, simple_ld_block=simple_ld_blocks),
        expand(complex_finemapping_pattern, complex_ld_block=complex_ld_blocks),
        expand(coloc_pattern, simple_ld_block=simple_ld_blocks),
        expand(complex_coloc_pattern, complex_ld_block=complex_ld_blocks),
        raw_coloc_results,
        coloc_results,
        all_study_regions,
        results_metadata

rule extract_regions_from_studies:
    params: lambda wildcards: list(filter(bool, wildcards.study_location.split("/")))[-1]
    output: extracted_study_pattern
    threads: 1
    run:
        study = studies_to_process[studies_to_process.study_name == str(params)]
        if (len(study) != 1): raise ValueError(f'More than 1 study found for {str(params)}')
        study = study.iloc[0]

        if study.data_format == 'opengwas':
            command = f'./extract_regions_from_opengwas.sh \
                {study.study_name} \
                {study.ancestry} \
                {study.sample_size} \
                {study.category} \
                {study.study_location} \
                {study.extracted_location} \
                {study.p_value_threshold}'
        elif study.data_format == 'besd':
            command = f'Rscript extract_regions_from_besd.R \
                --extracted_study_location {study.extracted_location} \
                --extracted_output_file {output}'
        else:
            raise ValueError(f'Cant ingest unknown data format: {study.data_format}')

        subprocess.run(command, shell=True)

rule organise_extracted_studies_into_ld_regions:
    input: expand(extracted_study_pattern, study_location=extracted_studies)
    output: temporary(ld_blocks_to_process)
    threads: 1
    shell:
        """
        Rscript organise_extracted_regions_into_ld_regions.R --output_file {output}
        """

def impute_rule(defined_pattern, name):
    rule:
        name: f'{name}_impute_per_ld_block'
        input: ld_blocks_to_process
        output: temporary(defined_pattern)
        threads: 56 if name == 'complex' else 24
        priority: 1 if name == 'complex' else 0
        params:
            ld_dir=lambda wildcards, output: os.path.dirname(output[0])
        run:
            ld_block = params.ld_dir.replace(LD_BLOCK_DATA_DIR, '')
            ld_blocks = pd.read_csv(ld_blocks_to_process, sep='\t')
            skip_block = len(ld_blocks[ld_blocks.data_dir == params.ld_dir]) == 0

            if name == 'complex':
                #env_vars = "export LD_PRELOAD= && export OMP_NUM_THREADS=32 && export MKL_NUM_THREADS=32 && NUMEXPR_NUM_THREADS=32"
                env_vars = "export LD_PRELOAD="
            else:
                env_vars = "export LD_PRELOAD= && export OMP_NUM_THREADS=16 && export MKL_NUM_THREADS=16 && NUMEXPR_NUM_THREADS=16"

            if skip_block:
                command = f"mkdir -p $(dirname {output}) && touch {output}"
            else:
                command = f"{env_vars} && python3 standardise_and_impute_region.py \
                    --ld_block {ld_block} \
                    --completed_output_file {output}"
            subprocess.run(command, shell=True)

def finemap_rule(imputation_pattern, finemaping_pattern, name):
    rule:
        name: f'{name}_finemap_per_ld_block'
        input: imputation_pattern 
        output: temporary(finemaping_pattern)
        threads: 4
        params:
            ld_dir=lambda wildcards, output: os.path.dirname(output[0])
        run:
            ld_block = params.ld_dir.replace(LD_BLOCK_DATA_DIR, '')
            ld_blocks = pd.read_csv(ld_blocks_to_process, sep='\t')
            skip_block = len(ld_blocks[ld_blocks.data_dir == params.ld_dir]) == 0

            if skip_block:
                command = f"mkdir -p $(dirname {output}) && touch {output}"
            else:
                command = f"Rscript finemap_extracted_regions.R \
                    --ld_block {ld_block} \
                    --completed_output_file {output}"
            subprocess.run(command, shell=True)

def coloc_rule(finemapping_pattern, coloc_pattern, name):
    rule:
        name: f'{name}_coloc_per_ld_block'
        input:
            finemap = finemapping_pattern
        output: temporary(coloc_pattern)
        params:
            ld_dir=lambda wildcards, output: os.path.dirname(output[0])
        run:
            ld_block = params.ld_dir.replace(LD_BLOCK_RESULTS_DIR, '')
            ld_blocks = pd.read_csv(ld_blocks_to_process, sep='\t')
            skip_block = len(ld_blocks[ld_blocks.results_dir == params.ld_dir]) == 0

            if skip_block:
                command = f"mkdir -p $(dirname {output}) && touch {output}"
            else:
                command = f"Rscript colocalise_studies_per_ld_block.R \
                    --ld_block {ld_block} \
                    --coloc_result_file {output}"

            subprocess.run(command, shell=True)


impute_rule(complex_imputation_pattern,'complex')
impute_rule(imputation_pattern,'simple')

finemap_rule(complex_imputation_pattern, complex_finemapping_pattern, 'complex')
finemap_rule(imputation_pattern, finemapping_pattern, 'simple')

coloc_rule(complex_finemapping_pattern, complex_coloc_pattern, 'complex')
coloc_rule(finemapping_pattern, coloc_pattern, 'simple')

rule compile_results:
    input: expand(coloc_pattern, simple_ld_block=simple_ld_blocks), expand(complex_coloc_pattern, complex_ld_block=complex_ld_blocks)
    threads: 2
    output:
        coloc_results = coloc_results,
        raw_coloc_results = raw_coloc_results,
        all_study_regions = all_study_regions,
        results_metadata = results_metadata
    shell:
       """
       mkdir -p $(dirname {output})
       Rscript compile_results.R \
           --studies_to_process {studies_to_process_file} \
           --studies_processed {studies_processed_file} \
           --all_study_regions_file {output.all_study_regions} \
           --raw_coloc_results_file {output.raw_coloc_results} \
           --coloc_results_file {output.coloc_results} \
           --compiled_results_metadata_file {output.results_metadata}
       """

# rule perform_mr_analysis:
#     input:
#         all_study_regions: all_study_regions,
#         coloc_results: coloc_results,
#     output: mr_results
#     shell:
#         """
#         Rscript perform_mr_analysis.R --all_study_regions {input.all_study_regions} --coloc_results {input.coloc_results} --mr_result_file {output}
#         """

onsuccess:
    print('Yay!  Please look here:')
    print(coloc_results)
    print(all_study_regions)
    print(results_metadata)

onerror:
    print(':(')
