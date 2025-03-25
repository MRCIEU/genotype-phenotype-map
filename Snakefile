import os
import numpy as np
import pandas as pd
import subprocess
import sys

os.chdir('pipeline_steps')
TEST_RUN = os.getenv('TEST_RUN')
DATA_DIR = os.getenv('DATA_DIR')
RESULTS_DIR = os.getenv('RESULTS_DIR')
TIMESTAMP = os.getenv('TIMESTAMP')

PIPELINE_METADATA = DATA_DIR + 'pipeline_metadata/'
STUDY_DIR = DATA_DIR + 'study/'
LD_BLOCK_DATA_DIR = DATA_DIR + 'ld_blocks/'

### INPUT DATA FILES
studies_to_process_file = PIPELINE_METADATA + 'studies_to_process.tsv'
studies_to_process = pd.read_csv(studies_to_process_file , sep='\t')

if len(studies_to_process) == 0:
    print('No studies to process, exiting.')
    sys.exit()

ld_blocks = pd.read_csv('data/ld_blocks.tsv', sep='\t')

relevant_ancestries = np.isin(ld_blocks['ancestry'], studies_to_process['ancestry'].unique())
ld_blocks = ld_blocks[relevant_ancestries]
ld_blocks = [f'{ld.ancestry}/{ld.chr}/{ld.start}-{ld.stop}' for i, ld in ld_blocks.iterrows()]

#TODO: if imputation goes well, we can remove the idea of simple and complex ld blocks, and just have ld blocks
complex_ld_blocks = []
simple_ld_blocks = [block for block in ld_blocks if block not in complex_ld_blocks]
# if TEST_RUN == 'test':
#     complex_ld_blocks = ['EUR/8/80453471-82816871']
# TODO: this is just removing badly behaving regions, this shouldn't be kept around later
complex_ld_blocks = []

extracted_studies = [s["extracted_location"] for i,s in studies_to_process.iterrows()]
extracted_study_pattern = '{study_location}extracted_snps.tsv'

standardisation_pattern = LD_BLOCK_DATA_DIR + '{simple_ld_block}/standardisation_complete'
complex_standardisation_pattern = LD_BLOCK_DATA_DIR + '{complex_ld_block}/complex_standardisation_complete'
imputation_pattern = LD_BLOCK_DATA_DIR + '{simple_ld_block}/imputation_complete'
complex_imputation_pattern = LD_BLOCK_DATA_DIR + '{complex_ld_block}/complex_imputation_complete'
finemapping_pattern = LD_BLOCK_DATA_DIR + '{simple_ld_block}/finemapping_complete'
complex_finemapping_pattern = LD_BLOCK_DATA_DIR + '{complex_ld_block}/complex_finemapping_complete'
coloc_pattern = LD_BLOCK_DATA_DIR + '{simple_ld_block}/coloc_complete'
complex_coloc_pattern = LD_BLOCK_DATA_DIR + '{complex_ld_block}/complex_coloc_complete'
compare_rare_pattern = LD_BLOCK_DATA_DIR + '{simple_ld_block}/compare_rare_complete'
complex_compare_rare_pattern = LD_BLOCK_DATA_DIR + '{complex_ld_block}/complex_compare_rare_complete'

### OUTPUT DATA FILES
studies_processed_file = RESULTS_DIR + 'studies_processed.tsv'
ld_blocks_to_process = f'{PIPELINE_METADATA}updated_ld_blocks_to_colocalise.tsv'
current_results_dir = f'{RESULTS_DIR}{TIMESTAMP}'

raw_coloc_results = f'{current_results_dir}/raw_coloc_results.tsv'
rare_results = f'{current_results_dir}/rare_results.tsv'
study_extractions = f'{current_results_dir}/study_extractions.tsv'
results_metadata = f'{current_results_dir}/results_metadata.tsv'
variant_annotations = f'{current_results_dir}/variant_annotations.tsv'
pipeline_summary_output = f'{current_results_dir}/pipeline_summary.html'
gpm_db_file = f'{RESULTS_DIR}gpm.db'
backup_done_file = '/tmp/backup_done'

rule all:
    input: expand(extracted_study_pattern, study_location=extracted_studies),
        ld_blocks_to_process,
        expand(imputation_pattern, simple_ld_block=simple_ld_blocks),
        expand(complex_imputation_pattern, complex_ld_block=complex_ld_blocks),
        expand(finemapping_pattern, simple_ld_block=simple_ld_blocks),
        expand(complex_finemapping_pattern, complex_ld_block=complex_ld_blocks),
        expand(coloc_pattern, simple_ld_block=simple_ld_blocks),
        expand(complex_coloc_pattern, complex_ld_block=complex_ld_blocks),
        expand(compare_rare_pattern, simple_ld_block=simple_ld_blocks),
        expand(complex_compare_rare_pattern, complex_ld_block=complex_ld_blocks),
        raw_coloc_results,
        rare_results,
        study_extractions,
        results_metadata
        # gpm_db_file
        # backup_done_file,
        # pipeline_summary_output

rule extract_regions_from_studies:
    params: lambda wildcards: list(filter(bool, wildcards.study_location.split("/")))[-1]
    output: extracted_study_pattern
    threads: 1
    retries: 1
    run:
        study = studies_to_process[studies_to_process.study_name == str(params)]
        if (len(study) != 1): raise ValueError(f'More than 1 study found for {str(params)}')
        study = study.iloc[0]

        if study.data_format == 'opengwas':
            command = f'Rscript extract_regions_from_opengwas.R \
                --extracted_study_location {study.extracted_location} \
                --extracted_output_file {output}'
        elif study.data_format == 'besd':
            command = f'Rscript extract_regions_from_besd.R \
                --extracted_study_location {study.extracted_location} \
                --extracted_output_file {output}'
        elif study.data_format == 'tsv' and study.variant_type != 'common':
            command = f'Rscript extract_regions_from_rare_tsv.R \
                --extracted_study_location {study.extracted_location} \
                --extracted_output_file {output}'
        else:
            raise ValueError(f'Cant ingest unknown data format: {study.data_format}')

        subprocess.run(command, shell=True)

rule organise_extracted_studies_into_ld_blocks:
    input: expand(extracted_study_pattern, study_location=extracted_studies)
    output: temporary(ld_blocks_to_process)
    threads: 1
    shell:
        """
        Rscript organise_extracted_regions_into_ld_blocks.R --output_file {output}
        """

def standardise_rule(standardisation_pattern, name):
    rule:
        name: f'{name}_standardise_per_ld_block'
        input: ld_blocks_to_process
        output: temporary(standardisation_pattern)
        retries: 1
        threads: 1
        params:
            ld_dir=lambda wildcards, output: os.path.dirname(output[0])
        run:
            ld_blocks = pd.read_csv(ld_blocks_to_process, sep='\t')
            skip_block = len(ld_blocks[ld_blocks.data_dir == params.ld_dir]) == 0

            if skip_block:
                command = f"mkdir -p $(dirname {output}) && touch {output}"
            else:
                ld_block = params.ld_dir.replace(LD_BLOCK_DATA_DIR, '')
                command = f"Rscript standardise_studies_in_ld_block.R \
                    --ld_block {ld_block} \
                    --completed_output_file {output}"
            subprocess.run(command, shell=True)


def impute_rule(standardisation_pattern, imputation_pattern, name):
    rule:
        name: f'{name}_impute_per_ld_block'
        input: standardisation_pattern
        output: temporary(imputation_pattern)
        retries: 1
        threads: 5
        params:
            ld_dir=lambda wildcards, output: os.path.dirname(output[0])
        run:
            ld_blocks = pd.read_csv(ld_blocks_to_process, sep='\t')
            skip_block = len(ld_blocks[ld_blocks.data_dir == params.ld_dir]) == 0

            if skip_block:
                command = f"mkdir -p $(dirname {output}) && touch {output}"
            else:
                ld_block = params.ld_dir.replace(LD_BLOCK_DATA_DIR, '')
                command = f"Rscript impute_studies_in_ld_block.R \
                    --ld_block {ld_block} \
                    --completed_output_file {output}"
            subprocess.run(command, shell=True)

def finemap_rule(imputation_pattern, finemaping_pattern, name):
    rule:
        name: f'{name}_finemap_per_ld_block'
        input: imputation_pattern 
        output: temporary(finemaping_pattern)
        retries: 1
        threads: 3
        params:
            ld_dir=lambda wildcards, output: os.path.dirname(output[0])
        run:
            ld_blocks = pd.read_csv(ld_blocks_to_process, sep='\t')
            skip_block = len(ld_blocks[ld_blocks.data_dir == params.ld_dir]) == 0

            if skip_block:
                command = f"mkdir -p $(dirname {output}) && touch {output}"
            else:
                ld_block = params.ld_dir.replace(LD_BLOCK_DATA_DIR, '')
                command = f"Rscript finemap_studies_in_ld_block.R \
                    --ld_block {ld_block} \
                    --completed_output_file {output} \
                    --complex_block {name == 'complex'}"
            subprocess.run(command, shell=True)

def coloc_rule(finemapping_pattern, coloc_pattern, name):
    rule:
        name: f'{name}_coloc_per_ld_block'
        retries: 2
        threads: 12
        # threads: 20
        input:
            finemap = finemapping_pattern
        output: temporary(coloc_pattern)
        params:
            ld_dir=lambda wildcards, output: os.path.dirname(output[0])
        run:
            ld_blocks = pd.read_csv(ld_blocks_to_process, sep='\t')
            skip_block = len(ld_blocks[ld_blocks.data_dir == params.ld_dir]) == 0

            if skip_block:
                command = f"mkdir -p $(dirname {output}) && touch {output}"
            else:
                ld_block = params.ld_dir.replace(LD_BLOCK_DATA_DIR, '')
                command = f"Rscript colocalise_studies_in_ld_block.R \
                    --ld_block {ld_block} \
                    --completed_output_file {output}"
            subprocess.run(command, shell=True)

def compare_rare_rule(standardisation_pattern, compare_rare_pattern, name):
    rule:
        name: f'{name}_compare_rare_per_ld_block'
        retries: 2
        threads: 2
        input:
            standardisation = standardisation_pattern 
        output: temporary(compare_rare_pattern)
        params:
            ld_dir=lambda wildcards, output: os.path.dirname(output[0])
        run:
            ld_blocks = pd.read_csv(ld_blocks_to_process, sep='\t')
            skip_block = len(ld_blocks[ld_blocks.data_dir == params.ld_dir]) == 0

            if skip_block:
                command = f"mkdir -p $(dirname {output}) && touch {output}"
            else:
                ld_block = params.ld_dir.replace(LD_BLOCK_DATA_DIR, '')
                command = f"Rscript compare_rare_studies_in_ld_block.R \
                    --ld_block {ld_block} \
                    --completed_output_file {output}"
            subprocess.run(command, shell=True)

standardise_rule(complex_standardisation_pattern, 'complex')
standardise_rule(standardisation_pattern, 'simple')

impute_rule(complex_standardisation_pattern, complex_imputation_pattern,'complex')
impute_rule(standardisation_pattern, imputation_pattern,'simple')

finemap_rule(complex_imputation_pattern, complex_finemapping_pattern, 'complex')
finemap_rule(imputation_pattern, finemapping_pattern, 'simple')

coloc_rule(complex_finemapping_pattern, complex_coloc_pattern, 'complex')
coloc_rule(finemapping_pattern, coloc_pattern, 'simple')

compare_rare_rule(complex_standardisation_pattern, complex_compare_rare_pattern, 'complex')
compare_rare_rule(standardisation_pattern, compare_rare_pattern, 'simple')

#TODO: tempoararily removing all post processing steps, for speed

rule compile_results:
    input: expand(coloc_pattern, simple_ld_block=simple_ld_blocks), expand(complex_coloc_pattern, complex_ld_block=complex_ld_blocks),
        expand(compare_rare_pattern, simple_ld_block=simple_ld_blocks), expand(complex_compare_rare_pattern, complex_ld_block=complex_ld_blocks)
    threads: 1
    output:
        raw_coloc_results = raw_coloc_results,
        rare_results = rare_results,
        study_extractions = study_extractions,
        results_metadata = results_metadata,
        variant_annotations = variant_annotations,
        # pipeline_summary = pipeline_summary_output
    shell:
        """
        mkdir -p $(dirname {output})
        Rscript compile_results.R \
            --studies_to_process {studies_to_process_file} \
            --studies_processed {studies_processed_file} \
            --study_extractions_file {output.study_extractions} \
            --rare_results_file {output.rare_results} \
            --coloc_results_file {output.raw_coloc_results} \
            --compiled_results_metadata_file {output.results_metadata}

#         rsync -Lavzh $RESULTS_DIR $BACKUP_DIR/results/ --exclude=".*"
#         """

# rule backup_data_dir:
#     input: raw_coloc_results, rare_results, study_extractions, results_metadata, variant_annotations
#     threads: 1
#     output: temporary(backup_done_file)
#     shell:
#         """
#         rsync -Lavzh --exclude='*cached*' $DATA_DIR/ld_blocks $BACKUP_DIR/data/
#         rsync -Lavzh --ignore-missing-args $DATA_DIR/study $BACKUP_DIR/data/
#         touch {output}
#         """

#rule create_results_db:
#    input: raw_coloc_results, rare_results, study_extractions, results_metadata, variant_annotations
#    threads: 1
#    output:
#        gpm_db = gpm_db_file,
#    shell:
#        """
#        Rscript create_db_from_results.R \
#            --results_dir {current_results_dir} \
#            --gpm_db_file {output.gpm_db}
#        """

onsuccess:
    print('Yay!  Please look here:')
    # print(pipeline_summary_output)
    print(raw_coloc_results)
    print(rare_results)
    print(study_extractions)
    print(results_metadata)

onerror:
    print(':(')
