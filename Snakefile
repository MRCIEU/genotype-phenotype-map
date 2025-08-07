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

extracted_studies = [s["extracted_location"] for i,s in studies_to_process.iterrows()]
extracted_study_pattern = '{study_location}extracted_snps.tsv'

standardisation_pattern = LD_BLOCK_DATA_DIR + '{ld_block}/standardisation_complete'
imputation_pattern = LD_BLOCK_DATA_DIR + '{ld_block}/imputation_complete'
finemapping_pattern = LD_BLOCK_DATA_DIR + '{ld_block}/finemapping_complete'
coloc_pattern = LD_BLOCK_DATA_DIR + '{ld_block}/coloc_complete'
compare_rare_pattern = LD_BLOCK_DATA_DIR + '{ld_block}/compare_rare_complete'

### OUTPUT DATA FILES
studies_processed_file = RESULTS_DIR + 'latest/studies_processed.tsv.gz'
traits_processed_file = RESULTS_DIR + 'latest/traits_processed.tsv.gz'
ld_blocks_to_process = f'{PIPELINE_METADATA}updated_ld_blocks_to_colocalise.tsv'
current_results_dir = RESULTS_DIR + 'current'
timestamped_results_dir = f'{RESULTS_DIR}{TIMESTAMP}'

coloc_pairwise_results = f'{current_results_dir}/coloc_pairwise_results.tsv.gz'
coloc_clustered_results = f'{current_results_dir}/coloc_clustered_results.tsv.gz'
rare_results = f'{current_results_dir}/rare_results.tsv.gz'
study_extractions = f'{current_results_dir}/study_extractions.tsv.gz'
new_studies_processed = f'{current_results_dir}/studies_processed.tsv.gz'
new_traits_processed = f'{current_results_dir}/traits_processed.tsv.gz'
pipeline_summary_output = f'{current_results_dir}/pipeline_summary.html'

studies_db_file = f'{current_results_dir}/studies.db'
associations_db_file = f'{current_results_dir}/associations.db'
coloc_pairs_db_file = f'{current_results_dir}/coloc_pairs.db'
ld_db_file = f'{current_results_dir}/ld.db'
gwas_upload_db_file = f'{current_results_dir}/gwas_upload.db'
create_dbs_done_file = f'{current_results_dir}/create_dbs_done'

backup_done_file = f'{current_results_dir}/backup_done'
sync_done_file = f'{current_results_dir}/sync_done'
svg_files_ready_file = f'{current_results_dir}/svg_files_ready'

rule all:
    input: expand(extracted_study_pattern, study_location=extracted_studies),
        ld_blocks_to_process,
        expand(imputation_pattern, ld_block=ld_blocks),
        expand(finemapping_pattern, ld_block=ld_blocks),
        expand(coloc_pattern, ld_block=ld_blocks),
        expand(compare_rare_pattern, ld_block=ld_blocks),
        coloc_pairwise_results,
        coloc_clustered_results,
        rare_results,
        study_extractions,
        new_studies_processed,
        new_traits_processed,
        studies_db_file,
        associations_db_file,
        coloc_pairs_db_file,
        ld_db_file,
        gwas_upload_db_file,
        create_dbs_done_file,
        backup_done_file,
        svg_files_ready_file
        # sync_done_file,
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

rule standardise_rule:
    name: f'standardise_per_ld_block'
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


rule:
    name: 'impute_per_ld_block'
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

rule finemap_rule:
    name: f'finemap_per_ld_block'
    input: imputation_pattern 
    output: temporary(finemapping_pattern)
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
                --completed_output_file {output}"
        subprocess.run(command, shell=True)

rule coloc_rule:
    name: f'coloc_per_ld_block'
    retries: 2
    threads: 2
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
            command = f"Rscript coloc_and_cluster_studies_in_ld_block.R \
                --ld_block {ld_block} \
                --completed_output_file {output}"
        subprocess.run(command, shell=True)

rule compare_rare_rule:
    name: f'compare_rare_per_ld_block'
    retries: 2
    threads: 3
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

rule compile_results:
    input: expand(coloc_pattern, ld_block=ld_blocks), expand(compare_rare_pattern, ld_block=ld_blocks)
    threads: 1
    output:
        coloc_pairwise_results = coloc_pairwise_results,
        coloc_clustered_results = coloc_clustered_results,
        rare_results = rare_results,
        study_extractions = study_extractions,
        new_studies_processed = new_studies_processed,
        new_traits_processed = new_traits_processed
    shell:
        """
        mkdir -p $(dirname {output})
        Rscript compile_results.R \
            --studies_to_process {studies_to_process_file} \
            --studies_processed {studies_processed_file} \
            --traits_processed {traits_processed_file} \
            --new_studies_processed_file {output.new_studies_processed} \
            --new_traits_processed_file {output.new_traits_processed} \
            --study_extractions_file {output.study_extractions} \
            --coloc_pairwise_results_file {output.coloc_pairwise_results} \
            --coloc_clustered_results_file {output.coloc_clustered_results} \
            --rare_results_file {output.rare_results}

         rsync -Lavzh $RESULTS_DIR $BACKUP_DIR/results/ --exclude=".*"
         """

# rule backup_data_dir:
#     input: coloc_pairwise_results, coloc_clustered_results, rare_results, study_extractions
#     threads: 1
#     output: temporary(backup_done_file)
#     shell:
#         """
#         rsync -Lavzh --exclude='*cached*' $DATA_DIR/ld_blocks $BACKUP_DIR/data/
#         rsync -Lavzh --ignore-missing-args $DATA_DIR/study $BACKUP_DIR/data/
#         touch {output}
#         """

rule create_results_db:
   input: coloc_pairwise_results, coloc_clustered_results, rare_results, study_extractions, new_studies_processed, new_traits_processed
   threads: 1
   output:
       studies_db_file = studies_db_file,
       associations_db_file = associations_db_file,
       coloc_pairs_db_file = coloc_pairs_db_file,
       ld_db_file = ld_db_file,
       gwas_upload_db_file = gwas_upload_db_file,
       create_dbs_done_file = temporary(create_dbs_done_file)
   shell:
       """
       Rscript create_db_from_results.R \
           --studies_db_file {studies_db_file} \
           --associations_db_file {associations_db_file} \
           --coloc_pairs_db_file {coloc_pairs_db_file} \
           --ld_db_file {ld_db_file} \
           --gwas_upload_db_file {gwas_upload_db_file} \
           --completed_output_file {output.create_dbs_done_file}
       """

rule prepare_svg_files_for_use:
    input: studies_db_file, create_dbs_done_file
    threads: 1
    output: temporary(svg_files_ready_file)
    shell:
        """
        R -e "source('constants.R'); source('svg_helpers.R'); prepare_svg_files_for_use()"
        touch {output}
        """

rule copy_results_to_timestamped_dir:
    input: svg_files_ready_file
    threads: 1
    output: timestamped_results_dir
    shell:
        """
        mkdir -p {timestamped_results_dir}
        cp -r {current_results_dir} {timestamped_results_dir}
        """

onsuccess:
    print('Yay!  Please look here:')
    print(coloc_pairwise_results)
    print(coloc_clustered_results)
    print(rare_results)
    print(study_extractions)

onerror:
    print(':(')
