# genotype-phenotype-map

Scoping document: https://uob.sharepoint.com/:w:/r/teams/grp-brs-tds/Shared%20Documents/The%20Human%20Genotype-Phenotype%20Map/Kickoff%20Document.docx?d=w3c2c1890b3bf4ec6a8e3883780b05c94&csf=1&web=1&e=KFn1PY


## How to use:

1. Step 1. Clone the repo on to ieu-p1
    * `git clone git@github.com:MRCIEU/genotype-phenotype-map.git && cd genotype-phenotype-map`
2. Populate .env file
    * You will most likely want to use the default set in .env_example
    * `cp .env_example .env`
3. Add info in pipeline_steps/data/study_list.csv
    * This is the file that specifies which studies will be ingested by the pipeline, if you want to add a study or set of studies,
add a row to this csv indicating where they are on disk, and what script to ingest them with
4. `./run_pipeline.sh`
    * This will first calculate the studies that haven't been processed, then start the snakemake pipeline

Notes: Snakemake slows down if you try to run the pipeline with too many studies at once.  Please keep it to less than 200,000 studies.

## Concepts

* **Study:** Any individual summary statistic study to ingest.  A GWAS is a study (ex. BMI), a specific gene-tissue-expression is a study (ex. GTEx Whole Blood WASH7P)
* **LD Block**: The human genome has been chunked into regions which have minimal LD overlap, as calculated by [ldetect](https://github.com/jmacdon/LDblocks_GRCh38)
* **Docker / Apptainer**: All of the needed command line tools, R and python packages are installed in [a docker image](https://hub.docker.com/repository/docker/mrcieu/genotype-phenotype-map): 
* **Pipeline**: This is used to run the series of steps needed to ingest data and perform the calculations.  The pipeline uses [Snakemake version 7](https://snakemake.readthedocs.io/en/v7.25.4/index.html)

## How to contribute 

This has only been run on `ieu-p1.epi.bris.ac.uk`, please get access to there.

To make changes and test them, there is a `test pipeline`.
1. Make updates to the `data/test_list.csv` file to include the data you want, be sure to ingest data from a different place than the real data lives
2. Make your code changes
3. Run the `./test_pipeline.sh` script (with an optional `--delete` flag for deleting existing test data from previous runs)
    * This will save the data in `/local-scratch/projects/genotype-phenotype-map/test/data/`

## How it works

1. All the command line tools, python and R packages are installed in a docker image (see `docker` directory).
2. Each step in the pipeline is in the `pipeline_steps` directory
3. The pipeline is run using `Snakemake`.  You can see the `Snakefile` and `config.yaml` files to see how it runs.
4. The snakemake params are tuned to ieu-p1 specifically, namely that there are 256 CPUs, and maybe 200-300GB of RAM available (and not much else running on the box)

## Data

### Input

There are a few mandatory data directories
1. Precompiled ld matrices per ld region: `$DATA_DIR/ld_reference_panel_hg38/<ANCESTRY>/CHR/BP-RANGE.*`
3. Study data: can be anywhere on the box, specified in `data/study_list.csv`, however all study data is stored in `/local-scratch/data/`.  There does need to be a corresponding way to extract data from it though.
    * The folder and names of files matter.  The study id is constructed using the top level folder name and the name of the file 

### Output

Both data `DATA_DIR` and `RESULTS_DIR` directories will be populated

`DATA_DIR`
* `study`: each study will have `extracted_snps.tsv`, and a file in each of the 4 subdirectories related to the extraciton
  * `extracted`, `standardised`, `imputed`, and `finemapped`
* `ld_blocks`: each ld block will have data related to the extracted study regions inside the block, primarily for colocalisation
* `pipeline_metadata`:
  * `updated_ld_blocks_to_colocalise.tsv`: not all ld blocks need to be run, depending on the studies ingested.  This stores the list of ld regions that need to be run
  * `studies_to_process.tsv`: all studies that are currently being processed, if the pipeline is successful, these are appended to `results/studies_processed.tsv.gz`
  * `logs/snakemake.log`: output from the last time `./run_pipeline.sh` was run (or `./test_pipeline.sh`)

`RESULTS_DIR`
* `<timestamp> directory`: final results for the specific run of the pipeline. Latest timestamp is the latest result
* `latest/`: latest results for the pipeline, which include the following:
  * `studies_processed.tsv.gz`: latest iteration of studies that were processed
  * `traits_processed.tsv.gz`: latest iteration of traits that were processed
  * `associations.db`: associations database
  * `ld.db`: ld database
  * `gwas_upload.db`: gwas upload database
  * `studies.db`: studies database
