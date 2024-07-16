# genotype-phenotype-map

Scoping document: https://uob.sharepoint.com/:w:/r/teams/grp-brs-tds/Shared%20Documents/The%20Human%20Genotype-Phenotype%20Map/Kickoff%20Document.docx?d=w3c2c1890b3bf4ec6a8e3883780b05c94&csf=1&web=1&e=KFn1PY


## How to use:

### 1. Clone the repo on to ieu-p1

`git clone git@github.com:MRCIEU/genotype-phenotype-map.git && cd genotype-phenotype-map`

### 2. Populate .env file
You will most likely want to use the default set in .env_example
`cp .env_example .env`

### 3. Add info in pipeline_steps/data/study_list.csv

This is the file that specifies which studies will be ingested by the pipeline, if you want to add a study or set of studies,
add a row to this csv indicating where they are on disk, and what script to ingest them with

### 4. ./run_pipeline.sh



## How to contribute 

There is not a fantastic way to have a test environment for the pipeline, but to test any addition to the pipline I would suggest:
1. Test the individual pipeline_step you are adding, using the data you want to pass in
2. There is a `./test_pipeline.sh`, which loads a tiny subset of different data, and writes to a different place.

