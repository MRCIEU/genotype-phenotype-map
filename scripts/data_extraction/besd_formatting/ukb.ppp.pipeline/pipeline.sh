#!/bin/bash

# Define env variables
export INPUT_DIR="/path/to/data/ukb-ppp/raw/ukb.ppp.pGWAS/european_discovery"
export PROCESSED_DIR="/path/to/pipeline_outdat/processed/ukb.european"
export SNPS_UPDATED_DIR="/path/to/pipeline_outdat/snps.updated/ukb.european"
export ESD_DIR="/path/to/pipeline_outdat/esd/ukb.european"
export FLIST="/path/to/pipeline_input/ukb.ppp.european.flist"
export BESD_OUT="/path/to/data/ukb-ppp/besd_formatted/european/ukb.ppp.european"

# Run scripts sequentially
./01.process.raw.sh

Rscript 02.update.snps.R

./03.concatenate.sh

./04.make.besd.sh
