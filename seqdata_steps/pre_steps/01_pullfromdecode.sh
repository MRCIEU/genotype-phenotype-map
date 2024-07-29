#!/bin/bash

set -e

if [[$# -ne 1]]; then
    echo "Provide file listing URLs for download"
    exit 0
fi

LIST=$1
# This path is wrong
LD_REGIONS=/home/pipeline/pipeline_steps/data/ld_regions.tsv
ANCESTRY="EUR"

cd ${RAWDATA_DIR}/ukb-seq/downloads/halldorexwas
mkdir -p decode_data_filtered

for LINK in $(cat ${LIST}); do
        FILENAME=$(echo ${LINK} | cut -d"-" -f3)
        BASENAME=$(basename ${FILENAME} .txt.gz)

        echo "Downloading: ${FILENAME}"
        wget "${LINK}" | gunzip > WGS_tmp.txt

        # Lower MAF threshold (for 10 expected observations of minor allele in 431,805 white British)
        MIN_MAF=0.000012

        # Extract chr and pos of rare variants passing filter (MAF<=0.01, MAF>=0.000012, INFO>=0.5, P<=0.01)
        awk -v min_maf=${MIN_MAF} \
        '$9<=0.01 && $12>=0.5 && (($7<=0.01 && $7>=min_maf) || (1-$7<=0.01 && 1-$7>=min_maf)) \
        {print $1, $2}' WGS_tmp.txt > top_tmp.txt

        sed -i '' 's/chr//g' top_tmp.txt

        # Identify regions containing at least one variant passing filter
        while IFS=' ' read -r CHR POS; do
            if [[ -z $CHR]] || [[ -z $POS]]; then
                continue
            fi

        REGION=$(awk -v chr=$CHR -v bp=$POS -v an=$ANCESTRY \
            '{ if ($1 == chr && $2 < bp && bp < $3 && $4 == an) print $1, $2, $3}' $LD_REGIONS)

        echo ${REGION} >> regions_tmp.txt
        done < top_tmp.txt

        ## CONTINUE






        





        # Filter for single row in given region (against reference set of regions)
        # De-duplicate list

        # Loop through variants and extract variants in 1MB surrounding region

        while IFS=' ' read -r chr bp; do
                awk -v CHR=${chr} -v BP=${bp} '
                BEGIN {min=BP-500000; max = BP+500000}
                $1==CHR && $2>=min && $2<=max' 
                WGS_tmp.txt >> window_tmp.txt
        done < top_tmp.txt

        # Excude duplicate lines (probs a much smarter way to do this)
        OUTFILE="${basename}_p1e5windows.txt"
        
        $(sort window_tmp.txt | uniq) > ${OUTFILE}

        gzip ${OUTFILE}

        # Remove temporary files
        rm WGS_tmp.txt top_tmp.txt
done