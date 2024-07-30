#!/bin/bash

set -e

DCODE_DATA=${RAWDATA_DIR}/ukb-seq/downloads/halldorexwas/decode_data
LD_REGIONS=${home}/seqdata_steps/pre_steps/data/pyrho_EUR_LD_blocks.bed

# Lower MAF threshold (for 10 expected observations of minor allele in 431,805 white British)
MIN_MAF=0.000012

CHR="Chrom"
POS="Pos"
P="pval"
MAF="effectAlleleFreq"
INFO="info"

OUT_DIR=${RAWDATA_DIR}/ukb-seq/downloads/halldorexwas/decode_data_filtered
mkdir -p ${OUT_DIR}

for FILE in $(ls ${DCODE_DATA}); do

    start_time=$(date +%s)

    BASENAME=$(basename ${FILE} .txt.gz)
    OUTNAME=${BASENAME}_filtered.txt

    echo "Processing: ${FILE}"

    if [[ ! -e "${OUT_DIR}/${OUTNAME}.gz" ]] ; then

        # Indicies of required columns (annoyingly the order is inconsistent)
        read -r header < ${DECODE_DATA}/${FILE}
        chr_index=$(echo $header | awk -v col="$CHR" '{for (i=1; i<=NF; i++) if ($i == col) print i}')
        pos_index=$(echo $header | awk -v col="$POS" '{for (i=1; i<=NF; i++) if ($i == col) print i}')
        p_index=$(echo $header | awk -v col="$P" '{for (i=1; i<=NF; i++) if ($i == col) print i}')
        maf_index=$(echo $header | awk -v col="$MAF" '{for (i=1; i<=NF; i++) if ($i == col) print i}')
        info_index=$(echo $header | awk -v col="$INFO" '{for (i=1; i<=NF; i++) if ($i == col) print i}')

        if [[ -z $chr_index ]] || [[ -z $pos_index ]] || [[ -z $p_index ]] || [[ -z $maf_index ]] || [[ -z $info_index ]] ; then
            echo "One or more columns not found in the header"
            exit 1
        fi
            
        # Filter 1 : MAF<=0.01, MAF>=0.000012, INFO>=0.5, P<=0.00005
        # Extract chr and pos of rare variants passing filter

        awk -v min_maf=${MIN_MAF} -v chr=${chr_index} -v pos=${pos_index} -v p=${p_index} -v maf=${maf_index} -v info=${info_index} \
        '$p <= 0.00005 && $info >= 0.5 && (($maf <= 0.01 && $maf >= min_maf) || (1-$maf <= 0.01 && 1-$maf >= min_maf)) \
        {print $chr, $pos}' ${DCODE_DATA}/${FILE} > ${OUT_DIR}/top_tmp.txt

        # Keep unique (some variant position duplication in file)   
        sort -n -k1,1 -k2,2 ${OUT_DIR}/top_tmp.txt | uniq > ${OUT_DIR}/top_tmp_uniq.txt

        # Loop through variants passing filter 1, extract 1MB surrounding region and apply filter 2
        # Filter 2 : MAF<=0.01, MAF>=0.000012, INFO>=0.5, P<=0.01

        while IFS=' ' read -r chr bp; do

            awk -v CHR=${chr} -v BP=${bp} \
            -v min_maf=${MIN_MAF} -v chr=${chr_index} -v pos=${pos_index} -v p=${p_index} -v maf=${maf_index} -v info=${info_index} \
            'BEGIN {min = BP-500000; max = BP+500000} \
            $chr == CHR && $pos >= min && $pos <= max && \
            $p <= 0.01 && $info >= 0.5 && (($maf <= 0.01 && $maf >= min_maf) || (1-$maf <= 0.01 && 1-$maf >= min_maf))' \
            ${DCODE_DATA}/${FILE} >> ${OUT_DIR}/window_tmp.txt

        done < ${OUT_DIR}/top_tmp_uniq.txt

        # Keep unique (duplicates from overlapping windows) and write to output
        sort -n -k1,1 -k2,2 window_tmp.txt | uniq > ${OUT_DIR}/${OUTNAME}
        gzip ${OUT_DIR}/${OUTNAME}

        rm ${OUT_DIR}/top_tmp.txt  ${OUT_DIR}/top_tmp_uniq.txt ${OUT_DIR}/window_tmp.txt

    else
        
        echo "File ${OUT_DIR}/${OUTNAME}.gz already exists"

    fi

    end_time=$(date +%s)
    echo "Duration: ($start_time-$end_time)/3600 hour(s)"

done