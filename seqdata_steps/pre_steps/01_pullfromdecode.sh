#!/bin/bash

set -e

# List of file names
LIST=$1

RAWDATA_DIR=/local-scratch/data/
DCODE_DATA=${RAWDATA_DIR}/ukb-seq/downloads/halldorexwas/decode_data

# Lower MAF threshold (for 10 expected observations of minor allele in 431,805 white British)
MIN_MAF=0.000012

CHR="Chrom"
POS="Pos"
P="pval"
MAF="effectAlleleFreq"
INFO="info"

OUT_DIR=${RAWDATA_DIR}/ukb-seq/downloads/halldorexwas/decode_data_filtered
mkdir -p ${OUT_DIR}

for FILE in $(cat ${LIST}); do

    start_time=$(date +%s)

    BASENAME=$(basename ${FILE} .txt.gz)
    OUTNAME=${BASENAME}_filtered.txt

    #echo "Checking integrity of ${FILE}"
    #CALC_CHEKSUM=$(md5sum ${DCODE_DATA}/${FILE})
    #CHECKSUM=$(cat ${DCODE_DATA}/${BASENAME}.md5sum)
    #[[ ${CALC_CHECKSUM} == ${CHECKSUM} ]] || exit 1

    echo "Unzipping ${FILE}..."
    
    gunzip ${DCODE_DATA}/${FILE}
    INFILE=${BASENAME}.txt

    echo "Processing file: ${INFILE}"

    if [[ ! -e "${OUT_DIR}/${OUTNAME}.gz" ]] ; then

	# Temporary directory for processing
    	TMP_DIR=${OUT_DIR}/tmp_${BASENAME}
    	mkdir -p ${TMP_DIR}

        # Indicies of required columns (annoyingly the order is inconsistent)
        read -r header < ${DCODE_DATA}/${INFILE}
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
        echo "Extracting chr and pos of rare variants passing filter 1"

        awk -v min_maf=${MIN_MAF} -v chr=${chr_index} -v pos=${pos_index} -v p=${p_index} -v maf=${maf_index} -v info=${info_index} \
        '$p <= 0.00005 && $info >= 0.5 && (($maf <= 0.01 && $maf >= min_maf) || (1-$maf <= 0.01 && 1-$maf >= min_maf)) \
        {print $chr, $pos}' ${DCODE_DATA}/${INFILE} > ${TMP_DIR}/top_tmp.txt

        # Keep unique (some variant position duplication in file)
        sort -V -k1,1 -k2,2 ${TMP_DIR}/top_tmp.txt | uniq > ${TMP_DIR}/top_tmp_uniq.txt

	echo "No. filter 1 hits: $(wc -l < ${TMP_DIR}/top_tmp_uniq.txt)"
	
	# Grep surrounding 200,000 lines each side of filter 1 variants and split files
	echo "Splitting windows surrounding filter 1 hits"
	
	i=1
	while IFS= read -r line; do
		keepchr=$(echo ${line} | awk '{print $1}')
		pasteline=$(echo ${line} | awk '{print $1":"$2":"}')

		rg ${pasteline} ${DCODE_DATA}/${INFILE} -C 200000 | awk -v chr=${chr_index} -v keep=${keepchr} \
		'$chr == keep' > "${TMP_DIR}/searchwindow_${i}_tmp.txt"
		
		i=$((i+1))
	done < ${TMP_DIR}/top_tmp_uniq.txt

        # Loop through variants passing filter 1, extract 1MB surrounding region and apply filter 2
	# Run in parallel over split files
        # Filter 2 : MAF<=0.01, MAF>=0.000012, INFO>=0.5, P<=0.01
	
	echo "Extracting filter 2 variants"

	i=$(wc -l < ${TMP_DIR}/top_tmp_uniq.txt)
	
	for ((j=1; j<=i; j++)); do
		
		# Position of the original filtered variant
		refpos=$(awk -v pos=${pos_index} -v J=${j} 'NR==J {print $pos}' ${TMP_DIR}/top_tmp_uniq.txt)

            	awk -v min_maf=${MIN_MAF} -v ref=${refpos} -v pos=${pos_index} \
		-v p=${p_index} -v maf=${maf_index} -v info=${info_index} \
            	'BEGIN {min = ref-500000; max = ref+500000} \
            	$pos >= min && $pos <= max && \
            	$p <= 0.01 && $info >= 0.5 && (($maf <= 0.01 && $maf >= min_maf) || (1-$maf <= 0.01 && 1-$maf >= min_maf))' \
            	${TMP_DIR}/searchwindow_${j}_tmp.txt > ${TMP_DIR}/serchwindow_${j}_filtered_tmp.txt
        done

	# Join split filtered files
	
	echo "Joining split filtered files"
	
	$(echo ${header}; cat ${TMP_DIR}/serchwindow_*_filtered_tmp.txt) > ${TMP_DIR}/finalvars_tmp.txt

        # Keep unique (duplicates from overlapping windows) add header and write to output

	echo "Writing output to: ${OUT_DIR}/${OUTNAME}"

	sort -V -k1,1 -k2,2 ${TMP_DIR}/finalvars_tmp.txt | uniq > ${OUT_DIR}/${OUTNAME}
        
	echo "Final variant count: $(wc -l < ${OUT_DIR}/${OUTNAME})"

	gzip ${OUT_DIR}/${OUTNAME}

	# Re-zip original file
	gzip ${DCODE_DATA}/${INFILE}

        #rm -r ${TMP_DIR}

    else
        
        echo "File ${OUT_DIR}/${OUTNAME}.gz already exists"

    fi

    end_time=$(date +%s)
    duration=$((end_time - start_time))
    echo "Duration: $((duration / 3600)) hours(s)"

done
