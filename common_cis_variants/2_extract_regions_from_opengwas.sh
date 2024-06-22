#!/bin/bash
set -e
GWASES_TO_EXTRACT=$(cat $DATA_DIR/pipeline_metadata/gwases_to_extract.tsv)
LD_REGIONS=/home/common_cis_variants/data/ld_regions.tsv

while IFS= read -r GWAS_STUDY; do
  STUDY_DIR=$(echo $GWAS_STUDY | awk '{print $1}')
  EXTRACTION_DIR=$(echo $GWAS_STUDY | awk '{print $2}')
  ANCESTRY=$(echo $GWAS_STUDY | awk '{print $3}')
  SAMPLE_SIZE=$(echo $GWAS_STUDY | awk '{print $4}')
  P_VALUE=$(echo $GWAS_STUDY | awk '{print $5}')
  DATA_TYPE=$(echo $GWAS_STUDY | awk '{print $6}')

  mkdir -p $EXTRACTION_DIR/original $EXTRACTION_DIR/imputed $EXTRACTION_DIR/finemapped
  STUDY=$(basename $STUDY_DIR)
  echo $STUDY_DIR
  cd $STUDY_DIR

  EXTRACTED_SNPS=$EXTRACTION_DIR/extracted_snps.tsv
  echo -e "CHR\tBP\tLOG_P\tANCESTRY\tLD_REGION\tFILE\tCIS_TRANS" > $EXTRACTED_SNPS

  cat clump.txt
  ALL_CHR_POS=$(/home/bcftools/bcftools query -i 'ID=@clump.txt' --format "chr%CHROM %POS [%LP]\n" $STUDY.vcf.gz)
  echo $ALL_CHR_POS

  while IFS= read -r CHR_POS; do
    CHR=$(echo $CHR_POS | awk '{print substr($1,4)}')
    POS=$(echo $CHR_POS | awk '{print $2}')
    LOG_P=$(echo $CHR_POS | awk '{print $3}')
    if [[ -z $CHR ]] || [[ -z $POS ]] || [[ -z $LOG_P ]]; then
      continue
    fi

    BEGINNING_END=$(awk -v chr=$CHR -v bp=$POS -v an=$ANCESTRY \
      '{ if ($1 == chr && $2 < bp && bp < $3 && $4 == an) print $2 "-" $3 }' $LD_REGIONS
    )

    if [[ -z $BEGINNING_END ]]; then
      echo "Can't find LD region match for $ANCESTRY $CHR:$BP, skipping..."
      continue
    fi

    REGION="$CHR:$BEGINNING_END"
    echo "Region to extract: $REGION"
    EXTRACTED_FILE=$EXTRACTION_DIR/original/${ANCESTRY}_${CHR}_${POS}.tsv

    echo -e "RSID\tCHR\tBP\tEA\tOA\tEAF\tBETA\tSE\tLP" > $EXTRACTED_FILE
    /home/bcftools/bcftools query --regions $REGION --format "[%ID]\t[%CHROM]\t[%POS]\t[%REF]\t[%ALT]\t[%AF]\t[%ES]\t[%SE]\t[%LP]" $STUDY.vcf.gz >> $EXTRACTED_FILE

    SPECIFIC_LD_REGION="${ANCESTRY}/${CHR}_${BEGINNING_END//-/_}"
    LD_REGION_SNPLIST=$DATA_DIR/ld_block_matrices/${SPECIFIC_LD_REGION}.snplist
    echo -e "${CHR}\t${POS}\t${LOG_P}\t${ANCESTRY}\t${SPECIFIC_LD_REGION}\t${EXTRACTED_FILE}\tNA" >> $EXTRACTED_SNPS

    #TODO: maybe move this to the imputation or finemapping step
    #Rscript  /home/common_cis_variants/2_standardise_gwas_for_finemap.R --gwas_filename $EXTRACTED_FILE --snp_list $LD_REGION_SNPLIST
  done <<< $ALL_CHR_POS

  jq -n --arg sample_size $SAMPLE_SIZE --arg ancestry $ANCESTRY --arg p_val $P_VALUE --arg data_type $DATA_TYPE \
    '{"ancestry": $ancestry, "sample_size": $sample_size, "p_value_threshold":$p_val, "data_type":$data_type}' \
    > $EXTRACTION_DIR/extraction_metadata.json

  cd -
done <<< $GWASES_TO_EXTRACT
