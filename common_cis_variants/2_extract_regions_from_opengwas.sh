#!/bin/bash
set -e
MEGABASE=1000000 # +/- 1 Mega base #TODO not using this anymore, just grabbing whole ld region
GWASES_TO_EXTRACT=$(cat $DATA_DIR/pipeline_metadata/gwases_to_extract.tsv)
LD_REGIONS=data/ld_regions.tsv

while IFS= read -r GWAS_STUDY; do
  STUDY_DIR=$(echo $GWAS_STUDY | awk '{print $1}')
  EXTRACTION_DIR=$(echo $GWAS_STUDY | awk '{print $2}')
  ANCESTRY=$(echo $GWAS_STUDY | awk '{print $3}')
  SAMPLE_SIZE=$(echo $GWAS_STUDY | awk '{print $4}')
  P_VALUE=$(echo $GWAS_STUDY | awk '{print $5}')

  STUDY=$(basename $STUDY_DIR)
  echo $STUDY_DIR
  cd $STUDY_DIR

  mkdir -p $EXTRACTION_DIR
  EXTRACTED_SNPS=$EXTRACTION_DIR/extracted_snps.tsv
  echo -e "CHR\tBP\tLOG_P\tANCESTRY\tLD_REGION\tFILE" > $EXTRACTED_SNPS

  echo "Step 1: Extract regions using bcftools"

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

    REGION="$CHR:$BEGINNING_END"
    echo "Region to extract: $REGION"
    EXTRACTED_FILE=$EXTRACTION_DIR/${ANCESTRY}_${CHR}_${POS}.z

    echo -e "rsid chromosome position allele1 allele2 maf beta se lp" > $EXTRACTED_FILE
    /home/bcftools/bcftools query --regions $REGION --format "[%ID] [%CHROM] [%POS] [%REF] [%ALT] [%AF] [%ES] [%SE] [%LP]" $STUDY.vcf.gz >> $EXTRACTED_FILE

    SPECIFIC_LD_REGION="${ANCESTRY}/${CHR}_${BEGINNING_END//-/_}}"
    LD_REGION_SNPLIST=$DATA_DIR/ld_block_matrices/${SPECIFIC_LD_REGION}.snplist
    echo -e "${CHR}\t${POS}\t${LOG_P}\t${ANCESTRY}\t${SPECIFIC_LD_REGION}\t${EXTRACTED_FILE}" >> $EXTRACTED_SNPS

    #TODO: maybe move this to the imputation or finemapping step
    Rscript  /home/common_cis_variants/2_standardise_gwas_for_finemap.R --gwas_filename $EXTRACTED_FILE --snp_list $LD_REGION_SNPLIST
  done <<< $ALL_CHR_POS

  jq -n --arg sample_size $SAMPLE_SIZE --arg ancestry $ANCESTRY --arg p_val $P_VALUE \
    '{"ancestry": $ancestry, "sample_size": $sample_size, "p_value_threshold":$p_val}' > $EXTRACTION_DIR/extraction_metadata.json

  cd -
done <<< $GWASES_TO_EXTRACT
