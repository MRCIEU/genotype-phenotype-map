#!/bin/bash
set -e
MEGABASE=1000000 # +/- 1 Mega base
GWASES_TO_EXTRACT=$(cat $DATA_DIR/pipeline_metadata/gwases_to_extract.tsv)

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
  echo -e "CHR\tBP\tLOG_P\tANCESTRY" > $EXTRACTED_SNPS

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
    BEGINNING=$(($POS-$MEGABASE))

    if [[ $BEGINNING -lt 0 ]]; then
      BEGINNING=0
    fi
    END=$(($POS+$MEGABASE))

    REGION="$CHR:$BEGINNING-$END"
    echo "Region to extract: $REGION"
    EXTRACTED_FILE=$EXTRACTION_DIR/${ANCESTRY}_${CHR}_${POS}.z

    echo -e "rsid chromosome position allele1 allele2 maf beta se lp" > $EXTRACTED_FILE
    /home/bcftools/bcftools query --regions $REGION --format "[%ID] [%CHROM] [%POS] [%REF] [%ALT] [%AF] [%ES] [%SE] [%LP]" $STUDY.vcf.gz >> $EXTRACTED_FILE
    echo -e "${CHR}\t${POS}\t${LOG_P}\t${ANCESTRY}" >> $EXTRACTED_SNPS
  done <<< $ALL_CHR_POS

  jq -n --arg sample_size $SAMPLE_SIZE --arg ancestry $ANCESTRY --arg p_val $P_VALUE \
    '{"ancestry": $ancestry, "sample_size": $sample_size, "p_value_threshold":$p_val}' > $EXTRACTION_DIR/extraction_metadata.json

  echo "Step 2: Finemap extracted regions"

  FINEMAP_DIR=$EXTRACTION_DIR/finemap
  PLINK_DIR=$EXTRACTION_DIR/plink
  IN_FILE=$FINEMAP_DIR/finemap_infile.txt
  mkdir -p $FINEMAP_DIR
  mkdir -p $PLINK_DIR
  echo "z;ld;snp;config;cred;n_samples" > $IN_FILE

  cd $EXTRACTION_DIR
  for FILE in $(ls | grep $ANCESTRY); do
    echo $FILE
    FILE_BASENAME=$(basename $FILE .z)
    INFO=(${FILE_BASENAME//_/ })
    FINEMAP_OUTPUT=$FINEMAP_DIR/$FILE_BASENAME
    PLINK_OUTPUT=$PLINK_DIR/$FILE_BASENAME

    SNP_LIST=$(grep -v rsid $FILE | awk '{print $1}')
    SNP_LIST_FILE=${FINEMAP_OUTPUT}.snplist
    echo -e $SNP_LIST > $SNP_LIST_FILE

    plink1.9 --bfile $THOUSAND_GENOMES/$ANCESTRY --chr ${INFO[1]} --extract $SNP_LIST_FILE --r square spaces --out $PLINK_OUTPUT --write-snplist --keep-allele-order
    Rscript  /home/common_cis_variants/2_standardise_gwas_for_finemap.R --gwas_filename $FILE --snp_list $PLINK_OUTPUT.snplist

    echo "$EXTRACTION_DIR/$FILE;$PLINK_OUTPUT.ld;$FINEMAP_OUTPUT.snp;$FINEMAP_OUTPUT.config;$FINEMAP_OUTPUT.cred;$SAMPLE_SIZE" >> $IN_FILE
  done

  #TODO: looking at other finemapping options...
  #finemap --sss --in-files $IN_FILE

  cd -
done <<< $GWASES_TO_EXTRACT
