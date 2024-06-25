#!/bin/bash
set -e
if [[ $# -ne 6 ]] ; then
  echo "Incorrect number of arguments: need 6"
  exit 0
fi

ANCESTRY=$1
SAMPLE_SIZE=$2
DATA_TYPE=$3
ORIG_STUDY_DIR=$4
EXTRACTED_STUDY_DIR=$5
P_VALUE=$6

LD_REGIONS=/home/common_cis_variants/data/ld_regions.tsv
STUDY=$(basename $ORIG_STUDY_DIR)
echo $ORIG_STUDY_DIR
cd $ORIG_STUDY_DIR

mkdir -p $EXTRACTED_STUDY_DIR/original $EXTRACTED_STUDY_DIR/imputed $EXTRACTED_STUDY_DIR/finemapped

EXTRACTED_SNPS=$EXTRACTED_STUDY_DIR/extracted_snps.tsv
echo -e "CHR\tBP\tLOG_P\tANCESTRY\tLD_REGION\tFILE\tCIS_TRANS" > $EXTRACTED_SNPS

echo "RSIDs to extract:"
cat clump.txt
ALL_CHR_POS=$(/home/bcftools/bcftools query -i 'ID=@clump.txt' --format "chr%CHROM %POS [%LP]\n" $STUDY.vcf.gz)

while IFS=' ' read -r CHR POS LOG_P; do
  CHR=$(echo $CHR| awk '{print substr($1,4)}')
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
  EXTRACTED_FILE=$EXTRACTED_STUDY_DIR/original/${ANCESTRY}_${CHR}_${POS}.tsv

  echo -e "RSID\tCHR\TBP\tEA\tOA\tEAF\tBETA\tSE\tLP" > $EXTRACTED_FILE
  /home/bcftools/bcftools query --regions $REGION --format "[%ID]\t[%CHROM]\t[%POS]\t[%REF]\t[%ALT]\t[%AF]\t[%ES]\t[%SE]\t[%LP]" $STUDY.vcf.gz >> $EXTRACTED_FILE

  SPECIFIC_LD_REGION="${ANCESTRY}/${CHR}_${BEGINNING_END//-/_}}"
  LD_REGION_SNPLIST=$DATA_DIR/ld_block_matrices/${SPECIFIC_LD_REGION}.snplist
  echo -e "${CHR}\t${POS}\t${LOG_P}\t${ANCESTRY}\t${SPECIFIC_LD_REGION}\t${EXTRACTED_FILE}\tNA" >> $EXTRACTED_SNPS
done <<< $ALL_CHR_POS

jq -n --arg sample_size $SAMPLE_SIZE --arg ancestry $ANCESTRY --arg p_val $P_VALUE --arg data_type $DATA_TYPE \
  '{"ancestry": $ancestry, "sample_size": $sample_size, "p_value_threshold":$p_val, "data_type":$data_type}' \
  > $EXTRACTED_STUDY_DIR/extraction_metadata.json
