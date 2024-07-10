#!/bin/bash
set -e

if [[ $# -ne 8 ]] ; then
  echo "Incorrect number of arguments: need 8"
  exit 0
fi

STUDY_NAME=$1
ANCESTRY=$2
SAMPLE_SIZE=$3
DATA_TYPE=$4
ORIG_STUDY_DIR=$5
EXTRACTED_STUDY_DIR=$6
P_VALUE=$7
GENE=$8

LD_REGIONS=/home/common_cis_variants/data/ld_regions.tsv
STUDY=$(basename $ORIG_STUDY_DIR)
cd $ORIG_STUDY_DIR

mkdir -p $EXTRACTED_STUDY_DIR/original $EXTRACTED_STUDY_DIR/imputed $EXTRACTED_STUDY_DIR/finemapped

EXTRACTED_SNPS=$EXTRACTED_STUDY_DIR/extracted_snps.tsv
echo -e "chr\tbp\tlog_p\tancestry\tld_region\tfile\tcis_trans" > $EXTRACTED_SNPS

echo "RSIDs to extract for ${STUDY_NAME}: $(wc -l < clump.txt)"
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
  EXTRACTED_FILE=$EXTRACTED_STUDY_DIR/original/${ANCESTRY}_${CHR}_${POS}.tsv

  echo -e "RSID\tCHR\tBP\tEA\tOA\tEAF\tBETA\tSE\tLP" > $EXTRACTED_FILE
  /home/bcftools/bcftools query --regions $REGION --format "[%ID]\t[%CHROM]\t[%POS]\t[%REF]\t[%ALT]\t[%AF]\t[%ES]\t[%SE]\t[%LP]" $STUDY.vcf.gz >> $EXTRACTED_FILE

  gzip $EXTRACTED_FILE
  EXTRACTED_FILE="$EXTRACTED_FILE.gz"

  SPECIFIC_LD_REGION="${ANCESTRY}/${CHR}_${BEGINNING_END//-/_}"

  #TODO: delete later, once cis_trans stuff is in
  CIS_TRANS=NA
  if [[ $STUDY =~ 'eqtl'  ]]; then
    CIS_TRANS=cis
  fi
  echo -e "${CHR}\t${POS}\t${LOG_P}\t${ANCESTRY}\t${SPECIFIC_LD_REGION}\t${EXTRACTED_FILE}\t${CIS_TRANS}" >> $EXTRACTED_SNPS
done <<< $ALL_CHR_POS

