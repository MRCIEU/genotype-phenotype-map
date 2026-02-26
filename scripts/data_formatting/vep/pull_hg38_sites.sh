#!/bin/bash

source ~/config.sh

#dx upload update_varid.R --destination="${project}:/Scripts/"

# Extract sites list from UKB TopMed imputation files with filters
imp_dir="/mnt/project/Bulk/Imputation/Imputation\ from\ genotype\ \(TOPmed\)/helper_files"

cmd="

cp ${imp_dir}/*.sites.vcf.gz .

for chr in {1..22}; do

# Get the list of variants that pass the threshold
bcftools filter -e 'INFO/R2<0.8' ukb21007_c\${chr}_b0_v1.sites.vcf.gz | \
bcftools filter -e 'INFO/AF<0.001' | \
bcftools filter -e 'INFO/AF>0.999' | \
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%AF\t%R2\n' > variant_info_chr\${chr}.txt
done

cat variant_info_chr*.txt > variant_info_hg38.txt

# Format for VEP input
awk '{print substr(\$1,4), \$2, \$2, \$4\"/\"\$5, \"+\"}' variant_info_hg38.txt > variant_info_hg38_forVEP.txt

# Add alphabetised identifiers
Rscript update_varid.R variant_info_hg38_forVEP.txt

rm variant_info_chr*.txt
rm *.sites.vcf.gz
"

dx run swiss-army-knife \
    -icmd="${cmd}" \
	-iin="${project}:/Scripts/update_varid.R" \
    --tag="hg38 sites" \
	--destination="${project}:/Data/UKBGenotypes/Imputed" \
    --brief \
    --yes \
    --instance-type="mem2_ssd1_v2_x8"
