# https://www.cog-genomics.org/plink/1.9/input#make_set
#extract snp list of range, then run the --r command

plink1.9 --bfile $THOUSAND_GENOMES/$ANCESTRY --chr ${INFO[1]} --extract $SNP_LIST_FILE --r square spaces --out $PLINK_OUTPUT --write-snplist --keep-allele-order
