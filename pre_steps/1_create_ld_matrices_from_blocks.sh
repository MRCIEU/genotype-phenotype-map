#!/bin/bash

# location of plink format UKB European data
d="/local-scratch/projects/Lifecourse-GWAS/ukb/geno_input/"

# sample 50k randomly 
shuf -n 50000 $d/data.chr01.fam | cut -f 1,2 > ukb50k.ids


ldmat_region () {
    chr=$1
    chrl=$(printf "%02d" $chr)
    p1=$2
    p2=$3
    datadir=$4
    keeplist=$5
    plink=$6
    outdir=$7
    out="${outdir}/${chr}/${p1}-${p2}"
    f="$out.unphased.vcor1"
    if [ -e "$f" ]; then
        echo "done $out"
    else
        bfile="${datadir}/data.chr${chrl}"
        tfile=$(mktemp)
        echo "$chr $p1 $p2 $chr_$p2_$p2" > $tfile
        mkdir -p ${outdir}/${chr}
        $plink --bfile $bfile --chr $chr --extract range $tfile --r-unphased square --out $out --keep-allele-order
        $plink --bfile $bfile --chr $chr --extract range $tfile  --out $out  --keep-allele-order --freq
        $plink --bfile $bfile --chr $chr --extract range $tfile  --out $out  --keep-allele-order --make-just-bim
    fi
}

ldmat_region 1 1892607 3582736 $d ukb50k.ids /local-scratch/projects/genotype-phenotype-map/bin/plink2 temp
ldmat_region 22 44995308 46470495 $d ukb50k.ids /local-scratch/projects/genotype-phenotype-map/bin/plink2 temp

{
    while IFS=$'\t' read -r CHR START_BP END_BP; do
        ldmat_region $CHR $START_BP $END_BP $d ukb50k.ids /local-scratch/projects/genotype-phenotype-map/bin/plink2 temp
    done
} < eur_ldregions.txt

# parallel -j 60 -a eur_ldregions.txt ./ldmat.sh {} $d ukb50k.ids /local-scratch/projects/genotype-phenotype-map/bin/plink2 temp

ls -1 temp/*/*vcor1 | wc -l
wc -l eur_ldregions.txt
