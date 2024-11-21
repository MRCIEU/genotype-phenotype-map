# Variant annotation

Variant annotation was performed for all variants used to generate the hg38 build LD reference panel using the Ensembl Variant Effect Predictor [VEP](https://www.ensembl.org/info/docs/tools/vep/index.html)

Variants were pulled from the UK Biobank imputed genotype TOPmed helper `.site.vcf.gz` files with filtering for inputation quality score and allele frequency as in `pull_hg38_sites.sh`. 

`update_varid.R` was used to reformat variant information to VEP ensembl input:\
`CHR POS_START POS_END REF/ALT STRAND NAME`\
All imputed variants are reported on the top strand (output file: `variant_info_hg38_forVEP.txt`). Variant name is in the format `chr_bp_A1_A2` with alleles ordered alphabetically.

`annotate_vep.sh` first extracts 9,681,625 variants used in the generation of the hg38 build LD reference panel (`variant_info_hg38_forVEP_runset.txt`) then calls `update_indelpositions.R` to update the chromosomal positions and ref/alt coding for short insertions and deletions as specified [here](https://www.ensembl.org/info/docs/tools/vep/vep_formats.html#default).

VEP was run with the parameters:

```
vep -i variant_info_hg38_forVEP_runset_posupdate.txt \
--cache \
--assembly GRCh38 \
--af_1kg \
--sift b \
--polyphen b \
--regulatory \
--show_ref_allele \
--symbol \
--protein \
--canonical \
--mane \
--biotype \
--check_ref \
--skipped_variants_file vep_skipped.txt \
--pick \
--no_stats \
--force_overwrite \
-v \
-o vep_variantannotations_hg38.txt
```

The `--pick` argument returns one consequence per variant as per the criteria described [here](https://www.ensembl.org/info/docs/tools/vep/script/vep_other.html#pick). Variants are preferentially annotated relative to the MANE select transcript, variant consequenced can be ranked according to [this table](https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html#consequences). The diagram below shows the location of each class of variant relative the transcript structure.

![](https://www.ensembl.org/info/genome/variation/prediction/consequences.svg)









