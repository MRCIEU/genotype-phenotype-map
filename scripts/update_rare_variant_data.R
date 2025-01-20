backman_files <- Sys.glob('/local-scratch/data/ukb-seq/downloads/backmanexwas/*.tsv.gz')
backman_columns <- 
Name    chromosome      base_pair_location      other_allele    effect_allele   Trait   Cohort  Model   odds_ratio      ci_lower        ci_upper        p_value effect_allele_frequency standard_error
lapply(backman_files, function(file) {
  rv <- vroom::vroom(file, show_column_types=F) |>
    dplyr::rename(CHR=chromosome, BP=base_pair_location, P=p_value, EA=effect_allele, OA=other_allele, trait=Trait, OR=odds_ratio, CI_LOWER=ci_lower, CI_UPPER=ci_upper, EAF=effect_allele_frequency, SE=standard_error)
  vroom::vroom_write(rv, file)
}

genebass_files <- Sys.glob('/local-scratch/data/ukb-seq/downloads/genebass/rare/*.tsv')
locus	alleles	markerID	gene	annotation	call_stats	n_cases	n_controls	heritability	inv_normalized	trait_type	phenocode	pheno_sex	coding	modifier	n_cases_defined	n_cases_both_sexes	n_cases_females	n_cases_males	description	coding_description	AC	AF	BETA	SE	AF.Cases	AF.Controls	Pvalue
lapply(genebass_files, function(file) {
  rv <- vroom::vroom(file, show_column_types=F) |>
    tidyr::separate(col=markerID, into=c("CHR", "BP", "EA", "OA"), sep='chr:_/') |>
    dplyr::rename(P=Pvalue, EAF=AF)
  file <- sub('.tsv$', '.tsv.gz', file)
  vroom::vroom_write(rv, file)
})

az_dir <- Sys.glob('/local-scratch/data/ukb-seq/downloads/azexwas/*.csv.gz')
Variant,Variant type,Phenotype,Category,Model,Consequence type,Gene,Transcript,cDNA change,Amino acid change,
Exon rank,No. samples,No. AA genotypes,No. AB genotypes,No. BB genotypes,AAF,% AB or BB genotypes,% BB genotypes,
p-value,Effect size,Effect size standard error,Effect size LCI,Effect size UCI,Median value (cases),
Median value (controls),MAF
lapply(genebass_files, function(file) {
  rv <- vroom::vroom(file, show_column_types=F) |>
    tidyr::separate(col=Variant, into=c("CHR", "BP", "EA", "OA"), sep='-') |>
    dplyr::rename(P='p-value', EAF=MAF, BETA='Effect size', SE='Effect size standard error')
  file <- sub('csv.gz', 'tsv.gz', file)
  vroom::vroom_write(rv, file)
}

