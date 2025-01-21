backman_files <- Sys.glob('/local-scratch/data/ukb-seq/downloads/backmanexwas/*.tsv.gz')
backman_dir <- '/local-scratch/data/hg38/backman'
# backman_columns <- 
# Name    chromosome      base_pair_location      other_allele    effect_allele   Trait   Cohort  Model   odds_ratio      ci_lower        ci_upper        p_value effect_allele_frequency standard_error
lapply(backman_files, function(file) {
  new_file_name <- paste0(backman_dir, tools::file_path_sans_ext(file), 'tsv.gz')
  rv <- vroom::vroom(file)
  if ("beta" %in% colnames(rv)) {
    rv <- dplyr::rename(rv, CHR=chromosome, BP=base_pair_location, P=p_value, EA=effect_allele, OA=other_allele, trait=Trait, BETA=beta, EAF=effect_allele_frequency, SE=standard_error)
  } else if ("odds_ratio" %in% colnames(rv)) {
    rv <- dplyr::rename(rv, CHR=chromosome, BP=base_pair_location, P=p_value, EA=effect_allele, OA=other_allele, trait=Trait, OR=odds_ratio, CI_LOWER=ci_lower, CI_UPPER=ci_upper, EAF=effect_allele_frequency)
  }
  vroom::vroom_write(rv, new_file_name)
}

genebass_files <- Sys.glob('/Users/wt23152/Documents/Projects/scratch/021/data/study/rare/genebass*.tsv')
genebass_dir <- '/Users/wt23152/Documents/Projects/scratch/021/data/study/rare/new/'
# genebass_dir <- '/local-scratch/data/hg38/genebass_filtered'
lapply(genebass_files, function(file) {
  new_file_name <- paste0(genebass_dir, basename(tools::file_path_sans_ext(file)), 'tsv.gz')
  print(new_file_name)
  rv <- vroom::vroom(file) |>
    tidyr::separate(col=markerID, into=c("CHR", "BP", "OA", "EA"), sep='[:_/]') |>
    dplyr::mutate(CHR=sub('chr', '', CHR)) |>
    dplyr::rename(P=Pvalue, EAF=AF)
  file <- sub('.tsv$', '.tsv.gz', file)
  vroom::vroom_write(rv, new_file_name)
})

#if you have binary they give you 
# if they give you odds ratio column present use control AAF
# if there is'nt odds ratio column, use AAF?
az_files <- Sys.glob('/local-scratch/data/ukb-seq/downloads/azexwas/*.csv.gz')
Variant,Variant type,Phenotype,Category,Model,Consequence type,Gene,Transcript,cDNA change,Amino acid change,
Exon rank,No. samples,No. AA genotypes,No. AB genotypes,No. BB genotypes,AAF,% AB or BB genotypes,% BB genotypes,
p-value,Effect size,Effect size standard error,Effect size LCI,Effect size UCI,Median value (cases),
Median value (controls),MAF

az_dir <- '/local-scratch/data/hg38/azexwas'
lapply(az_files, function(file) {
  new_file_name <- paste0(az_dir, tools::file_path_sans_ext(file), 'tsv.gz')
  rv <- vroom::vroom(file)
  if ("beta" %in% colnames(rv)) {
    rv <- tidyr::separate(rv, col=Variant, into=c("CHR", "BP", "OA", "EA"), sep='-') |>
      dplyr::rename(P='p-value', EAF=AAF, BETA='Effect size', SE='Effect size standard error')
  } else if () {
    rv <- tidyr::separate(rv, col=Variant, into=c("CHR", "BP", "OA", "EA"), sep='-') |>
      dplyr::rename(P='p-value', EAF=MAF, BETA='Effect size', SE='Effect size standard error')
  }
  vroom::vroom_write(rv, new_file_name)
}
