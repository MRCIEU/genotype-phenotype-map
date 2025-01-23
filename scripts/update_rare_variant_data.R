# backman_files <- Sys.glob('/local-scratch/data/ukb-seq/downloads/backmanexwas/*.tsv.gz')
# backman_dir <- '/local-scratch/data/hg38/backman'
# file.copy('/local-scratch/data/ukb-seq/downloads/backmanexwas/metadata.tsv', backman_dir)
backman_files <- Sys.glob('/Users/wt23152/Documents/Projects/scratch/021/data/study/rare/GCST*')
backman_dir <- '/Users/wt23152/Documents/Projects/scratch/021/data/study/rare/new/'
file.copy('/Users/wt23152/Documents/Projects/scratch/021/data/study/rare/backman_metadata.tsv', backman_dir)
lapply(backman_files, function(file) {
  new_basename <- basename(sub("\\..*$", "", file))
  new_basename <- sub('_buildGRCh38', '', new_basename)
  new_file_name <- paste0(backman_dir, new_basename, '.tsv.gz')
  rv <- vroom::vroom(file)
  if ("beta" %in% colnames(rv)) {
    rv <- dplyr::rename(rv, CHR=chromosome, BP=base_pair_location, P=p_value, EA=effect_allele, OA=other_allele, trait=Trait, BETA=beta, EAF=effect_allele_frequency, SE=standard_error)
  } else if ("odds_ratio" %in% colnames(rv)) {
    rv <- dplyr::rename(rv, CHR=chromosome, BP=base_pair_location, P=p_value, EA=effect_allele, OA=other_allele, trait=Trait, OR=odds_ratio, CI_LOWER=ci_lower, CI_UPPER=ci_upper, EAF=effect_allele_frequency)
  }
  vroom::vroom_write(rv, new_file_name)
})

# genebass_files <- Sys.glob('/local-scratch/data/ukb-seq/downloads/genebass/*.tsv.gz')
# genebass_dir <- '/local-scratch/data/hg38/genebass_filtered'
# file.copy('/local-scratch/data/ukb-seq/downloads/genebass/metadata.tsv', genebass_dir)
genebass_files <- Sys.glob('/Users/wt23152/Documents/Projects/scratch/021/data/study/rare/genebass*.tsv')
genebass_dir <- '/Users/wt23152/Documents/Projects/scratch/021/data/study/rare/new/'
file.copy('/Users/wt23152/Documents/Projects/scratch/021/data/study/rare/backman_metadata.tsv', backman_dir)
lapply(genebass_files, function(file) {
  new_basename <- basename(sub("\\..*$", "", file))
  new_basename <- sub('genebass_ukbwes_', '', new_basename)
  new_basename <- sub('_rare_filtered', '', new_basename)
  new_file_name <- paste0(genebass_dir, new_basename, '.tsv.gz')

  print(new_file_name)
  rv <- vroom::vroom(file) |>
    tidyr::separate(col=markerID, into=c("CHR", "BP", "OA", "EA"), sep='[:_/]', remove = F) |>
    dplyr::mutate(CHR=sub('chr', '', CHR)) |>
    dplyr::rename(P=Pvalue, EAF=AF)
  file <- sub('.tsv$', '.tsv.gz', file)
  vroom::vroom_write(rv, new_file_name)
})

# az_files <- Sys.glob('/local-scratch/data/ukb-seq/downloads/azexwas/*.csv.gz')
# az_dir <- '/local-scratch/data/hg38/azexwas'
# file.copy('/local-scratch/data/ukb-seq/downloads/azexwas/metadata.tsv', az_dir)
az_metadata <- vroom::vroom('/local-scratch/data/ukb-seq/downloads/azexwas/metadata.tsv', show_col_types = F)
az_metadata$trait <- sub('.*#', '', az_metadata$trait)

az_files <- Sys.glob('/Users/wt23152/Documents/Projects/scratch/021/data/study/rare/[1u]*')
az_dir <- '/Users/wt23152/Documents/Projects/scratch/021/data/study/rare/new/'
lapply(az_files, function(file) {
  new_basename <- basename(sub("\\..*$", "", file))
  new_basename <- sub('_buildGRCh38', '', new_basename)
  new_file_name <- paste0(genebass_dir, new_basename, '.tsv.gz')

  rv <- vroom::vroom(file)
  if ("Effect size" %in% colnames(rv)) {
    rv <- tidyr::separate(rv, col=Variant, into=c("CHR", "BP", "OA", "EA"), sep='-', remove = F) |>
      dplyr::rename(P='p-value', EAF='AAF', BETA='Effect size', SE='Effect size standard error') |>
      dplyr::mutate(Phenotype = sub(".*#", "", Phenotype))
  } else if ("Odds ratio" %in% colnames(rv)) {
    rv <- tidyr::separate(rv, col=Variant, into=c("CHR", "BP", "OA", "EA"), sep='-', remove = F) |>
      dplyr::rename(P='p-value', EAF='Control AAF', OR='Odds ratio', CI_LOWER='Odds ratio LCI', CI_UPPER='Odds ratio UCI') |>
      dplyr::mutate(Phenotype = sub(".*#", "", Phenotype))
  }
  vroom::vroom_write(rv, new_file_name)
})
