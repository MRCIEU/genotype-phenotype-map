options(error = function() traceback(20))
Sys.setenv('VROOM_CONNECTION_SIZE' = 500000)
data_dir <- Sys.getenv('DATA_DIR')
results_dir <- Sys.getenv('RESULTS_DIR')
TEST_RUN <- Sys.getenv('TEST_RUN', NA)

genome_wide_p_value_threshold <- 5e-8
lowest_p_value_threshold <- 1.5e-4

pipeline_metadata_dir <- glue::glue('{data_dir}pipeline_metadata/')
ld_block_data_dir <- glue::glue('{data_dir}ld_blocks/')
ld_reference_panel_dir <- glue::glue('{data_dir}ld_reference_panel_hg38/')
liftover_dir <- glue::glue('{data_dir}liftover/')
extracted_study_dir <- glue::glue('{data_dir}study/')

bespoke_parsing_options <- list(none='none', gtex_sqtl='gtex_sqtl')

#This is an intentionally ordered list
#splice_variant -> transcript -> gene_expression -> protein -> metabolome -> phenotype... methylation goes where?
ordered_data_types <- list(splice_variant='splice_variant',
                           transcript='transcript',
                           gene_expression='gene_expression',
                           protein='protein',
                           metabolome='metabolome',
                           phenotype='phenotype'
)
study_categories <- list(binary='Binary', continuous='Continuous')
data_formats <- list(opengwas='opengwas', besd='besd', hail='hail')
cis_trans <- list(cis_only='cis', trans_only='trans', cis_trans='cis_and_trans')
variant_type <- list(common='common', rare='rare')
ancestry_map <- list(EUR='European', EAS='East Asian', AFR='African', SAS='South Asian')
reverse_ancestry_map <- setNames(names(ancestry_map), ancestry_map)

reference_builds <- list(GRCh36="GRCh36", GRCh37="GRCh37", GRCh38="GRCh38")
available_liftover_conversions <- list(
  'GRCh36GRCh37' = glue::glue('{liftover_dir}hg18ToHg19.over.chain.gz'),
  'GRCh38GRCh37' = glue::glue('{liftover_dir}hg38ToHg19.over.chain.gz'),
  'GRCh37GRCh38' = glue::glue('{liftover_dir}/hg19ToHg38.over.chain.gz')
)

standardised_column_types <- vroom::cols(
  chr = vroom::col_character(),
  bp = vroom::col_number(),
  sample_size = vroom::col_number(),
  p_value_threshold = vroom::col_number(),
  snps_removed_by_reference_panel=vroom::col_number(),
  eaf_from_reference_panel=vroom::col_logical(),
  time_taken=vroom::col_character()
)

imputed_column_types <- vroom::cols(
  chr = vroom::col_character(),
  bp = vroom::col_number(),
  sample_size = vroom::col_number(),
  p_value_threshold = vroom::col_number(),
  time_taken=vroom::col_character()
)

finemapped_column_types <- vroom::cols(
  chr = vroom::col_character(),
  bp = vroom::col_number(),
  min_p = vroom::col_number(),
  sample_size = vroom::col_number(),
  p_value_threshold = vroom::col_number(),
  first_finemap_num_results = vroom::col_number(),
  second_finemap_num_results = vroom::col_number(),
  qc_step_run = vroom::col_logical(),
  snps_removed_by_qc = vroom::col_number(),
  time_taken = vroom::col_character()
)

file_prefix <- function(file_path) {
  file_name <- basename(file_path)
  file_prefix <- sub('\\..*', '', file_name)
  return(file_prefix)
}

ld_block_dirs <- function(block) {
  ld_info <- data.frame(block = block,
                        ld_block_data = glue::glue('{ld_block_data_dir}{block}'),
                        ld_reference_panel_prefix=glue::glue('{ld_reference_panel_dir}{block}')
  )
  return(ld_info)
}

construct_ld_block <- function(ancestry, chr, start, stop) {
  block <- ld_block_string(ancestry, chr, start, stop)
  return(ld_block_dirs(block))
}

ld_block_string <- function(ancestry, chr, start, stop) {
  return(glue::glue('{ancestry}/{chr}/{start}-{stop}'))
}

gwas_health_check <- function(gwas) {
  if (any(gwas$P < 0 | gwas$P > 1)) {
    stop("GWAS has some P values outside accepted range.  Please fix GWAS or remove it from pipeline")
  }
  if (any(gwas$EAF < 0 | gwas$EAF > 1)) {
    stop("GWAS has some EAF values outside accepted range.  Please fix GWAS or remove it from pipeline")
  }
  if (any(gwas$SE < 0) {
    stop("GWAS has some EAF values outside accepted range.  Please fix GWAS or remove it from pipeline")
  }
}

filter_gwas <- function(gwas, common=T) {
  gwas <- dplyr::filter(gwas,
    (is.na(EAF) | (EAF < 0.99 & EAF > 0.01)) &
    !is.na(CHR) & !is.na(CHR) & !is.na(BP) & !is.na(EA) & !is.na(OA) &
    !is.na(P) & !is.na(BETA)
  )
  return(gwas)
}