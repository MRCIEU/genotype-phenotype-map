options(error = function() traceback(20))
Sys.setenv('VROOM_CONNECTION_SIZE' = 500000)
data_dir <- Sys.getenv('DATA_DIR')
results_dir <- Sys.getenv('RESULTS_DIR')
TEST_RUN <- Sys.getenv('TEST_RUN', NA)

genome_wide_p_value_threshold <- 5e-8
lowest_p_value_threshold <- 1.5e-4

pipeline_metadata_dir <- paste0(data_dir, 'pipeline_metadata/')
ld_block_data_dir <- paste0(data_dir, 'ld_blocks/')
ld_reference_panel_dir <- paste0(data_dir, 'ld_reference_panel/')
liftover_dir <- paste0(data_dir, 'liftover/')
extracted_study_dir <- paste0(data_dir, 'study/')

ld_block_results_dir <- paste0(results_dir, 'ld_blocks/')
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

file_prefix <- function(file_path) {
  file_name <- basename(file_path)
  file_prefix <- sub('\\..*', '', file_name)
  return(file_prefix)
}

ld_block_dirs <- function(block) {
  ld_info <- data.frame(block = block,
                        ld_block_data = paste0(ld_block_data_dir, block),
                        ld_block_results = paste0(ld_block_results_dir, block),
                        ld_reference_panel_prefix=paste0(ld_reference_panel_dir, block)
  )
  return(ld_info)
}

construct_ld_block <- function(ancestry, chr, start, stop) {
  block <- ld_block_string(ancestry, chr, start, stop)
  return(ld_block_dirs(block))
}

ld_block_string <- function(ancestry, chr, start, stop) {
  return(paste0(ancestry, '/', chr, '/', start, '_', stop))
}
