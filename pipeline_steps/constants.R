options(error = function() traceback(20))
Sys.setenv('VROOM_CONNECTION_SIZE' = 500000)

data_dir <- Sys.getenv('DATA_DIR')
gwas_upload_dir <- Sys.getenv('GWAS_UPLOAD_DIR')
results_dir <- Sys.getenv('RESULTS_DIR')
oracle_api_server <- Sys.getenv('ORACLE_SERVER')
oracle_pipeline_server <- Sys.getenv('ORACLE_SERVER')
TEST_RUN <- Sys.getenv('TEST_RUN', NA)

genome_wide_p_value_threshold <- 5e-8
lowest_p_value_threshold <- 1.5e-4
minimum_extraction_size <- 150
posterior_prob_h4_threshold <- 0.8

#TODO: look into removing these
posterior_prob_threshold <- 0.5
posterior_prob_thresholds <- list(
  strong=0.8,
  moderate=0.6
)

gpm_website_data <- list(
  url = 'https://gpm.opengwas.io',
  contact = 'https://gpm.opengwas.io/contact',
  name = 'The Genotype-Phenotype Map Team'
)

latest_results_dir <- glue::glue('{results_dir}latest/')
current_results_dir <- glue::glue('{results_dir}current/')
results_analysis_dir <- glue::glue('{latest_results_dir}analysis/')
pipeline_metadata_dir <- glue::glue('{data_dir}pipeline_metadata/')
ld_block_data_dir <- glue::glue('{data_dir}ld_blocks/')
ld_reference_panel_dir <- glue::glue('{data_dir}ld_reference_panel_hg38/')
liftover_dir <- glue::glue('{data_dir}liftover/')
extracted_study_dir <- glue::glue('{data_dir}study/')
variant_annotation_dir <- glue::glue('{data_dir}variant_annotation/')
svg_dir <- glue::glue('{data_dir}svgs/')

server_sync_dir <- file.path(data_dir, 'rsync_to_server')
oracle_data_dir <- '/oradiskvdb1/data/'

bespoke_parsing_options <- list(none='none', gtex_sqtl='gtex_sqtl')

data_types <- list(splice_variant='splice_variant',
                           transcript='transcript',
                           gene_expression='gene_expression',
                           protein='protein',
                           methylation='methylation',
                           metabolome='metabolome',
                           cell_trait='cell_trait',
                           plasma_protein='plasma_protein',
                           phenotype='phenotype'
)
data_type_names <- list(splice_variant='sQTL',
                           transcript='tQTL',
                           gene_expression='eQTL',
                           protein='pQTL',
                           methylation='methQTL',
                           metabolome='metaQTL',
                           cell_trait='Cell Trait',
                           plasma_protein='Plasma Protein',
                           phenotype='Phenotype'
)
study_categories <- list(continuous='continuous', categorical='categorical')
data_formats <- list(opengwas='opengwas', besd='besd', tsv='tsv')
cis_trans <- list(cis_only='cis', trans_only='trans', cis_trans='cis_trans')
variant_types <- list(common='common', rare_exome='rare_exome', rare_wgs='rare_wgs')
ancestry_map <- list(EUR='European', EAS='East Asian', AFR='African', SAS='South Asian')
reverse_ancestry_map <- setNames(names(ancestry_map), ancestry_map)

reference_builds <- list(GRCh37="GRCh37", GRCh38="GRCh38")
available_liftover_conversions <- list(
  'GRCh36GRCh37' = glue::glue('{liftover_dir}hg18ToHg19.over.chain.gz'),
  'GRCh38GRCh37' = glue::glue('{liftover_dir}hg38ToHg19.over.chain.gz'),
  'GRCh37GRCh38' = glue::glue('{liftover_dir}/hg19ToHg38.over.chain.gz')
)
extraction_file_types <- list(vcf='vcf', csv='csv')
coverage_types <- list(dense='dense', sparse='sparse')

standardised_gwas_columns <- c('CHR','BP','EA','OA','EAF','BETA','SE','P','SNP','Z','GENE')
required_columns <- c("CHR","BP","EA","OA","EAF","BETA","SE","P")

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
  time_taken = vroom::col_character(),
  cis_trans = vroom::col_character(),
  ignore = vroom::col_logical()
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
  # ld_info <- dplyr::bind_cols(ld_info, ld_block_components(block))
  return(ld_info)
}

construct_ld_block <- function(ancestry, chr, start, stop) {
  block <- ld_block_string(ancestry, chr, start, stop)
  ld_info <- ld_block_dirs(block)
  ld_info <- dplyr::bind_cols(ld_info, data.frame(ancestry = ancestry, chr = chr, start = start, stop = stop))
  return(ld_info)
}

ld_block_string <- function(ancestry, chr, start, stop) {
  return(glue::glue('{ancestry}/{chr}/{start}-{stop}'))
}

ld_block_components <- function(ld_block) {
  components <- strsplit(ld_block, '[/-]')[[1]]
  return(data.frame(
    ancestry = components[1], 
    chr = components[2], 
    start = as.numeric(components[3]), 
    stop = as.numeric(components[4])
  ))
}

flattened_ld_block_name <- function(ld_block_string) {
  return(gsub('[/-]', '_', ld_block_string))
}

update_directories_for_worker <- function(worker_guid) {
  ld_block_data_dir <<- glue::glue('{gwas_upload_dir}ld_blocks/gwas_upload/{worker_guid}/')
  extracted_study_dir <<- glue::glue('{gwas_upload_dir}study/gwas_upload/{worker_guid}/')
  pipeline_metadata_dir <<- glue::glue('{gwas_upload_dir}pipeline_metadata/gwas_upload/{worker_guid}/')
}

diff_time_taken <- function(start_time) {
  return(hms::as_hms(difftime(Sys.time(), start_time)))
}

safe_lapply <- function(X, FUN, ...) {
  lapply(X, function(x) {
    tryCatch(
      FUN(x, ...),
      error = function(e) {
        stop("An error occurred: ", conditionMessage(e))
      }
    )
  })
}
