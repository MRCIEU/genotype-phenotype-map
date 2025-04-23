options(error = function() traceback(20))
Sys.setenv('VROOM_CONNECTION_SIZE' = 500000)
data_dir <- Sys.getenv('DATA_DIR')
results_dir <- Sys.getenv('RESULTS_DIR')
oracle_server <- Sys.getenv('ORACLE_SERVER')
TEST_RUN <- Sys.getenv('TEST_RUN', NA)

genome_wide_p_value_threshold <- 5e-8
lowest_p_value_threshold <- 1.5e-4
lowest_rare_p_value_threshold <- 1.5e-4

gpm_website_data <- list(
  url = 'https://gpm.opengwas.io',
  contact = 'https://gpm.opengwas.io/contact',
  name = 'The Genotype-Phenotype Map Team'
)

latest_results_dir <- glue::glue('{results_dir}latest/')
pipeline_metadata_dir <- glue::glue('{data_dir}pipeline_metadata/')
ld_block_data_dir <- glue::glue('{data_dir}ld_blocks/')
ld_reference_panel_dir <- glue::glue('{data_dir}ld_reference_panel_hg38/')
liftover_dir <- glue::glue('{data_dir}liftover/')
extracted_study_dir <- glue::glue('{data_dir}study/')
variant_annotation_dir <- glue::glue('{data_dir}variant_annotation/')

server_sync_dir <- file.path(data_dir, 'rsync_to_server')
oracle_data_dir <- '/oradiskvdb1/data/'

bespoke_parsing_options <- list(none='none', gtex_sqtl='gtex_sqtl')

#This is an intentionally ordered list
#splice_variant -> transcript -> gene_expression -> protein -> metabolome -> phenotype... methylation goes where?
data_types <- list(splice_variant='splice_variant',
                           transcript='transcript',
                           gene_expression='gene_expression',
                           protein='protein',
                           metabolome='metabolome',
                           phenotype='phenotype'
)
data_type_names <- list(splice_variant='sQTL',
                           transcript='tQTL',
                           gene_expression='eQTL',
                           protein='pQTL',
                           metabolome='mQTL',
                           phenotype='GWAS'
)
study_categories <- list(binary='binary', continuous='continuous', categorical='categorical')
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

standardised_gwas_columns <- c('CHR','BP','EA','OA','EAF','BETA','SE','P','SNP','Z', 'GENE')
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
  cis_trans = vroom::col_character()
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

update_directories_for_worker <- function(worker_guid) {
  ld_block_data_dir <<- glue::glue('{data_dir}ld_blocks/gwas_upload/{worker_guid}/')
  extracted_study_dir <<- glue::glue('{data_dir}study/gwas_upload/{worker_guid}/')
  pipeline_metadata_dir <<- glue::glue('{data_dir}pipeline_metadata/gwas_upload/{worker_guid}/')
}