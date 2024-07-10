options(error = function() traceback(20))
Sys.setenv('VROOM_CONNECTION_SIZE' = 500000)
data_dir <- Sys.getenv('DATA_DIR')
results_dir <- Sys.getenv('RESULTS_DIR')

pipeline_metadata_dir <- paste0(data_dir, 'pipeline_metadata/')
ld_block_data_dir <- paste0(data_dir, 'ld_blocks/')
ld_block_matrices_dir <- paste0(data_dir, 'ld_block_matrices/')
extracted_study_dir <- paste0(data_dir, 'study/')

ld_block_results_dir <- paste0(results_dir, 'ld_blocks/')

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
databases <- list(opengwas='opengwas')
data_source <- list(ukb='UK Biobank', gtex='GTEx', finngen='Finn Gen', bbj='Biobank Japan', eqtl_gen='eQTL Gen', ukb_ppp='UK Biobank Proteomic...')
ancestry_map <- list(EUR='European', EAS='East Asian', AFR='African')
reverse_ancestry_map <- setNames(names(ancestry_map), ancestry_map)

file_prefix <- function(file_path) {
  file_name <- basename(file_path)
  file_prefix <- sub('\\..*', '', file_name)
  return(file_prefix)
}
