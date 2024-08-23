options(error = function() traceback(20))
Sys.setenv('VROOM_CONNECTION_SIZE' = 500000)
data_dir <- Sys.getenv('DATA_DIR')
results_dir <- Sys.getenv('RESULTS_DIR')
TEST_RUN <- Sys.getenv('TEST_RUN', NA)

DEFAULT_P_VALUE_THRESHOLD <- 5e-8

MINIMUM_STUDY_REGION_SIZE <- 200

pipeline_metadata_dir <- paste0(data_dir, 'pipeline_metadata/')
ld_block_data_dir <- paste0(data_dir, 'ld_blocks/')
ld_block_matrices_dir <- paste0(data_dir, 'ld_block_matrices/allele_flip/')
thousand_genomes_dir <- paste0(data_dir, '1000genomes/')
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
data_formats <- list(opengwas='opengwas', besd='besd')
cis_trans <- list(cis_only='cis', trans_only='trans', cis_trans='cis_and_trans')
ancestry_map <- list(EUR='European', EAS='East Asian', AFR='African')
reverse_ancestry_map <- setNames(names(ancestry_map), ancestry_map)

reference_builds <- list(GRCh36="GRCh36", GRCh37="GRCh37", GRCh38="GRCh38")
available_liftover_conversions <- list(
  'GRCh36GRCh37' = 'hg18ToHg19.over.chain.gz',
  'GRCh38GRCh37' = 'hg18ToHg19.over.chain.gz',
  'GRCh37GRCh38' = 'hg18ToHg19.over.chain.gz',
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
                        ld_matrix_prefix=paste0(ld_block_matrices_dir, block)
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

standardise_alleles <- function(gwas) {
  gwas$EA <- toupper(gwas$EA)
  gwas$OA <- toupper(gwas$OA)
  gwas$EAF <- as.numeric(gwas$EAF)

  to_flip <- (gwas$EA > gwas$OA) & (!gwas$EA %in% c("D", "I"))
  if (any(to_flip)) {
    gwas$EAF[to_flip] <- 1 - gwas$EAF[to_flip]
    if ('BETA' %in% names(gwas)) {
      gwas$BETA[to_flip] <- -1 * gwas$BETA[to_flip]
    }
    if ('Z' %in% names(gwas)) {
      gwas$Z[to_flip] <- -1 * gwas$Z[to_flip]
    }

    temp <- gwas$OA[to_flip]
    gwas$OA[to_flip] <- gwas$EA[to_flip]
    gwas$EA[to_flip] <- temp
  }

  gwas$SNP <- toupper(paste0(gwas$CHR, ":", format(gwas$BP, scientific = F, trim = T), "_", gwas$EA, "_", gwas$OA))
  return(gwas)
}
