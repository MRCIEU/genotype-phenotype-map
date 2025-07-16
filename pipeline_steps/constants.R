options(error = function() traceback(20))
Sys.setenv('VROOM_CONNECTION_SIZE' = 500000)

data_dir <- Sys.getenv('DATA_DIR')
gwas_upload_dir <- Sys.getenv('GWAS_UPLOAD_DIR')
results_dir <- Sys.getenv('RESULTS_DIR')
oracle_server <- Sys.getenv('ORACLE_SERVER')
TEST_RUN <- Sys.getenv('TEST_RUN', NA)

genome_wide_p_value_threshold <- 5e-8
lowest_p_value_threshold <- 1.5e-4
lowest_rare_p_value_threshold <- 1.5e-4
minimum_extraction_size <- 150

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
                           phenotype='phenotype'
)
data_type_names <- list(splice_variant='sQTL',
                           transcript='tQTL',
                           gene_expression='eQTL',
                           protein='pQTL',
                           methylation='methQTL',
                           metabolome='metaQTL',
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

#' Generate log Bayes Factor from Z-score
#'
#' @param z Z-score
#' @param se Standard error
#' @param prior_v Prior variance
#'
#' @return Log Bayes Factor
convert_z_to_lbf <- function(z, se, eaf, sample_size, study_type, effect_priors=c(continuous=0.15, categorical=0.2)) {
  estimated_sd <- estimate_variance(se, eaf, sample_size)
  if (study_type == study_categories$continuous) {
    sd_prior <- effect_priors[study_categories$continuous] * estimated_sd
  } else {
    sd_prior <- effect_priors[study_categories$categorical]
  }
  r <- sd_prior^2 / (sd_prior^2 + se^2)
  lbf = (log(1 - r) + (r * z^2)) / 2
  return(lbf)
}

##' Estimate trait standard deviation given vectors of variance of coefficients,  MAF and sample size
##'
##' Estimate is based on var(beta-hat) = var(Y) / (n * var(X))
##' var(X) = 2*maf*(1-maf)
##' so we can estimate var(Y) by regressing n*var(X) against 1/var(beta)
##' 
##' @title Estimate trait variance, internal function
##' @param SE vector of standard errors
##' @param EAF vector of MAF (same length as SE)
##' @param n sample size
##' @return estimated standard deviation of Y
##' 
estimate_variance <- function(se, eaf, n) {
    oneover <- 1/se^2
    nvx <- 2 * n * eaf * (1-eaf)
    m <- lm(nvx ~ oneover - 1)
    cf <- coef(m)[['oneover']]
    if(cf < 0)
        stop("estimated sdY is negative - this can happen with small datasets, or those with errors.  A reasonable estimate of sdY is required to continue.")
    return(sqrt(cf))
}

#' Convert log Bayes Factor to abs(Z-score)
#'
#' @param lbf Log Bayes Factor
#' @param se Standard error
#' @param prior_v Prior variance
#'
#' @return abs(Z-score): Note that this is the absolute value of the Z-score, not the Z-score itself
convert_lbf_to_abs_z <- function(lbf, se, prior_v = 50) {
  r <- prior_v / (prior_v + se^2)
  z <- sqrt((2 * lbf - log(sqrt(1-r)))/r)
  return(z)
}

#' Convert log Bayes Factor to summary stats
#'
#' @param gwas of summary statistics, with EAF as a mandatory column (allele frequencies for each SNP)
#' @param lbf p-vector of log Bayes Factors for each SNP
#' @param n Overall sample size
#' @param prior_v Variance of prior distribution. SuSiE uses 50
#'
#' @return tibble with altered BETA, SE, P, and Z
update_gwas_with_log_bayes_factor <- function(gwas, lbf, sample_size, prior_v = 50) {
  se <- sqrt(1 / (2 * sample_size * gwas$EAF * (1-gwas$EAF)))
  r <- prior_v / (prior_v + se^2)
  z <- sqrt((2 * lbf - log(sqrt(1-r)))/r)
  beta <- z * se
  p <- abs(2 * pnorm(abs(z), lower.tail = F))

  gwas <- dplyr::mutate(gwas, BETA = beta, SE = se, P = p, Z = z) |>
    dplyr::filter(!is.na(BETA) & !is.na(SE) & BETA != Inf & SE != Inf)

  return(gwas)
}