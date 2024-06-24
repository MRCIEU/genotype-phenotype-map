source("constants.R")

parser <- argparser::arg_parser("Finemap studies per region")
parser <- argparser::add_argument(parser, "--ld_region_prefix", help = "GWAS filename", type = "character")
parser <- argparser::add_argument(parser, "--ld_block_dir", help = "LD block that the ", type = "character")
args <- argparser::parse_args(parser)

# args <- list(ld_block_dir="/Users/wt23152/Documents/Projects/scratch/011/data/ld_blocks/EUR/10/93335047_95396367/",
#              ld_region_prefix="/Users/wt23152/Documents/Projects/scratch/011/data/ld_block_matrices/EUR/10_93335047_95396367"
# )
#
# args <- list(ld_block_dir=paste0(data_dir, "/ld_blocks/EUR/10/93335047_95396367/"),
#              ld_region_prefix=paste0(data_dir, "/ld_block_matrices/EUR/10_93335047_95396367")
# )


main <- function(args) {
  ld_region_finemap_dir <- paste0(args$ld_block_dir, '/finemapped/')
  ld_region <- vroom::vroom(paste0(args$ld_region_prefix, '.ld'), col_names=F, show_col_types = F)
  ld_region <- ld_region[, 1:(ncol(ld_region)-1)]
  ld_region_from_reference_panel <- vroom::vroom(paste0(args$ld_region_prefix, '.tsv'), show_col_types = F)

  extracted_study_file <- paste0(args$ld_block_dir, '/extracted_studies.tsv')
  extracted_studies <- vroom::vroom(extracted_study_file, show_col_types = F)

  finemapping_results <- apply(extracted_studies, 1, function (study) {
    finemap_info <- data.frame(finemap_successful=F, number_of_finemap_gwases=1)
    finemap_file_prefix <- sub('imputed', 'finemapped', study[['file']]) |>
      stringr::str_sub(end=-5)
    gwas <- vroom::vroom(study[['file']], show_col_types = F)

    if (file.exists(paste0(finemap_file_prefix, "_1.tsv"))) {
      print('Already finemapped, skipping.')
      finemap_info$finemap_successful <- study[['finemap_successful']]
      finemap_info$number_of_finemap_gwases <- study[['number_of_finemap_gwases']]
      return(finemap_info)
    }

    if (typeof(gwas$EAF) == 'character') {
      warning('EAF is not populated, cant split results.  Skipping.')
      finemap_file <- paste0(finemap_file_prefix, '_1.tsv')
      file.copy(study[['file']], finemap_file)
      file.symlink(finemap_file, ld_region_finemap_dir, overwrite=T)
      return(finemap_info)
    }

    gwas <- dplyr::filter(gwas, RSID %in% ld_region_from_reference_panel$RSID)
    keep <- ld_region_from_reference_panel$RSID %in% gwas$RSID
    ld_for_gwas <- ld_region[keep, keep]
    ld_matrix <- matrix(as.vector(data.matrix(ld_for_gwas)), nrow=nrow(ld_for_gwas), ncol=nrow(ld_for_gwas))

    tryCatch(expr = {
      # carma_results <- carma(study, gwas, ld_for_gwas)
      conditioned_gwases <- susie(study[['file']], gwas, ld_matrix, as.numeric(study[['sample_size']]))

      if (length(conditioned_gwases) > 1) {
        #maybe clump and find new lead snp to tag it with?
        #create finemapped_studies.tsv
      }

      for (i in seq_along(conditioned_gwases)) {
        finemap_file <- paste0(finemap_file_prefix, '_', i,'.tsv')
        vroom::vroom_write(conditioned_gwases[[i]], finemap_file)
        file.symlink(finemap_file, ld_region_finemap_dir, overwrite=T)
      }
      finemap_info$finemap_succesful<- T
      finemap_info$number_of_finemap_gwases <- length(conditioned_gwases)
      return(finemap_info)
    }, warning = function(w) {
      print(w)
      finemap_file <- paste0(finemap_file_prefix, '_1.tsv')
      file.copy(study[['file']], finemap_file)
      file.symlink(finemap_file, ld_region_finemap_dir, overwrite=T)
      return(finemap_info)
    })
  }) |> dplyr::bind_rows()

  extracted_studies <- dplyr::mutate(extracted_studies,
    finemap_successful = finemapping_results$finemap_successful,
    number_of_finemap_gwases = finemapping_results$number_of_finemap_gwases
  )
  vroom::vroom_write(extracted_studies, extracted_study_file)
}

#lapply(all_conditioned_gwases, plot_gwas)
#TODO: delete later
plot_gwas <- function(gwas) {
  range <- c(min(gwas$position), max(gwas$position))
  if(!"P" %in% colnames(gwas)) return()
  gwas$P[gwas$P == 0] <- .Machine$double.xmin
  gwas <- gwas[!is.na(gwas$P), ]
  qqman::manhattan(gwas, p='P', chr = 'CHR', bp = 'BP', snp = 'RSID', xlim=range)
}

# flatten_list_of_lists <- function(x) {
#   if (is.data.frame(x)) return(list(x))
#   if (!is.list(x)) return(x)
#   unlist(lapply(x, flattenMixed), FALSE)
# }
#
# convert_negative_log_p_to_p <- function(gwas) {
#   gwas$P <- 10^-gwas$LOG_P
#   return(gwas)
# }

susie <- function(study, gwas, ld_matrix, n) {
  start_time <- Sys.time()
  susie_results <- susieR::susie_rss(z=gwas$Z, R=ld_matrix, n=n, L=5)
  print(Sys.time() - start_time)

  saveRDS(susie_results, paste0(data_dir, 'finemap_tests/susie_', basename(study), '.rds'))

  conditioned_gwases <- apply(susie_results$lbf_variable, 1, function(lbf) {
    conditioned_gwas <- lbf_to_z_cont(gwas$RSID, lbf, n, gwas$EAF)
    gwas <- dplyr::mutate(gwas, BETA = conditioned_gwas$BETA, SE = conditioned_gwas$SE, P = conditioned_gwas$P, Z = conditioned_gwas$Z)
    return(gwas)
  })
  return(conditioned_gwases)
}

#' Convert log Bayes Factor to summary stats
#'
#' @param lbf p-vector of log Bayes Factors for each SNP
#' @param n Overall sample size
#' @param af p-vector of allele frequencies for each SNP
#' @param prior_v Variance of prior distribution. SuSiE uses 50
#'
#' @return tibble with lbf, af, beta, se, z
#' @export
lbf_to_z_cont <- function(RSID, lbf, n, af, prior_v = 50) {
  SE <- sqrt(1 / (2 * n * af * (1-af)))
  r <- prior_v / (prior_v + SE^2)
  Z <- sqrt((2 * lbf - log(sqrt(1-r)))/r)
  BETA <- Z * SE
  P <- abs(2 * pnorm(abs(Z), lower.tail = F))
  return(data.frame(RSID, BETA, SE, Z, P))
}

carma <- function(study, gwas, ld_matrix) {
  #input_columns <- list(BETA="beta", SE="se", RSID="rsid", CHR="chromosome", BP="position", EA="allele1", OA="allele2", EAF="maf", LOG_P="lp")
  #gwas <- standardise_gwas(gwas, input_columns=input_columns, N=args$gwas_n)
  gwas <- convert_beta_and_se_to_z_score(gwas)
  lambda_list<-list()
  z.list<-list()
  ld.list<-list()
  lambda.list<-list()
  z.list[[1]] <- gwas$z
  ld.list[[1]] <- ld_matrix
  lambda.list[[1]] <- 1
  print(dim(ld_matrix))
  print(typeof(ld_matrix))

  #outlier_switch=T because of out sample LD matrix
  carma_result <- CARMA::CARMA(z.list, ld.list, lambda.list=lambda.list, outlier.switch=TRUE)
  print(paste0(data_dir, '/finemap_tests/carma_', basename(study), '.rds'))
  saveRDS(carma_result, paste0(data_dir, '/finemap_tests/carma_', basename(study), '.rds'))
}

main(args)
