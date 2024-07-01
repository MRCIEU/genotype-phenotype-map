source("constants.R")

parser <- argparser::arg_parser("Finemap studies per region")
parser <- argparser::add_argument(parser, "--ld_region_prefix", help = "GWAS filename", type = "character")
parser <- argparser::add_argument(parser, "--ld_block_dir", help = "LD block that the ", type = "character")
args <- argparser::parse_args(parser)

main <- function(args) {
  ld_region_finemap_dir <- paste0(args$ld_block_dir, '/finemapped/')
  dir.create(ld_region_finemap_dir, recursive=T, showWarnings=F)

  ld_region <- vroom::vroom(paste0(args$ld_region_prefix, '.ld'), col_names=F, show_col_types = F)
  ld_region <- ld_region[, 1:(ncol(ld_region)-1)]
  ld_region_from_reference_panel <- vroom::vroom(paste0(args$ld_region_prefix, '.tsv'), show_col_types = F)

  imputed_studies_file <- paste0(args$ld_block_dir, '/imputed_studies.tsv')
  imputed_studies <- vroom::vroom(imputed_studies_file, show_col_types = F)

  finemapped_results <- apply(imputed_studies, 1, function (study) {
    sample_size <- as.numeric(study['sample_size'])
    print(paste('Finemapping', study['file']))
    finemap_file_prefix <- sub('imputed', 'finemapped', study[['file']])
    finemap_file_prefix <- sub("\\..*", "", finemap_file_prefix)
    dir.create(dirname(finemap_file_prefix), recursive=T, showWarnings=F)

    failed_finemap_file <- paste0(finemap_file_prefix, '_1.tsv')
    failed_finemap_symlink <- paste0(ld_region_finemap_dir, study['study'], "_", study['chr'], "_", study['bp'], "_1.tsv")
    failed_finemap_info <- data.frame(study=study[['study']], file=failed_finemap_file, chr=study[['chr']],
        bp=study['bp'], p_value_threshold=study['p_value_threshold'], category=study['category'], sample_size=sample_size, finemap_suceeded=F
    )

    gwas <- vroom::vroom(study[['file']], show_col_types = F)

    if (file.exists(paste0(finemap_file_prefix, "_1.tsv"))) {
      print('Already finemapped, skipping.')
      return()
    }

    if (typeof(gwas$EAF) == 'character') {
      warning('EAF is not populated, cant split results, skipping.')
      file.copy(study[['file']], failed_finemap_file)
      file.symlink(failed_finemap_file, failed_finemap_symlink)
      return(failed_finemap_info)
    }


    bp <- as.numeric(study['bp'])
    range <- 1000000
    gwas <- dplyr::filter(gwas, BP > (bp - floor(range/2)) & BP < (bp + floor(range/2)))
    
    print(nrow(gwas))
    gwas <- dplyr::filter(gwas, RSID %in% ld_region_from_reference_panel$RSID)
    print(nrow(gwas))
    keep <- ld_region_from_reference_panel$RSID %in% gwas$RSID
    ld_for_gwas <- ld_region[keep, keep]
    ld_matrix <- matrix(as.vector(data.matrix(ld_for_gwas)), nrow=nrow(ld_for_gwas), ncol=ncol(ld_for_gwas))

<<<<<<< Updated upstream
    start_time <- Sys.time()
    carma_result <- carma(study, gwas, ld_for_gwas)
    print(paste('carma time', Sys.time() - start_time))
    saveRDS(carma_result, paste0(data_dir, 'finemap_tests/carma_small_', file_prefix(study[['file']]), '.rds'))
=======

    tryCatch(expr = {
      # carma_results <- carma(study, gwas, ld_for_gwas)
      conditioned_gwases <- susie(study[['file']], gwas, ld_matrix, as.numeric(study[['sample_size']]))
>>>>>>> Stashed changes

    start_time <- Sys.time()
    susie_result <- susieR::susie_rss(z=gwas$Z, R=ld_matrix, n=sample_size, L=5)
    print(paste('susie time', Sys.time() - start_time))
    saveRDS(susie_result, paste0(data_dir, 'finemap_tests/susie_small_', file_prefix(study[['file']]), '.rds'))

    if (susie_result$converged == F) {
      print(paste('susie result failed to converge for', study[['file']], ', skipping..'))
      file.copy(study[['file']], failed_finemap_file)
      file.symlink(failed_finemap_file, failed_finemap_symlink)
      return(failed_finemap_info)
    }

    #TODO: do some sort of filtering based on susie output if the results are above certain thresholds

    conditioned_gwases <- apply(susie_result$lbf_variable, 1, function(lbf) {
      conditioned_gwas <- lbf_to_z_cont(gwas$RSID, lbf, sample_size, gwas$EAF)
      gwas <- dplyr::mutate(gwas, BETA = conditioned_gwas$BETA, SE = conditioned_gwas$SE, P = conditioned_gwas$P, Z = conditioned_gwas$Z)
      return(gwas)
    })

    new_bps <- c()
    new_files <- c()
    for (i in seq_along(conditioned_gwases)) {
      finemap_file <- paste0(finemap_file_prefix, '_', i, '.tsv')
      finemap_symlink <- paste0(ld_region_finemap_dir, study['study'], "_", study['chr'], "_", study['bp'], "_", i, ".tsv")
      vroom::vroom_write(conditioned_gwases[[i]], finemap_file)
      file.symlink(finemap_file, finemap_symlink)

      new_bps <- c(new_bps, find_new_top_snp(finemap_file, study[['bp']]))
      new_files <- c(new_files, finemap_file)
    }
    succeeded_finemap_info <- data.frame(study=study[['study']], file=new_files, chr=study[['chr']],
      bp=new_bps, p_value_threshold=study['p_value_threshold'], category=study['category'], sample_size=sample_size, finemap_suceeded=T
    )
    return(succeeded_finemap_info)
  }) |> dplyr::bind_rows()

  finemapped_results_file <- paste0(args$ld_block_dir, '/finemapped_studies.tsv')
  vroom::vroom_write(finemapped_results, finemapped_results_file, append = file.exists(finemapped_results_file))
  vroom::vroom_write(data.frame(), file=paste0(args$ld_block_dir, '/finemapping_complete'))
}

#lapply(all_conditioned_gwases, plot_gwas)
#TODO: delete later
plot_finemapped_gwas <- function(gwas) {
  range <- c(min(gwas$BP), max(gwas$BP))
  if(!"P" %in% colnames(gwas)) return()
  gwas$P[gwas$P == 0] <- .Machine$double.xmin
  gwas <- gwas[!is.na(gwas$P), ]
  qqman::manhattan(gwas, p='P', chr = 'CHR', bp = 'BP', snp = 'RSID', xlim=range)
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

#TODO: maybe use plink --clump to find it?  And filter if the top hit isn't 5e-8?
find_new_top_snp <- function(gwas_file, bp) {
  return(bp)
}

carma <- function(study, gwas, ld_matrix) {
  lambda_list<-list()
  z.list<-list()
  ld.list<-list()
  lambda.list<-list()
  z.list[[1]] <- gwas$Z
  ld.list[[1]] <- ld_matrix
  lambda.list[[1]] <- 1

  #outlier_switch=T because of out sample LD matrix
  carma_result <- CARMA::CARMA(z.list, ld.list, lambda.list=lambda.list, outlier.switch=TRUE)
  return(carma_result)
}

main(args)
