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
    finemap_file_prefix <- sub('imputed', 'finemapped', study[['file']])
    finemap_file_prefix <- sub("\\..*", "", finemap_file_prefix)
    if (file.exists(paste0(finemap_file_prefix, "_1.tsv"))) {
      return()
    }
    cat(paste0('Finemapping ', study['file'], ": "))

    dir.create(dirname(finemap_file_prefix), recursive=T, showWarnings=F)
    gwas <- vroom::vroom(study[['file']], show_col_types = F)

    if (typeof(gwas$EAF) == 'character') {
      message('EAF is not populated, cant split results, skipping.')
      failed_finemap_info <- process_unfinemapped_gwas(gwas)
      return(failed_finemap_info)
    }

    keep <- ld_region_from_reference_panel$RSID %in% gwas$RSID
    ld_for_gwas <- ld_region[keep, keep]
    ld_matrix <- matrix(as.vector(data.matrix(ld_for_gwas)), nrow=nrow(ld_for_gwas), ncol=ncol(ld_for_gwas))
    testthat::expect_true(nrow(gwas) == nrow(ld_for_gwas), 'gwas and ld matrix should match size')
    susie_result <- susieR::susie_rss(z=gwas$Z, R=ld_matrix, n=sample_size, L=5)

    if (susie_result$converged == F || is.null(susie_result$sets$cs_index) || length(susie_result$sets$cs_index) == 0) {
      message('susie either didnt converge or has no credible sets, skipping.')
      failed_finemap_info <- process_unfinemapped_gwas(gwas)
      return(failed_finemap_info)
    }

    new_bps <- c()
    new_files <- c()
    for (i in susie_result$sets$cs_index) {
      conditioned_gwas <- update_gwas_with_log_bayes_factor(gwas, susie_result$lbf_variable[i, ], sample_size)
      #plot_finemapped_gwas(conditioned_gwas)
      finemap_file <- paste0(finemap_file_prefix, '_', i, '.tsv')
      finemap_symlink <- paste0(ld_region_finemap_dir, study['study'], "_", study['chr'], "_", study['bp'], "_", i, ".tsv")
      vroom::vroom_write(conditioned_gwas, finemap_file)
      file.symlink(finemap_file, finemap_symlink)

      #find lead SNP in new credible set
      important_row <- susie_result$sets$cs[paste0('L', i)][[1]][[1]]
      new_bp <- conditioned_gwas[important_row, ]$BP
      new_bps <- c(new_bps, new_bp)
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

process_unfinemapped_gwas <- function(gwas) {
  failed_finemap_file <- paste0(finemap_file_prefix, '_1.tsv')
  failed_finemap_symlink <- paste0(ld_region_finemap_dir, study['study'], "_", study['chr'], "_", study['bp'], "_1.tsv")
  failed_finemap_info <- data.frame(study=study[['study']],
                                    file=failed_finemap_file,
                                    chr=study[['chr']],
                                    bp=study['bp'],
                                    p_value_threshold=study['p_value_threshold'],
                                    category=study['category'],
                                    sample_size=sample_size,
                                    finemap_suceeded=F
  )

  gwas <- populate_beta_with_known_z_scores(gwas)
  vroom::vroom_write(gwas, failed_finemap_file)
  file.symlink(failed_finemap_file, failed_finemap_symlink)
  return(failed_finemap_info)
}

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
#' @param gwas of summary statistics, with EAF as a mandatory column (allele frequencies for each SNP)
#' @param lbf p-vector of log Bayes Factors for each SNP
#' @param n Overall sample size
#' @param prior_v Variance of prior distribution. SuSiE uses 50
#'
#' @return tibble with altered BETA, SE, P, LP, and Z
update_gwas_with_log_bayes_factor <- function(gwas, lbf, sample_size, prior_v = 50) {
  SE <- sqrt(1 / (2 * sample_size * gwas$EAF * (1-gwas$EAF)))
  r <- prior_v / (prior_v + SE^2)
  Z <- sqrt((2 * lbf - log(sqrt(1-r)))/r)
  BETA <- Z * SE
  P <- abs(2 * pnorm(abs(Z), lower.tail = F))
  LP <- -log10(P)

  gwas <- dplyr::mutate(gwas, BETA = BETA, SE = SE, P = P, LP = LP, Z = Z)
  return(gwas)
}

#' Calculates missing BETA and SE values, given a full set of Z-scores
#'
#' @param gwas of summary statistics, with partially populated BETA and SE columns, and fully popualted Z and EAF
#' @param sample_size of GWAS
#'
#' @return gwas with fully populated BETA and SE columns
populate_beta_with_known_z_scores <- function(gwas, sample_size) {
  gwas$SE_new <- 1 / sqrt(2 * gwas$EAF * (1 - gwas$EAF) * sample_size)
  gwas$BETA_new <- gwas$Z * gwas$SE_new

  correction_gwas <- dplyr::filter(gwas, !is.na(BETA) & !is.null(BETA))
  correction <- lm(correction_gwas$BETA_new ~ correction_gwas$BETA)$coef[2]
  gwas$BETA_new <- gwas$BETA_new / correction
  gwas$SE_new <- gwas$SE_new / correction

  gwas <- dplyr::mutate(gwas, BETA = dplyr::if_else(is.na(BETA, BETA_new, BETA)),
                              SE = dplyr::if_else(is.na(SE, SE_new, SE))) |>
    dplyr::select(-BETA_new, -SE_new)

  return(gwas)
}

main(args)
