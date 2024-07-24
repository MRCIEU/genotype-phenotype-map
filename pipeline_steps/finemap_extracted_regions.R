source("constants.R")

parser <- argparser::arg_parser('Finemap studies per region')
parser <- argparser::add_argument(parser, '--ld_block', help = 'LD block that the ', type = 'character')
parser <- argparser::add_argument(parser, '--completed_output_file', help = 'Completed output file', type = 'character')
args <- argparser::parse_args(parser)

main <- function(args) {
  ld_info <- ld_block_dirs(args$ld_block)
  ld_region <- vroom::vroom(paste0(ld_info$ld_matrix_prefix, '.ld'), col_names=F, show_col_types = F)
  ld_region_from_reference_panel <- vroom::vroom(paste0(ld_info$ld_matrix_prefix, '.tsv'), show_col_types = F)

  imputed_studies_file <- paste0(ld_info$ld_block_data, '/imputed_studies.tsv')
  imputed_studies <- vroom::vroom(imputed_studies_file, show_col_types = F)

  if (nrow(imputed_studies) == 0) {
    finemapped_results <- data.frame()
  } else {
    finemapped_results <- apply(imputed_studies, 1, function (study) {
      finemap_file_prefix <- sub('imputed', 'finemapped', study[['file']])
      finemap_file_prefix <- sub('\\..*', '', finemap_file_prefix)
      if (file.exists(paste0(finemap_file_prefix, "_1.tsv.gz"))) {
        return()
      }
      cat(paste0('Finemapping ', study['file'], ": "))

      dir.create(dirname(finemap_file_prefix), recursive=T, showWarnings=F)
      gwas <- vroom::vroom(study[['file']], show_col_types = F)

      if (typeof(gwas$EAF) == 'character') {
        message('EAF is not populated, cant split results, skipping.')
        failed_finemap_info <- process_unfinemapped_gwas(gwas, study, finemap_file_prefix, 'no_eaf')
        return(failed_finemap_info)
      }

      keep <- ld_region_from_reference_panel$RSID %in% gwas$RSID
      ld_for_gwas <- ld_region[keep, keep]
      ld_matrix <- matrix(as.vector(data.matrix(ld_for_gwas)), nrow=nrow(ld_for_gwas), ncol=ncol(ld_for_gwas))
      testthat::expect_true(nrow(gwas) == nrow(ld_for_gwas), 'gwas and ld matrix should match size')

      sample_size <- as.numeric(study['sample_size'])
      susie_result <- susieR::susie_rss(z=gwas$Z, R=ld_matrix, n=sample_size)

      if (susie_result$converged == F || is.null(susie_result$sets$cs) || length(susie_result$sets$cs) <= 1) {
        message('susie either didnt converge or has no more than 1 credible sets, no need to adjust.')
        failed_finemap_info <- process_unfinemapped_gwas(gwas, study, finemap_file_prefix)
        if (susie_result$converged == T) failed_finemap_info$message <- 'less_than_2_cs'
        return(failed_finemap_info)
      }

      message(paste('found', length(susie_result$sets$cs_index), 'credible sets!'))

      new_bps <- c()
      new_files <- c()
      min_ps <- c()
      unique_ids <- c()
      for (i in susie_result$sets$cs_index) {
        finemap_num <- which(i == susie_result$sets$cs_index)
        conditioned_gwas <- update_gwas_with_log_bayes_factor(gwas, susie_result$lbf_variable[i, ], sample_size)
        finemap_file <- paste0(finemap_file_prefix, '_', finemap_num, '.tsv.gz')
        unique_id <- paste0(study['study'], "_", study['chr'], "_", study['bp'], "_", i)

        vroom::vroom_write(conditioned_gwas, finemap_file)

        #this finds the lead SNP in new credible set
        important_row <- susie_result$sets$cs[paste0('L', i)][[1]][[1]]
        new_bps <- c(new_bps, as.numeric(gwas[important_row, ]$BP))
        new_files <- c(new_files, finemap_file)
        min_ps <- c(min_ps, min(gwas$P))
        unique_ids <- c(unique_ids, unique_id)
      }

      succeeded_finemap_info <- data.frame(study=study[['study']],
                                           unique_study_id=unique_ids,
                                           file=new_files,
                                           chr=as.character(study[['chr']]),
                                           bp=new_bps,
                                           p_value_threshold=as.numeric(study['p_value_threshold']),
                                           min_p=min_ps,
                                           category=study['category'],
                                           sample_size=sample_size,
                                           cis_trans=study['cis_trans'],
                                           message='success'
      )
      return(succeeded_finemap_info)
    }) |> dplyr::bind_rows()
  }

  finemapped_results_file <- paste0(ld_info$ld_block_data, '/finemapped_studies.tsv')
  if (file.exists(finemapped_results_file)) {
    existing_finemapped_results <- vroom::vroom(finemapped_results_file,
      show_col_types = F,
      col_types = vroom::cols(
        chr = vroom::col_character(),
        bp = vroom::col_number(),
        min_p = vroom::col_number(),
        p_value_threshold = vroom::col_number()
      )
    )
    finemapped_results <- dplyr::bind_rows(existing_finemapped_results, finemapped_results) |>
      dplyr::distinct()
  }

  if (nrow(finemapped_results) == 0) finemapped_results <- empty_finemapped_info()

  vroom::vroom_write(finemapped_results, finemapped_results_file)
  vroom::vroom_write(data.frame(), args$completed_output_file)
}

empty_finemapped_info <- function() {
  return(data.frame(study=character(),
                    unique_study_id=character(),
                    file=character(),
                    chr=character(),
                    bp=numeric(),
                    p_value_threshold=numeric(),
                    min_p=numeric(),
                    category=character(),
                    sample_size=numeric(),
                    cis_trans=character(),
                    message=character()
    )
  )
}

process_unfinemapped_gwas <- function(gwas, study, finemap_file_prefix, message='failed') {
  sample_size <- as.numeric(study['sample_size'])
  gwas <- populate_beta_with_known_z_scores(gwas, sample_size)
  min_p <- min(gwas$P)

  failed_finemap_file <- paste0(finemap_file_prefix, '_1.tsv.gz')
  unique_id <- paste0(study['study'], "_", study['chr'], "_", study['bp'], "_1")
  failed_finemap_info <- data.frame(study=study[['study']],
                                    unique_study_id=unique_id,
                                    file=failed_finemap_file,
                                    chr=as.character(study[['chr']]),
                                    bp=as.numeric(study['bp']),
                                    p_value_threshold=as.numeric(study['p_value_threshold']),
                                    min_p=min_p,
                                    category=study['category'],
                                    sample_size=sample_size,
                                    cis_trans=study['cis_trans'],
                                    message=message
  )

  vroom::vroom_write(gwas, failed_finemap_file)
  return(failed_finemap_info)
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
  #TODO: option 1: take allele frequency from the reference panel
  #TODO: option 2: update hyprcoloc to be able to input lbf, so we don't have to do this conversion at all
  se <- sqrt(1 / (2 * sample_size * gwas$EAF * (1-gwas$EAF)))
  r <- prior_v / (prior_v + se^2)
  z <- sqrt((2 * lbf - log(sqrt(1-r)))/r)
  beta <- z * se
  p <- abs(2 * pnorm(abs(z), lower.tail = F))

  gwas <- dplyr::mutate(gwas, BETA = beta, SE = se, P = p, Z = z) |>
    dplyr::filter(!is.na(BETA) & !is.na(SE) & BETA != Inf & SE != Inf)

  return(gwas)
}

#' Calculates missing BETA, SE, and P values, given a full set of Z-scores
#'
#' @param gwas of summary statistics, with partially populated BETA and SE columns, and fully popualted Z and EAF
#' @param sample_size of GWAS
#'
#' @return gwas with fully populated BETA and SE columns
populate_beta_with_known_z_scores <- function(gwas, sample_size) {
  gwas$SE_new <- 1 / sqrt(2 * gwas$EAF * (1 - gwas$EAF) * sample_size)
  gwas$BETA_new <- gwas$Z * gwas$SE_new

  correction_gwas <- dplyr::filter(gwas, !is.na(BETA) & !is.null(BETA))
  correction_gwas$BETA_new[which(!is.finite(correction_gwas$BETA_new))] <- NA
  correction <- lm(correction_gwas$BETA_new ~ correction_gwas$BETA, na.action=na.omit)$coef[2]

  gwas$BETA_new <- gwas$BETA_new / correction
  gwas$SE_new <- gwas$SE_new / correction
  gwas$P_new <- abs(2 * pnorm(abs(gwas$Z), lower.tail = F))

  gwas <- dplyr::mutate(gwas, BETA = dplyr::if_else(is.na(BETA), BETA_new, BETA),
                              SE = dplyr::if_else(is.na(SE), SE_new, SE),
                              P = dplyr::if_else(is.na(P), P_new, P) ) |>
    dplyr::select(-BETA_new, -SE_new, -P_new) |>
    dplyr::filter(!is.na(BETA) & !is.na(SE) & BETA != Inf & SE != Inf)

  return(gwas)
}

main(args)
