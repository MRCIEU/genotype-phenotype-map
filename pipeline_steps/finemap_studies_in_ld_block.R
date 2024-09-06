source("constants.R")

parser <- argparser::arg_parser('Finemap studies per region')
parser <- argparser::add_argument(parser, '--ld_block', help = 'LD block that the ', type = 'character')
parser <- argparser::add_argument(parser, '--complex_block', help = 'Is the ld block complex', default = F, type = 'logical')
parser <- argparser::add_argument(parser, '--completed_output_file', help = 'Completed output file', type = 'character')
args <- argparser::parse_args(parser)

main <- function(args) {
  ld_info <- ld_block_dirs(args$ld_block)
  ld_matrix <- vroom::vroom(paste0(ld_info$ld_reference_panel_prefix, '.unphased.vcor1'), col_names=F, show_col_types = F)
  ld_region_from_reference_panel <- vroom::vroom(paste0(ld_info$ld_reference_panel_prefix, '.tsv'), show_col_types = F)

  imputed_studies_file <- paste0(ld_info$ld_block_data, '/imputed_studies.tsv')
  imputed_studies <- vroom::vroom(imputed_studies_file, show_col_types = F)

  finemapped_results_file <- paste0(ld_info$ld_block_data, '/finemapped_studies.tsv')
  existing_finemapped_results <- load_existing_finemapped_results(finemapped_results_file)

  if (nrow(imputed_studies) == 0) {
    finemapped_results <- data.frame()
  } else {
    finemapped_results <- apply(imputed_studies, 1, function (study) {
      start_time <- Sys.time()
      sample_size <- as.numeric(study['sample_size'])
      finemap_file_prefix <- sub('imputed', 'finemapped', study[['file']])
      finemap_file_prefix <- sub('\\..*', '', finemap_file_prefix)
      if (any(grepl(finemap_file_prefix, existing_finemapped_results$file))) {
        return()
      }
      gwas <- vroom::vroom(study[['file']], show_col_types = F)

      results <- run_susie_finemapping(gwas, study, ld_region_from_reference_panel, ld_matrix, finemap_file_prefix, sample_size, start_time)
      if (!is.null(results$failed_finemap_info)) {
        results$failed_finemap_info <- dplyr::bind_cols(results$failed_finemap_info, data.frame(
          first_finemap_num_results = 0,
          second_finemap_num_results = NA,
          qc_step_run = F,
          snps_removed_by_qc = 0
        ))
        return(results$failed_finemap_info)
      }
      study['first_finemap_num_results'] <- length(results$susie_result$sets$cs_index)

      #if there are a lot of susie results, run DENTIST, to see if there are any bad SNPs, then rerun susie if any SNPs are removed
      if (length(results$susie_result$sets$cs_index) > 3) {
        message('performing qc')
        qc_results <- perform_qc(gwas, study, ld_info$ld_reference_panel_prefix)
        study <- qc_results$study
        if (study['snps_removed_by_qc'] > 0) {
          message('doing second finemapping...')
          results <- run_susie_finemapping(qc_results$gwas, study, ld_region_from_reference_panel, ld_matrix, finemap_file_prefix, sample_size, start_time)
          if (!is.null(results$failed_finemap_info)) {
            results$failed_finemap_info <- dplyr::bind_cols(results$failed_finemap_info, data.frame(
              first_finemap_num_results = as.numeric(study['first_finemap_num_results']),
              second_finemap_num_results = 0,
              qc_step_run = T,
              snps_removed_by_qc = as.numeric(study['snps_removed_by_qc'])
            ))
            return(results$failed_finemap_info)
          }
          study['second_finemap_num_results'] <- length(results$susie_result$sets$cs_index)
        }
      } else {
        study['qc_step_run'] <- F
        study['snps_removed_by_qc'] <- 0
      }

      succeeded_finemap_info <- split_susie_result_into_conditional_gwases(results$susie_result, gwas, study, sample_size, finemap_file_prefix, start_time)
      return(succeeded_finemap_info)
    }) |> dplyr::bind_rows()
  }

  finemapped_results <- dplyr::bind_rows(existing_finemapped_results, finemapped_results) |> dplyr::distinct()
  if (nrow(finemapped_results) == 0) finemapped_results <- empty_finemapped_info()

  vroom::vroom_write(finemapped_results, finemapped_results_file)
  vroom::vroom_write(data.frame(), args$completed_output_file)
}

load_existing_finemapped_results <- function(finemapped_results_file) {
  if (file.exists(finemapped_results_file)) {
    return(vroom::vroom(finemapped_results_file,
                        show_col_types = F,
                        col_types = vroom::cols(
                          chr = vroom::col_character(),
                          bp = vroom::col_number(),
                          min_p = vroom::col_number(),
                          sample_size = vroom::col_number(),
                          p_value_threshold = vroom::col_number(),
                          first_finemap_num_results = vroom::col_number(),
                          second_finemap_num_results = vroom::col_number(),
                          qc_step_run = vroom::col_logical(),
                          snps_removed_by_qc = vroom::col_numeric()
                        )
    ))
  } else {
    return(empty_finemapped_info())
  }
}

empty_finemapped_info <- function() {
  return(data.frame(study=character(),
                    unique_study_id=character(),
                    file=character(),
                    ancestry=character(),
                    chr=character(),
                    bp=numeric(),
                    p_value_threshold=numeric(),
                    min_p=numeric(),
                    category=character(),
                    sample_size=numeric(),
                    cis_trans=character(),
                    finemap_message=character(),
                    first_finemap_num_results=numeric(),
                    second_finemap_num_results=numeric(),
                    qc_step_run=logical(),
                    snps_removed_by_qc=numeric(),
                    time_taken=character()
    )
  )
}

run_susie_finemapping <- function(gwas, study, ld_region_from_reference_panel, ld_matrix, finemap_file_prefix, sample_size, start_time) {
  keep <- ld_region_from_reference_panel$SNP %in% gwas$SNP
  ld_matrix_subset <- ld_matrix[keep, keep]
  ld_matrix_subset <- matrix(as.vector(data.matrix(ld_matrix_subset)), nrow=nrow(ld_matrix_subset), ncol=ncol(ld_matrix_subset))
  if (nrow(gwas) != nrow(ld_matrix_subset)) {
    stop(paste('Error:', study[['file']], 'GWAS', nrow(gwas), 'and ld matrix', nrow(ld_matrix_subset), 'should match size'))
  }
  failed_finemap_info <- NULL
  susie_result <- list(converged=F)

  tryCatch(expr = {
    susie_result <- susieR::susie_rss(z=gwas$Z, R=ld_matrix_subset, n=sample_size)
  }, error = function(e) {
    message(e)
  })

  if (susie_result$converged == F || is.null(susie_result$sets$cs) || length(susie_result$sets$cs) <= 1) {
    message(paste('Finemapping:', study['file'], 'susie either didnt converge'))
    failed_finemap_info <- process_unfinemapped_gwas(gwas, study, finemap_file_prefix, start_time)
    if (susie_result$converged == T) failed_finemap_info$finemap_message <- 'less_than_2_cs'
  }
  else {
    message(paste('Finemapping:', study['file'], 'found', length(susie_result$sets$cs_index), 'credible sets in', susie_result$niter, 'iterations!'))
  }
  return(list(susie_result=susie_result, failed_finemap_info=failed_finemap_info))
}

process_unfinemapped_gwas <- function(gwas, study, finemap_file_prefix, start_time, message='failed') {
  sample_size <- as.numeric(study['sample_size'])
  gwas <- populate_beta_with_known_z_scores(gwas, sample_size)
  min_p <- min(gwas$P)

  failed_finemap_file <- paste0(finemap_file_prefix, '_1.tsv.gz')
  unique_id <- paste0(study['study'], "_", study['ancestry'], '_', study['chr'], "_", trimws(study['bp']), "_1")
  failed_finemap_info <- data.frame(study=study[['study']],
                                    unique_study_id=unique_id,
                                    file=failed_finemap_file,
                                    ancestry=study[['ancestry']],
                                    chr=as.character(study[['chr']]),
                                    bp=as.numeric(study['bp']),
                                    p_value_threshold=as.numeric(study['p_value_threshold']),
                                    min_p=min_p,
                                    category=study['category'],
                                    sample_size=sample_size,
                                    cis_trans=study['cis_trans'],
                                    finemap_message=message,
                                    time_taken=as.character(hms::as_hms(difftime(Sys.time(), start_time)))
  )

  vroom::vroom_write(gwas, failed_finemap_file)
  return(failed_finemap_info)
}

split_susie_result_into_conditional_gwases <- function(susie_result, gwas, study, sample_size, finemap_file_prefix, start_time) {
      new_bps <- c()
      new_files <- c()
      min_ps <- c()
      unique_ids <- c()
      for (i in susie_result$sets$cs_index) {
        finemap_num <- which(i == susie_result$sets$cs_index)
        conditioned_gwas <- update_gwas_with_log_bayes_factor(gwas, susie_result$lbf_variable[i, ], sample_size)
        finemap_file <- paste0(finemap_file_prefix, '_', finemap_num, '.tsv.gz')
        unique_id <- paste0(study['study'], '_', study['ancestry'], '_', study['chr'], '_', trimws(study['bp']), '_', finemap_num)
        vroom::vroom_write(conditioned_gwas, finemap_file)

        #this finds the lead SNP in new credible set
        important_row <- susie_result$sets$cs[paste0('L', i)][[1]][[1]]
        new_bps <- c(new_bps, as.numeric(gwas[important_row, ]$BP))
        new_files <- c(new_files, finemap_file)
        min_ps <- c(min_ps, min(conditioned_gwas$P, na.rm = F))
        unique_ids <- c(unique_ids, unique_id)
      }
      time_taken <- as.character(hms::as_hms(difftime(Sys.time(), start_time)))

      succeeded_finemap_info <- data.frame(study=study[['study']],
                                           unique_study_id=unique_ids,
                                           file=new_files,
                                           ancestry=study[['ancestry']],
                                           chr=as.character(study[['chr']]),
                                           bp=new_bps,
                                           p_value_threshold=as.numeric(study['p_value_threshold']),
                                           min_p=min_ps,
                                           category=study['category'],
                                           sample_size=sample_size,
                                           cis_trans=study['cis_trans'],
                                           finemap_message='success',
                                           qc_step_run=as.logical(study[['qc_step_run']]),
                                           snps_removed_by_qc=as.numeric(study[['snps_removed_by_qc']]),
                                           time_taken=time_taken
      )
      return(succeeded_finemap_info)
}

#' perform_qc: 
#' 
perform_qc <- function(gwas, study, bfile) {
  study['qc_step_run'] <- T
  gwas <- populate_beta_with_known_z_scores(gwas, as.numeric(study['sample_size']))
  dentist_gwas <- dplyr::rename(gwas, A1='EA', A2='OA', freq='EAF', beta='BETA', se='SE', p='P') |>
    dplyr::mutate(N = as.numeric(study[['sample_size']])) |>
    dplyr::select(SNP, A1, A2, freq, beta, se, p, N)

  dentist_tmp_file <- tempfile(basename(study[['file']]))
  vroom::vroom_write(dentist_gwas, dentist_tmp_file)

  dentist_command <- paste('DENTIST --bfile', bfile,
                           '--gwas-summary', dentist_tmp_file,
                           '--chrID', study[['chr']],
                           '--out', dentist_tmp_file
  )
  system(dentist_command, wait = T)

  dentist_file_to_remove <- paste0(dentist_tmp_file, '.DENTIST.short.txt')
  dentist_full_file <- paste0(dentist_tmp_file, '.DENTIST.full.txt')
  if (!file.exists(dentist_file_to_remove)) {
    message('DENTIST command failed')
    study['snps_removed_by_qc'] <- 0
  }
  else {
    dentist_to_remove <- vroom::vroom(dentist_file_to_remove, col_names = F, delim = ' ', show_col_types = F)
    if (nrow(dentist_to_remove) > 0) {
      gwas <- dplyr::filter(gwas, !SNP %in% dentist_to_remove$X1) |>
        dplyr::select(SNP, RSID, dplyr::everything())
      dentist_file <- sub('.tsv.gz', '_dentist_removed.tsv', study['file'])
      dentist_full_remove <- vroom::vroom(dentist_full_file, col_names = F, show_col_types = F) |>
        dplyr::filter(X1 %in% dentist_to_remove$X1) |>
        dplyr::rename(SNP='X1', chisq='X2', nlogp='X3', dup='X4')
      print(head(dentist_full_remove))
      vroom::vroom_write(dentist_full_remove, dentist_file)
    }

    study['snps_removed_by_qc'] <- nrow(dentist_to_remove)
  }
  message(paste('keeping:', nrow(gwas), 'deleting:', study['snps_removed_by_qc']))
  return(list(gwas=gwas, study=study))
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
