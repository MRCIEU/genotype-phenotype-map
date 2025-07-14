source("constants.R")
source("create_svgs_from_gwas.R")

parser <- argparser::arg_parser('Finemap studies per region')
parser <- argparser::add_argument(parser, '--ld_block', help = 'LD block that the ', type = 'character')
parser <- argparser::add_argument(parser, '--completed_output_file', help = 'Completed output file', type = 'character')
parser <- argparser::add_argument(parser, '--worker_guid', help = 'Worker GUID', type = 'character', default = NA)
args <- argparser::parse_args(parser)

discard_gwas_size <- 150
minimum_gwas_size <- 700
number_finemapped_results_threshold <- 3

snp_annotations <- vroom::vroom(file.path(variant_annotation_dir, "vep_annotations_hg38.tsv.gz"), show_col_types =  F) |>
  dplyr::rename_with(tolower) |>
  dplyr::select(chr, bp, snp)

main <- function() {
  if (!is.na(args$worker_guid)) {
    update_directories_for_worker(args$worker_guid)
  }
  ld_info <- ld_block_dirs(args$ld_block)

  imputed_studies_file <- glue::glue('{ld_info$ld_block_data}/imputed_studies.tsv')
  if (!file.exists(imputed_studies_file)) {
    vroom::vroom_write(data.frame(), args$completed_output_file)
    return()
  }
  imputed_studies <- vroom::vroom(imputed_studies_file, show_col_types = F) |>
    dplyr::filter(variant_type == variant_types$common)

  #TODO: gzip all unphased.vcor1 files on ieup1, then remove this
  if (!is.na(args$worker_guid)) {
    ld_matrix_file <- glue::glue('{ld_info$ld_reference_panel_prefix}.unphased.vcor1.gz')
  } else {
    ld_matrix_file <- glue::glue('{ld_info$ld_reference_panel_prefix}.unphased.vcor1')
  }
  ld_matrix <- vroom::vroom(ld_matrix_file, col_names=F, show_col_types = F)
  ld_matrix_info <- vroom::vroom(glue::glue('{ld_info$ld_reference_panel_prefix}.tsv'), show_col_types = F)


  finemapped_results_file <- glue::glue('{ld_info$ld_block_data}/finemapped_studies.tsv')
  existing_finemapped_results <- load_existing_finemapped_results(finemapped_results_file)

  if (nrow(imputed_studies) == 0) {
    finemapped_results <- data.frame()
  } else {
    finemapped_results_list <- lapply(seq_len(nrow(imputed_studies)), function(i) {
      gc()
      tryCatch({
        study <- imputed_studies[i, , drop = FALSE]
        start_time <- Sys.time()
        sample_size <- as.numeric(study['sample_size'])
        flattened_block_name <- flattened_ld_block_name(args$ld_block)
        finemap_file_prefix <- glue::glue('{extracted_study_dir}/{study$study}/finemapped/{flattened_block_name}')

        study_already_finemapped <- any(study[['study']] == existing_finemapped_results$study, na.rm = TRUE)
        if (!is.na(study_already_finemapped) && study_already_finemapped) {
          return(NULL)
        }

        #we don't want to use extracted regions with too few SNPs, due to poor finemapping results
        #and spurious coloc results due to needing to harmonise SNPs 
        gwas <- vroom::vroom(study[['file']], show_col_types = F)
        if (nrow(gwas) < discard_gwas_size) {
          return(NULL)
        }

        # GodMC methylation studies are sparsley populated, so finemapping is not useful
        if (grepl('godmc-methylation', study['study'])) {
          unfinemapped_results <- process_unfinemapped_gwas(gwas, study, finemap_file_prefix, start_time, message='sparse_population')
          unfinemapped_results <- dplyr::bind_cols(unfinemapped_results, data.frame(
            first_finemap_num_results = NA,
            second_finemap_num_results = NA,
            qc_step_run = F,
            snps_removed_by_qc = NA
          ))
          return(unfinemapped_results)
        }

        results <- run_susie_finemapping(gwas, study, ld_matrix_info, ld_matrix, finemap_file_prefix, sample_size, start_time)
        if (!is.null(results$failed_finemap_info)) {
          results$failed_finemap_info <- dplyr::bind_cols(results$failed_finemap_info, data.frame(
            first_finemap_num_results = 0,
            second_finemap_num_results = NA,
            qc_step_run = F,
            snps_removed_by_qc = NA
          ))
          return(results$failed_finemap_info)
        }
        study['first_finemap_num_results'] <- length(results$susie_result$sets$cs_index)

        #if there are a lot of susie results, run DENTIST, to see if there are any bad SNPs, then rerun susie if any SNPs are removed
        if (length(results$susie_result$sets$cs_index) > number_finemapped_results_threshold) {
          message('performing qc')
          qc_results <- perform_qc(gwas, study, ld_info$ld_reference_panel_prefix)
          study <- qc_results$study

          if (study['snps_removed_by_qc'] > 0) {
            results <- run_susie_finemapping(qc_results$gwas, study, ld_matrix_info, ld_matrix, finemap_file_prefix, sample_size, start_time)
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
            gwas <- qc_results$gwas
          }
          else {
            study['second_finemap_num_results'] <- NA
          }
        } else {
          study['qc_step_run'] <- F
          study['snps_removed_by_qc'] <- NA
          study['second_finemap_num_results'] <- NA
        }

        susie_result_to_save <- list(
          converged = results$susie_result$converged,
          cs_index = results$susie_result$sets$cs_index,
          cs = results$susie_result$sets$cs
        )
        saveRDS(susie_result_to_save, glue::glue('{finemap_file_prefix}_results.rds'))

        succeeded_finemap_info <- split_susie_result_into_conditional_gwases(results$susie_result, gwas, study, sample_size, finemap_file_prefix, start_time)
          return(succeeded_finemap_info)
      }, error = function(e) {
        message(paste('Error finemapping:', study$study, e))
        stop(e)
      })
    })

    saveRDS(finemapped_results_list, glue::glue('{finemapped_results_file}_list.rds'))
    
    # Filter out NULL results and bind rows
    finemapped_results_list <- finemapped_results_list[!sapply(finemapped_results_list, is.null)]
    if (length(finemapped_results_list) > 0) {
      finemapped_results <- dplyr::bind_rows(finemapped_results_list)
    } else {
      finemapped_results <- data.frame()
    }
  }

  finemapped_results <- dplyr::bind_rows(existing_finemapped_results, finemapped_results) |> 
    dplyr::distinct(unique_study_id, .keep_all = TRUE)

  if (nrow(finemapped_results) == 0) finemapped_results <- empty_finemapped_info()

  vroom::vroom_write(finemapped_results, finemapped_results_file)
  vroom::vroom_write(data.frame(), args$completed_output_file)
}

load_existing_finemapped_results <- function(finemapped_results_file) {
  if (file.exists(finemapped_results_file)) {
    return(vroom::vroom(finemapped_results_file,
                        show_col_types = F,
                        col_types = finemapped_column_types 
    ))
  } else {
    return(empty_finemapped_info())
  }
}

empty_finemapped_info <- function() {
  return(data.frame(study=character(),
                    unique_study_id=character(),
                    ld_block=character(),
                    variant_type=character(),
                    file=character(),
                    ancestry=character(),
                    chr=character(),
                    bp=numeric(),
                    snp=character(),
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
                    time_taken=character(),
                    svg_file=character(),
                    file_with_lbfs=character(),
                    ignore=logical()
    )
  )
}

run_susie_finemapping <- function(gwas, study, ld_matrix_info, ld_matrix, finemap_file_prefix, sample_size, start_time) {
  failed_finemap_info <- NULL
  susie_result <- list(converged=F)

  if (nrow(gwas) < minimum_gwas_size) {
    failed_finemap_info <- process_unfinemapped_gwas(gwas, study, finemap_file_prefix, start_time)
    failed_finemap_info$finemap_message <- 'too_small_to_finemap'
    return(list(susie_result=susie_result, failed_finemap_info=failed_finemap_info))
  }

  keep <- ld_matrix_info$SNP %in% gwas$SNP
  ld_matrix_subset <- ld_matrix[keep, keep]
  ld_matrix_subset <- matrix(as.vector(data.matrix(ld_matrix_subset)), nrow=nrow(ld_matrix_subset), ncol=ncol(ld_matrix_subset))
  if (nrow(gwas) != nrow(ld_matrix_subset)) {
    stop(paste('Error:', study[['file']], 'GWAS', nrow(gwas), 'and ld matrix', nrow(ld_matrix_subset), 'should match size'))
  }

  tryCatch(expr = {
    susie_result <- susieR::susie_rss(z=gwas$Z, R=ld_matrix_subset, n=sample_size)
  }, error = function(e) {
    message(paste('Susie error for', study[['file']], ':', e$message))
    susie_result <<- list(converged=F, sets=list(cs_index=c(), cs=list()))
  })

  if (susie_result$converged == F || is.null(susie_result$sets$cs) || length(susie_result$sets$cs) <= 1) {
    new_bp <- NA
    new_snp <- NA
    if (susie_result$converged == T && length(susie_result$sets$cs_index) == 1) {
      #this finds the lead SNP in new credible set
      important_row <- susie_result$sets$cs[1][[1]][[1]]
      new_bp <- as.numeric(gwas[important_row, ]$BP)
      new_snp <- gwas[important_row, ]$SNP
    } else {
      new_bp <- gwas[which.min(gwas$P), ]$BP
    }

    failed_finemap_info <- process_unfinemapped_gwas(gwas, study, finemap_file_prefix, start_time, new_bp=new_bp, new_snp=new_snp)
    if (susie_result$converged == F) {
      message(paste('Finemapping:', study['file'], 'susie didnt converge'))
    }
    else if (susie_result$converged == T) {
      message(paste('Finemapping:', study['file'], 'susie found 1 credible set'))
      failed_finemap_info$finemap_message <- 'less_than_2_cs'
    }
  }
  else {
    message(paste('Finemapping:', study['file'], 'found', length(susie_result$sets$cs_index), 'credible sets in', susie_result$niter, 'iterations!'))
  }
  return(list(susie_result=susie_result, failed_finemap_info=failed_finemap_info))
}

process_unfinemapped_gwas <- function(gwas, study, finemap_file_prefix, start_time, message='failed', new_bp=NA, new_snp=NA) {
  sample_size <- as.numeric(study['sample_size'])
  min_p <- min(gwas$P)

  if (!is.na(new_bp)) {
    study['bp'] <- new_bp
  }
  if (!is.na(new_snp)) {
    study['snp'] <- new_snp
  } else {
    snp_entry <- snp_annotations |>
      dplyr::filter(chr == as.numeric(study[['chr']]) & bp == as.numeric(study['bp']))

    if (nrow(snp_entry) == 0) {
      message('finding new snp for: ', study[['chr']], ':', study['bp'])
      study['snp'] <- gwas[which.min(gwas$P), ]$SNP
    } else {
      study['snp'] <- snp_entry |>
        dplyr::slice_head(n = 1) |>
        dplyr::pull(snp)
    }
  }
  gwas <- dplyr::mutate(gwas, LBF = convert_z_to_lbf(Z, SE, EAF, study$sample_size, study$category))

  failed_finemap_file <- glue::glue('{finemap_file_prefix}_1.tsv.gz')
  unique_id <- glue::glue('{study["study"]}_{args$ld_block}_1')

  finemap_gwas <- dplyr::select(gwas, SNP, CHR, BP, SE, EAF, Z, IMPUTED, LBF)
  vroom::vroom_write(finemap_gwas, failed_finemap_file)

  file_with_lbfs <- glue::glue('{finemap_file_prefix}_with_lbf.tsv.gz')
  write_lbf_gwas(gwas, NULL, file_with_lbfs)

  block_name <- basename(failed_finemap_file) |> stringr::str_replace("\\.tsv\\.gz$", "")
  svg_file <- glue::glue("{extracted_study_dir}{study$study}/svgs/extractions/{block_name}.svg")
  create_svg_for_ld_block(finemap_gwas, study$study, svg_file)

  failed_finemap_info <- data.frame(study=study[['study']],
                                    unique_study_id=unique_id,
                                    ld_block = args$ld_block,
                                    variant_type=study[['variant_type']],
                                    file=failed_finemap_file,
                                    ancestry=study[['ancestry']],
                                    chr=as.character(study[['chr']]),
                                    bp=as.numeric(study['bp']),
                                    snp=study['snp'],
                                    p_value_threshold=as.numeric(study['p_value_threshold']),
                                    min_p=min_p,
                                    category=study['category'],
                                    sample_size=sample_size,
                                    cis_trans=study['cis_trans'],
                                    finemap_message=message,
                                    time_taken=as.character(hms::as_hms(difftime(Sys.time(), start_time))),
                                    svg_file=svg_file,
                                    file_with_lbfs=file_with_lbfs,
                                    ignore=F
  )

  return(failed_finemap_info)
}

split_susie_result_into_conditional_gwases <- function(susie_result, gwas, study, sample_size, finemap_file_prefix, start_time) {
  new_snps <- c()
  new_bps <- c()
  new_files <- c()
  min_ps <- c()
  unique_ids <- c()
  svg_files <- c()
  lbf_columns <- list()

  original_min_p <- min(gwas$P, na.rm = T)

  for (i in susie_result$sets$cs_index) {
    finemap_num <- which(i == susie_result$sets$cs_index)
    #occasionally, the new min p is smaller than the original min p, which is wrong, so we set it to the original min p
    p_value_gwas <- update_gwas_with_log_bayes_factor(gwas, susie_result$lbf_variable[i, ], sample_size)
    new_min_p <- min(p_value_gwas$P, na.rm = T)

    if (new_min_p < original_min_p) {
      new_min_p <- original_min_p
    }

    min_ps <- c(min_ps, new_min_p)

    conditioned_gwas <- dplyr::select(gwas, SNP, CHR, BP, SE, EAF, P, IMPUTED) |>
      dplyr::mutate(LBF = susie_result$lbf_variable[i, ])

    lbf_columns[[glue::glue('LBF_{finemap_num}')]] <- susie_result$lbf_variable[i, ]
    finemap_file <- glue::glue('{finemap_file_prefix}_{finemap_num}.tsv.gz')

    #this finds the lead SNP in new credible set
    important_row <- susie_result$sets$cs[paste0('L', i)][[1]][[1]]
    new_snps <- c(new_snps, gwas[important_row, ]$SNP)
    new_bp <- as.numeric(gwas[important_row, ]$BP)
    new_bps <- c(new_bps, new_bp)
    new_files <- c(new_files, finemap_file)


    unique_id <- glue::glue('{study["study"]}_{args$ld_block}_{finemap_num}')
    unique_ids <- c(unique_ids, unique_id)

    block_name <- basename(finemap_file) |> stringr::str_replace("\\.tsv\\.gz$", "")
    svg_file <- glue::glue("{extracted_study_dir}{study$study}/svgs/extractions/{block_name}.svg")
    svg_files <- c(svg_files, svg_file)
    create_svg_for_ld_block(conditioned_gwas, study$study, svg_file)

    # if the new credible set's bp is less than 2MB from the original bp, mark as cis, otherwise trans
    if (!is.na(study['cis_trans']) && study['cis_trans'] == cis_trans$cis_only) {
      if (abs(as.numeric(study['bp']) - new_bp) < 1000000) {
        study['cis_trans'] <- cis_trans$cis_only
      } else {
        study['cis_trans'] <- cis_trans$trans_only
      }
    }

    conditioned_gwas <- dplyr::select(conditioned_gwas, -P)
    vroom::vroom_write(conditioned_gwas, finemap_file)
  }

  lbf_columns <- as.data.frame(lbf_columns)
  file_with_lbfs <- glue::glue('{finemap_file_prefix}_with_lbf.tsv.gz')
  write_lbf_gwas(gwas, lbf_columns, file_with_lbfs)
  time_taken <- as.character(hms::as_hms(difftime(Sys.time(), start_time)))

  num_credible_sets <- length(unique_ids)

  succeeded_finemap_info <- data.frame(
    study = rep(study[['study']], num_credible_sets),
    unique_study_id = unique_ids,
    ld_block = rep(args$ld_block, num_credible_sets),
    variant_type = rep(study[['variant_type']], num_credible_sets),
    file = new_files,
    ancestry = rep(study[['ancestry']], num_credible_sets),
    snp = new_snps,
    chr = rep(as.character(study[['chr']]), num_credible_sets),
    bp = new_bps,
    p_value_threshold = rep(as.numeric(study[['p_value_threshold']]), num_credible_sets),
    min_p = min_ps,
    category = rep(study[['category']], num_credible_sets),
    sample_size = rep(sample_size, num_credible_sets),
    cis_trans = rep(study[['cis_trans']], num_credible_sets),
    finemap_message = rep('success', num_credible_sets),
    first_finemap_num_results = rep(as.numeric(study[['first_finemap_num_results']]), num_credible_sets),
    second_finemap_num_results = rep(as.numeric(study[['second_finemap_num_results']]), num_credible_sets),
    qc_step_run = rep(as.logical(study[['qc_step_run']]), num_credible_sets),
    snps_removed_by_qc = rep(as.numeric(study[['snps_removed_by_qc']]), num_credible_sets),
    time_taken = rep(time_taken, num_credible_sets),
    svg_file = svg_files,
    file_with_lbfs = rep(file_with_lbfs, num_credible_sets),
    ignore = rep(FALSE, num_credible_sets)
  )
  return(succeeded_finemap_info)
}

write_lbf_gwas <- function(gwas, lbf_columns, lbf_file) {
  lbf_gwas <- dplyr::select(gwas,
    dplyr::any_of(c('SNP', 'CHR', 'BP', 'EA', 'OA', 'EAF', 'Z', 'BETA', 'SE', 'IMPUTED', 'LBF'))
  )

  if (!is.null(lbf_columns) && nrow(lbf_columns) > 0) {
    lbf_gwas <- dplyr::bind_cols(lbf_gwas, lbf_columns)
  } else {
    lbf_gwas <- dplyr::rename(lbf_gwas, LBF_1 = LBF)
  }
  vroom::vroom_write(lbf_gwas, lbf_file)
}

#' perform_qc: 
#' 
perform_qc <- function(gwas, study, bfile) {
  study['qc_step_run'] <- T
  gwas <- populate_beta_with_known_z_scores(gwas, as.numeric(study['sample_size']))
  dentist_gwas <- dplyr::rename(gwas, A1='EA', A2='OA', freq='EAF', beta='BETA', se='SE', p='P') |>
    dplyr::mutate(N = as.numeric(study[['sample_size']])) |>
    dplyr::select(SNP, A1, A2, freq, beta, se, p, N)

  dentist_tmp_file <- withr::local_tempfile()
  vroom::vroom_write(dentist_gwas, dentist_tmp_file)

  dentist_command <- paste('DENTIST --bfile', bfile,
                           '--thread-num 8 ',
                           '--gwas-summary', dentist_tmp_file,
                           '--chrID', study[['chr']],
                           '--out', dentist_tmp_file
  )
  system(dentist_command, wait = T, ignore.stdout = T)

  dentist_file_to_remove <- glue::glue('{dentist_tmp_file}.DENTIST.short.txt')
  dentist_full_file <- glue::glue('{dentist_tmp_file}.DENTIST.full.txt')

  if (!file.exists(dentist_file_to_remove)) {
    message('DENTIST command failed')
    study['snps_removed_by_qc'] <- 0
  }
  else {
    dentist_to_remove <- vroom::vroom(dentist_file_to_remove, col_names = F, delim = ' ', show_col_types = F)
    if (nrow(dentist_to_remove) > 0) {
      message('removing ', nrow(dentist_to_remove), ' snps from gwas')
      gwas <- dplyr::filter(gwas, !SNP %in% dentist_to_remove$X1) |>
        dplyr::select(SNP, RSID, dplyr::everything())
      vroom::vroom_write(gwas, study[['file']])

      dentist_file <- sub('.tsv.gz', '_dentist_removed.tsv', study[['file']])
      dentist_full_remove <- vroom::vroom(dentist_full_file, col_names = F, show_col_types = F) |>
        dplyr::filter(X1 %in% dentist_to_remove$X1) |>
        dplyr::rename(SNP='X1', chisq='X2', nlogp='X3', dup='X4')

      vroom::vroom_write(dentist_full_remove, dentist_file)
    }

    study['snps_removed_by_qc'] <- nrow(dentist_to_remove)
  }
  message(paste('keeping:', nrow(gwas), 'deleting:', study[['snps_removed_by_qc']]))
  return(list(gwas=gwas, study=study))
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

main()