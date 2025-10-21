source('constants.R')
source('imputation_method.R')

imputation_correlation_threshold <- 0.7
imputation_significant_rows_overinflated_threshold <- 100

parser <- argparser::arg_parser('Impute GWASes for pipeline')
parser <- argparser::add_argument(parser, '--ld_block', help = 'LD block that the ', type = 'character')
parser <- argparser::add_argument(parser, '--completed_output_file', help = 'Completed output file', type = 'character')
parser <- argparser::add_argument(parser, '--worker_guid', help = 'Worker GUID', type = 'character', default = NA)
args <- argparser::parse_args(parser)

main <- function() {
  if (!is.na(args$worker_guid)) {
    update_directories_for_worker(args$worker_guid)
  }
  ld_info <- ld_block_dirs(args$ld_block)
  ld_matrix_info <- vroom::vroom(glue::glue('{ld_info$ld_reference_panel_prefix}.tsv'), show_col_types = F)
  #TODO: gzip all unphased.vcor1 files on ieup1, then remove this
  if (!is.na(args$worker_guid)) {
    ld_matrix_file <- glue::glue('{ld_info$ld_reference_panel_prefix}.unphased.vcor1.gz')
  } else {
    ld_matrix_file <- glue::glue('{ld_info$ld_reference_panel_prefix}.unphased.vcor1')
  }
  ld_matrix <- vroom::vroom(ld_matrix_file, col_names=F, show_col_types = F)
  ld_matrix_eig <- readRDS(glue::glue('{ld_info$ld_reference_panel_prefix}.ldeig.rds'))

  standardised_studies_file <- glue::glue('{ld_info$ld_block_data}/standardised_studies.tsv')
  standardised_studies  <- vroom::vroom(standardised_studies_file , show_col_types = F) |>
    dplyr::filter(variant_type == variant_types$common)

  imputed_studies_file <- glue::glue('{ld_info$ld_block_data}/imputed_studies.tsv')
  if (file.exists(imputed_studies_file)) {
    existing_imputed_studies <- vroom::vroom(imputed_studies_file, show_col_types = F, col_types = imputed_column_types)
  } else {
    existing_imputed_studies <- empty_imputed_studies()
  }

  if (nrow(standardised_studies) > 0) {
    imputed_studies <- apply(standardised_studies, 1, function (study) {
      start_time <- Sys.time()
      imputed_file <- sub('standardised', 'imputed', study[['file']])

      if (imputed_file %in% existing_imputed_studies$file) {
        return()
      }
      message('Imputing ', imputed_file)

      gwas <- vroom::vroom(study['file'], show_col_types = F)

      gwas_to_impute <- dplyr::left_join(
        dplyr::select(ld_matrix_info, -EAF),
        dplyr::select(gwas, -CHR, -BP, -EA, -OA),
        by=dplyr::join_by(SNP)
      ) 

      rows_to_impute <- !ld_matrix_info$SNP %in% gwas$SNP
      gwas_to_impute$EAF[rows_to_impute] <- ld_matrix_info$EAF[rows_to_impute]

      # Imputation is not useful or effective for sparsely populated studies, instead pad missing values
      if (study[['coverage']] == coverage_types$sparse) {
        filtered_results <- list(
          significant_rows_imputed = NA,
          significant_rows_filtered = NA
        )
        result <- pad_missing_values(gwas_to_impute)
        vroom::vroom_write(result$gwas, imputed_file)
      } else {
        result <- perform_imputation(imputed_file, gwas_to_impute, ld_matrix_eig)
        verify_imputation_results(result$gwas, imputed_file)
        filtered_results <- filter_imputation_results(result$gwas, ld_matrix, min(gwas$BP), max(gwas$BP))

        if((!is.na(result$b_cor) && result$b_cor >= imputation_correlation_threshold) || filtered_results$significant_rows_filtered < imputation_significant_rows_overinflated_threshold) {
          vroom::vroom_write(filtered_results$gwas, imputed_file)
        } else {
          vroom::vroom_write(gwas, imputed_file)
        }
      }

      time_taken <- as.character(hms::as_hms(difftime(Sys.time(), start_time)))

      imputation_info <- data.frame(
        study=study[['study']],
        file=imputed_file,
        ancestry=study['ancestry'],
        chr = as.character(study[['chr']]),
        bp = as.numeric(study[['bp']]),
        p_value_threshold = as.numeric(study[['p_value_threshold']]),
        category=study['category'],
        sample_size=as.numeric(study['sample_size']),
        cis_trans=study['cis_trans'],
        rows_imputed=result$rows_imputed,
        b_cor=result$b_cor,
        se_cor=result$se_cor,
        z_adj=result$z_adj_coef1,
        se_adj=result$se_adj_coef1,
        time_taken=time_taken,
        significant_rows_imputed=filtered_results$significant_rows_imputed,
        significant_rows_filtered=filtered_results$significant_rows_filtered,
        ld_block=ld_info$block,
        variant_type=study['variant_type'],
        coverage=study['coverage']
      )

      return(imputation_info)
    }) |> dplyr::bind_rows()

    if (nrow(imputed_studies) > 0) {
      imputed_studies <- dplyr::bind_rows(existing_imputed_studies, imputed_studies) |>
        dplyr::distinct(study, .keep_all = TRUE)
      vroom::vroom_write(imputed_studies, imputed_studies_file)
    }
  }

  vroom::vroom_write(data.frame(), args$completed_output_file)
}

empty_imputed_studies <- function() {
  return(
    data.frame(
      study = character(),
      file = character(),
      ancestry = character(),
      chr = character(),
      bp = numeric(),
      p_value_threshold = numeric(),
      category = character(),
      sample_size = numeric(),
      cis_trans = character(),
      rows_imputed = numeric(),
      b_cor = numeric(),
      se_cor = numeric(),
      z_adj = numeric(),
      se_adj = numeric(),
      time_taken = character(),
      ld_block = character(),
      variant_type = character(),
      coverage = character()
    )
  )
}

pad_missing_values <- function(gwas) {
  gwas$BETA[is.na(gwas$BETA)] <- 0
  gwas$SE[is.na(gwas$SE)] <- 1
  gwas$P[is.na(gwas$P)] <- 1
  gwas$Z[is.na(gwas$Z)] <- 0
  gwas$IMPUTED <- ifelse(gwas$BETA == 0 & gwas$SE == 1, T, F)

  return(list(
    gwas = gwas,
    rows_imputed = 0,
    b_cor = NA,
    se_cor = NA,
    z_adj_coef1 = NA,
    se_adj_coef1 = NA
  ))
}

verify_imputation_results <- function(gwas, imputed_file) {
  if (any(is.na(gwas$BETA))) {
    stop(glue::glue('BETA: NA in imputed study {imputed_file}'))
  }
  if (any(is.na(gwas$SE) | gwas$SE <= 0)) {
    stop(glue::glue('SE: NA or <= 0 in imputed study {imputed_file}'))
  }
  if (any(is.na(gwas$P) | gwas$P < 0 | gwas$P > 1)) {
    stop(glue::glue('P: NA or < 0 or > 1 in imputed study {imputed_file}'))
  }
}
main()
