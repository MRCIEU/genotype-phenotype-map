source('constants.R')
source('imputation_method.R')

imputation_correlation_threshold <- 0.7
p_value_filter_correlation_threshold <- 0.6

parser <- argparser::arg_parser('Impute GWASes for pipeline')
parser <- argparser::add_argument(parser, '--ld_block', help = 'LD block that the ', type = 'character')
parser <- argparser::add_argument(parser, '--completed_output_file', help = 'Completed output file', type = 'character')
args <- argparser::parse_args(parser)

main <- function() {
  ld_info <- ld_block_dirs(args$ld_block)
  ld_matrix_info <- vroom::vroom(glue::glue('{ld_info$ld_reference_panel_prefix}.tsv'), show_col_types = F)
  ld_matrix <- vroom::vroom(glue::glue('{ld_info$ld_reference_panel_prefix}.unphased.vcor1'), col_names=F, show_col_types = F)
  ld_matrix_eig <- readRDS(glue::glue('{ld_info$ld_reference_panel_prefix}.ldeig.rds'))

  standardised_studies_file <- glue::glue('{ld_info$ld_block_data}/standardised_studies.tsv')
  standardised_studies  <- vroom::vroom(standardised_studies_file , show_col_types = F)

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

      result <- perform_imputation(imputed_file, gwas_to_impute, ld_matrix_eig)
      filtered_results <- filter_imputation_results(result$gwas, ld_matrix, min(gwas$BP), max(gwas$BP))

      if(!is.na(result$b_cor) && result$b_cor >= imputation_correlation_threshold) {
        vroom::vroom_write(filtered_results$gwas, imputed_file)
      } else {
        vroom::vroom_write(gwas, imputed_file)
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
        z_adj=result$z_adj,
        se_adj=result$se_adj,
        time_taken=time_taken,
        significant_rows_imputed=filtered_results$significant_rows_imputed,
        significant_rows_filtered=filtered_results$significant_rows_filtered,
        ld_block=ld_info$block
      )

      return(imputation_info)
    }) |> dplyr::bind_rows()

    if (nrow(imputed_studies) > 0) {
      imputed_studies <- dplyr::bind_rows(existing_imputed_studies, imputed_studies) |>
        dplyr::distinct()
      vroom::vroom_write(imputed_studies, imputed_studies_file)
    }
  }

  vroom::vroom_write(data.frame(), args$completed_output_file)
}



#' filter_imputation_results: remove results with over-inflated p-vlaues
#' 
#' @param imputed gwas
#' @param ld correlation matrix
#' 
#' First, filters the imputation results to inside the range of the original gwas
#' Second, finds imputed variants with low pvalues, non-imputed rows that are correlated
#' If imputed variant is more significant the non-imputed rows, drop that variant
#' @return filtered gwas
filter_imputation_results <- function(gwas, ld_matrix, min_bp, max_bp) {
  only_keep_inside_gwas_range <- gwas$BP > min_bp & gwas$BP < max_bp 
  gwas <- gwas[only_keep_inside_gwas_range, ]

  snps_to_remove <- c()
  snps_to_investigate <- which(gwas$P < lowest_p_value_threshold & gwas$IMPUTED == T)
  for (snp_location in snps_to_investigate) {
    ld_correlations <- which(c(ld_matrix[snp_location, ]) > p_value_filter_correlation_threshold)
    ld_correlations <- ld_correlations[ld_correlations != snp_location]
    gwas_correlations <- gwas[(1:nrow(gwas) %in% ld_correlations) & gwas$IMPUTED == F, ]

    snp <- gwas[snp_location, ]

    if (length(gwas_correlations$P) > 0 && min(gwas_correlations$P * 0.1) > snp$P) {
      snps_to_remove <- c(snps_to_remove, snp_location)
    }
  }

  if (length(snps_to_remove) > 0) gwas <- gwas[-snps_to_remove, ]
  return(list(gwas = gwas,
    significant_rows_imputed = length(snps_to_investigate),
    significant_rows_filtered = length(snps_to_remove))
  )
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
      ld_block = character()
    )
  )
}

main()
