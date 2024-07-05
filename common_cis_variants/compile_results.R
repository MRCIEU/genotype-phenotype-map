source('constants.R')

POSTERIOR_PROB_THRESHOLD <- 0.8

parser <- argparser::arg_parser('Compile results')
parser <- argparser::add_argument(parser, "--studies_to_process", help = "Studies to process", type = "character")
parser <- argparser::add_argument(parser, "--studies_processed", help = "Current state of processed studies", type = "character")
parser <- argparser::add_argument(parser, "--coloc_result_files", help = "Coloc result files", type = "character", nargs = Inf)
parser <- argparser::add_argument(parser, "--compiled_results_file", help = "Compiled result file to save", type = "character")
args <- argparser::parse_args(parser)

main <- function(args) {
  all_studies_processed <- update_processed_study_metadata(args$studies_to_process, args$studies_processed)

  coloc_result_files <- split_string_into_vector(args$coloc_result_files)
  compiled_results <- compile_coloc_results(coloc_result_files, all_studies_processed)

  vroom::vroom_write(compiled_results, args$compiled_results_file)
}

update_processed_study_metadata <- function(studies_to_process_file, studies_processed_file) {
  studies_to_process <- vroom::vroom(studies_to_process_file, show_col_types=F)
  if (file.exists(studies_processed_file)) {
    studies_processed <- vroom::vroom(studies_processed_file, show_col_types=F)
    studies_processed <- rbind(studies_processed, studies_to_process) |> dplyr::distinct()
  } else {
    studies_processed <- studies_to_process
  }
  vroom::vroom_write(studies_processed, studies_processed_file)
}

compile_coloc_results <- function(coloc_result_files, studies_processed) {
  significant_results <- lapply(coloc_result_files, function(result_file) {
    if (!file.exists(result_file)) return()
    result <- vroom::vroom(result_file, delim='\t', show_col_types=F)
    #TODO: remove the last %in% names check once the bug has been fixed
    #if (is.null(result) || nrow(result) == 0 || !'posterior_prob' %in% names(result)) return ()

    significant_result <- dplyr::filter(result, !is.na(posterior_prob) & posterior_prob >= POSTERIOR_PROB_THRESHOLD)
    if (nrow(significant_result) > 0) return(significant_result)
  }) |> dplyr::bind_rows()

  pairwise_significant_results <- apply(significant_results, 1, function(result) {
    traits <- strsplit(result['traits'], ', ')[[1]]
    traits <- order_trait_by_type(traits, studies_processed)
    #TODO: also remove this now that we've changed coloc script
    #if (length(traits) < 2) {
    #  return()
    #}
    paired_results <- data.frame(t(utils::combn(traits, 2))) |>
      dplyr::rename(trait_a = X1, trait_b = X2) |>
      tidyr::separate(col = 'trait_a', into = c('study_a', 'ancestry_a', 'chr_a', 'bp_a', 'finemap_version_a'), sep='_', remove = F) |>
      tidyr::separate(col = 'trait_b', into = c('study_b', 'ancestry_b', 'chr_b', 'bp_b', 'finemap_version_b'), sep='_', remove = F) |>
      dplyr::mutate(posterior_prob = as.numeric(result['posterior_prob']),
                    candidate_snp = as.character(result['candidate_snp']),
                    posterior_explained_by_snp = as.numeric(result['posterior_explained_by_snp'])
      )
    paired_results$data_type_a <- studies_processed$data_type[match(paired_results$study_a, studies_processed$study_name)]
    paired_results$data_type_b <- studies_processed$data_type[match(paired_results$study_b, studies_processed$study_name)]

    #TODO: also remove this, like above
    #paired_results <- dplyr::filter(paired_results, study_a != study_b)

    return(paired_results)
  }) |> dplyr::bind_rows()
  return(pairwise_significant_results)
}

order_trait_by_type <- function(traits, studies_processed) {
  trait_studies <- unlist(lapply(traits, function(trait) strsplit(trait, '_')[[1]][1]))
  studies <- dplyr::filter(studies_processed, study_name %in% trait_studies) |>
    dplyr::select(study_name, data_type)
  traits <- traits[order(match(studies$data_type, ordered_data_types))]
  return(traits)
}

split_string_into_vector <- function(input_string) {
  return(unlist(strsplit(input_string, '[ ]')))
}

main(args)
