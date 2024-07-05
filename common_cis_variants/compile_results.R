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
  compiled_results <- compile_coloc_results(args$coloc_result_files, all_studies_processed)
  vroom::vroom_write(compiled_results, args$compiled_results_file)
}

update_processed_study_metadata <- function(studies_to_process, studies_processed) {
  gwases_to_extract <- vroom::vroom(studies_to_process)
  if (file.exists(studies_processed)) {
    studies_processed <- vroom::vroom(studies_processed)
    studies_processed <- rbind(studies_processed, gwases_to_extract) |> dplyr::distinct()
  } else {
    studies_processed <- gwases_to_extract
  }
  vroom::vroom_write(studies_processed, studies_processed)
}

compile_coloc_results <- function(coloc_result_files, studies_processed) {
  significant_results <- lapply(coloc_result_files, function(result_file) {
    result <- vroom::vroom(result_file)
    if (is.null(result)) return ()
    significant_result <- dplyr::filter(result, !is.na(posterior_prob) & posterior_prob >= POSTERIOR_PROB_THRESHOLD)
    if (nrow(significant_result) > 0) return(significant_result)
  }) |> dplyr::bind_rows()

  pairwise_significant_results <- apply(significant_results, 1, function(result) {
    traits <- strsplit(result['traits'], ', ')[[1]]
    traits <- order_trait_by_type(traits, studies_processed)

    paired_results <- data.frame(t(utils::combn(traits, 2))) |>
      dplyr::rename(trait_a = X1, trait_b = X2) |>
      tidyr::separate(col = 'trait_a', into = c('study_a', 'chr_a', 'bp_b', 'finemap_version_a'), sep='_', remove = F) |>
      tidyr::separate(col = 'trait_b', into = c('study_b', 'chr_b', 'bp_b', 'finemap_version_b'), sep='_', remove = F) |>
      dplyr::mutate(posterior_prob = as.numeric(result['posterior_prob']),
                    candidate_snp = as.character(result['candidate_snp']),
                    posterior_explained_by_snp = as.numeric(result['posterior_explained_by_snp'],)
      )
    paired_results$data_type_a <- studies_processed$data_type[match(paried_results$study_a, studies_processed$study)]
    paired_results$data_type_b <- studies_processed$data_type[match(paired_results$study_b, studies_processed$study)]
  }) |> dplyr::bind_rows()
  return(pairwise_significant_results)
}

order_trait_by_type <- function(traits, studies_processed) {
  studies <- dplyr::filter(studies_processed, study %in% traits) |>
    dplyr::select(study, data_type)
  traits <- studies$study[order(match(studies$data_type, ordered_data_types))]
  return(traits)
}

main(args)