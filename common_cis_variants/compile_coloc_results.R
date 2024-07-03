source('constants.R')

POSTERIOR_PROB_THRESHOLD <- 0.8

parser <- argparser::arg_parser('Compile results')
parser <- argparser::add_argument(parser, "--studies_to_process", help = "Studies to process", type = "character")
parser <- argparser::add_argument(parser, "--current_state", help = "Current state of processed studies", type = "character")
parser <- argparser::add_argument(parser, "--coloc_results", help = "Coloc result files", type = "character", nargs = Inf)
parser <- argparser::add_argument(parser, "--compiled_results_file", help = "Compiled result file to save", type = "character")
args <- argparser::parse_args(parser)

main <- function(args) {
  #Interesting example: ukb-b-7953_10_12300790
  gwases_to_extract <- vroom::vroom(args$studies_to_process)
  if (file.exists(args$current_state)) {
    current_state <- vroom::vroom(args$current_state)
    current_state <- rbind(current_state, gwases_to_extract) |> dplyr::distinct()
  } else {
    current_state <- gwases_to_extract
  }
  vroom::vroom_write(current_state, args$current_state)

  significant_results <- lapply(args$coloc_results, function(result_file) {
    result <- vroom::vroom(result_file)
    if (is.null(result)) return ()
    significant_result <- dplyr::filter(result, posterior_prob >= POSTERIOR_PROB_THRESHOLD)
    if (nrow(significant_result) > 0) return(significant_result)
  }) |> dplyr::bind_rows()

  pairwise_significant_results <- apply(significant_results, 1, function(result) {
    traits <- strsplit(result['traits'], ', ')[[1]]
    traits <- order_trait_by_type(traits, current_state)

    paired_results <- data.frame(t(utils::combn(traits, 2))) |>
      dplyr::rename(trait_a = X1, trait_b = X2) |>
      dplyr::mutate(posterior_prob = as.numeric(result['posterior_prob']),
                    candidate_snp = as.character(result['candidate_snp']),
                    posterior_explained_by_snp = as.numeric(result['posterior_explained_by_snp'])
      )
    #TODO: maybe add associated gene to paired_results, by finding current_state$associated_gene with trait_a on study name?
  }) |> dplyr::bind_rows()
}

order_trait_by_type <- function(traits, current_state) {
  studies <- dplyr::filter(current_state, study %in% traits) |>
    dplyr::select(study, data_type)
  traits <- studies$study[order(match(studies$data_type, ordered_data_types))]
  return(traits)
}


main(args)