source('constants.R')

POSTERIOR_PROB_THRESHOLD <- 0.8

parser <- argparser::arg_parser('Compile results from pipeline')
#INPUT
parser <- argparser::add_argument(parser, '--studies_to_process', help = 'Studies to process', type = 'character')
parser <- argparser::add_argument(parser, '--studies_processed', help = 'Current state of processed studies', type = 'character')
parser <- argparser::add_argument(parser, '--coloc_input_files', help = 'All files from the perform coloc pipeline step', type = 'character', nargs = Inf)
#OUTPUT
parser <- argparser::add_argument(parser, '--all_study_regions_file', help = 'Compiled result file to save', type = 'character')
parser <- argparser::add_argument(parser, '--coloc_results_file', help = 'Compiled result file to save', type = 'character')
parser <- argparser::add_argument(parser, '--compiled_results_metadata_file', help = 'Compiled result metadata file to save', type = 'character')

args <- argparser::parse_args(parser)

main <- function(args) {
  all_studies_processed <- update_processed_study_metadata(args$studies_to_process, args$studies_processed)
  coloc_results <- compile_coloc_results(args$coloc_input_files, all_studies_processed)

  all_study_regions <- compile_entire_list_of_extracted_study_regions()
  results_metadata <- aggregate_pipeline_metadata()

  vroom::vroom_write(coloc_results, args$coloc_results_file)
  vroom::vroom_write(all_study_regions, args$all_study_regions_file)
  vroom::vroom_write(results_metadata, args$compiled_results_metadata_file)
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

compile_entire_list_of_extracted_study_regions <- function() {
  ld_regions <- vroom::vroom('data/ld_regions.tsv', show_col_types = F)
  ld_block_dirs <- paste0(ld_block_data_dir, ld_regions$pop, '/', ld_regions$chr, '/', ld_regions$start, '_', ld_regions$stop, '/')

  all_finemapped_studies <- lapply(ld_block_dirs, function(ld_block_dir) {
    finemapped_studies <- vroom::vroom(paste0(ld_block_dir, 'finemapped_studies.tsv'), show_col_types = F)
  }) |>
    dplyr::bind_rows() |>
    dplyr::select(study, file, chr, bp)

  #TODO: populate both 'known' (from coloc results) and 'suspected 'genes associated here?
  return(all_finemapped_studies)
}

aggregate_pipeline_metadata <- function() {
  ld_regions <- vroom::vroom('data/ld_regions.tsv', show_col_types = F)
  ld_block_dirs <- paste0(ld_block_data_dir, ld_regions$pop, '/', ld_regions$chr, '/', ld_regions$start, '_', ld_regions$stop, '/')

  lapply(ld_block_dirs, function(ld_block_dir) {
    if (!file.exists(paste0(ld_block_dir, 'extracted_studies.tsv'))) {
      return()
    }
    extracted_studies <- vroom::vroom(paste0(ld_block_dir, 'extracted_studies.tsv'), show_col_types = F)
    imputed_studies <- vroom::vroom(paste0(ld_block_dir, 'imputed_studies.tsv'), show_col_types = F)
    finemapped_studies <- vroom::vroom(paste0(ld_block_dir, 'finemapped_studies.tsv'), show_col_types = F)
    deduplicated_finemapped_studies <- finemapped_studies[!duplicated(finemapped_studies$study),]

    return(data.frame(ld_region = sub(, ld_block_dir),
                      extracted_studies=nrow(extracted_studies),
                      studies_imputed=nrow(imputed_studies),
                      mean_snps_imputed=mean(imputed_studies$rows_imputed),
                      number_finemapped=nrow(finemapped_studies),
                      finemap_success=sum(deduplicated_finemapped_studies$message == 'success'),
                      finemap_faile=sum(deduplicated_finemapped_studies$message == 'failed'),
                      finemap_no_need=sum(deduplicated_finemapped_studies$message == 'less_than_2_cs')
    ))
  }) |> dplyr::bind_rows()
}

compile_coloc_results <- function(coloc_input_files, studies_processed) {
  significant_results <- lapply(coloc_input_files, function(coloc_file) {
    coloc <- vroom::vroom(coloc_file)
    if (is.null(coloc)) return ()
    significant_result <- dplyr::filter(coloc, !is.na(posterior_prob) & posterior_prob >= POSTERIOR_PROB_THRESHOLD)
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
                    posterior_explained_by_snp = as.numeric(result['posterior_explained_by_snp']),
                    cis_trans = NA
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
