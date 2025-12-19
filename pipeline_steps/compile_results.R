source('constants.R')

parser <- argparser::arg_parser('Compile results from pipeline')
#INPUT
parser <- argparser::add_argument(parser, '--studies_to_process', help = 'Studies to process', type = 'character')
parser <- argparser::add_argument(parser, '--studies_processed', help = 'Current state of processed studies', type = 'character')
parser <- argparser::add_argument(parser, '--traits_processed', help = 'Current state of processed traits', type = 'character')
#OUTPUT
parser <- argparser::add_argument(parser, '--new_studies_processed_file', help = 'New studies processed file to save', type = 'character')
parser <- argparser::add_argument(parser, '--new_traits_processed_file', help = 'New traits processed file to save', type = 'character')
parser <- argparser::add_argument(parser, '--study_extractions_file', help = 'Compiled result file to save', type = 'character')
parser <- argparser::add_argument(parser, '--coloc_pairwise_results_file', help = 'Compiled result file to save', type = 'character')
parser <- argparser::add_argument(parser, '--coloc_clustered_results_file', help = 'Compiled result file to save', type = 'character')
parser <- argparser::add_argument(parser, '--rare_results_file', help = 'Compiled result file to save', type = 'character')
parser <- argparser::add_argument(parser, '--compiled_results_metadata_file', help = 'Compiled result metadata file to save', type = 'character')

args <- argparser::parse_args(parser)

main <- function() {
  ld_blocks <- vroom::vroom('data/ld_blocks.tsv', show_col_types = F)
  ld_info <- construct_ld_block(ld_blocks$ancestry, ld_blocks$chr, ld_blocks$start, ld_blocks$stop) |>
    dplyr::filter(dir.exists(ld_block_data))

  message('Aggregating data produced by pipeline')
  pipeline_data <- aggregate_data_produced_by_pipeline(ld_info, args$studies_to_process, args$studies_processed, args$traits_processed)

  message('Creating study extractions')
  pipeline_data$study_extractions <- create_study_extractions(pipeline_data)

  # message('Aggregating pipeline metadata')
  # results_metadata <- aggregate_pipeline_metadata(pipeline_data, ld_info)

  message('Writing results')
  vroom::vroom_write(pipeline_data$studies_processed, args$new_studies_processed_file)
  vroom::vroom_write(pipeline_data$traits_processed, args$new_traits_processed_file)
  vroom::vroom_write(pipeline_data$study_extractions, args$study_extractions_file)
  vroom::vroom_write(pipeline_data$rare_results, args$rare_results_file)
  vroom::vroom_write(pipeline_data$coloc_clustered_results, args$coloc_clustered_results_file)
  vroom::vroom_write(pipeline_data$coloc_pairwise_results, args$coloc_pairwise_results_file)
  # vroom::vroom_write(results_metadata, args$compiled_results_metadata_file)
}

aggregate_data_produced_by_pipeline <- function(ld_info, studies_to_process_file, studies_processed_file, traits_processed_file) {
  #take all files in each ld_block, and concatenate the data into one dataframe each (for use in the rest of this script)
  extracted_studies_files <- Filter(function(file) file.exists(file), glue::glue('{ld_info$ld_block_data}/extracted_studies.tsv'))
  extracted_studies <- vroom::vroom(extracted_studies_files, show_col_types = F)

  standardised_studies_files <- Filter(function(file) file.exists(file),
    glue::glue('{ld_info$ld_block_data}/standardised_studies.tsv')
  )
  standardised_studies <- vroom::vroom(standardised_studies_files, show_col_types = F, col_types = standardised_column_types)

  imputed_studies_files <- Filter(function(file) file.exists(file), glue::glue('{ld_info$ld_block_data}/imputed_studies.tsv'))
  imputed_studies <- lapply(imputed_studies_files, function(file) {
    vroom::vroom(file, show_col_types = F, col_types = imputed_column_types)
  }) |> dplyr::bind_rows()

  finemapped_studies_files <- Filter(function(file) file.exists(file), glue::glue('{ld_info$ld_block_data}/finemapped_studies.tsv'))
  finemapped_studies <- vroom::vroom(finemapped_studies_files, show_col_types = F, col_types = finemapped_column_types)

  pairwise_coloc_input_files <- Filter(function(file) file.exists(file), glue::glue('{ld_info$ld_block_data}/coloc_pairwise_results.tsv.gz'))
  coloc_pairwise_results <- vroom::vroom(pairwise_coloc_input_files,
    delim='\t',
    show_col_types = F,
    col_select = c('unique_study_a', 'unique_study_b', 'PP.H0.abf', 'PP.H1.abf', 'PP.H2.abf', 'PP.H3.abf', 'PP.H4.abf', 'ld_block', 'false_positive', 'false_negative', 'ignore')
  )
  message('pairwise coloc results: ', nrow(coloc_pairwise_results))

  coloc_clustered_input_files <- Filter(function(file) file.exists(file), glue::glue('{ld_info$ld_block_data}/coloc_clustered_results.tsv.gz'))
  coloc_clustered_results <- vroom::vroom(coloc_clustered_input_files, delim='\t', show_col_types = F)
  message('clustered coloc results: ', nrow(coloc_clustered_results))

  coloc_clustered_results <- coloc_clustered_results |>
    dplyr::group_by(ld_block, component) |>
    dplyr::mutate(coloc_group_id = dplyr::cur_group_id()) |>
    dplyr::ungroup() |>
    dplyr::arrange(coloc_group_id)

  compare_rare_input_files <- Filter(function(file) file.exists(file), glue::glue('{ld_info$ld_block_data}/compare_rare_results.tsv'))
  if (length(compare_rare_input_files) == 0) {
    compare_rare_results <- data.frame()
  } else {
    compare_rare_results <- vroom::vroom(compare_rare_input_files, delim='\t', show_col_types = F)
  }

  studies_to_process <- vroom::vroom(studies_to_process_file, show_col_types = F)

  if (file.exists(studies_processed_file)) {
    studies_processed <- vroom::vroom(studies_processed_file, show_col_types=F)
    studies_processed <- dplyr::bind_rows(studies_processed, studies_to_process) |> dplyr::distinct()
  } else {
    studies_processed <- studies_to_process
  }

  if (file.exists(traits_processed_file)) {
    traits_processed <- vroom::vroom(traits_processed_file, show_col_types=F)
  } else {
    traits_processed <- data.frame(data_type=c(), study_name=c(), trait=c())
  }

  traits_unprocessed <- studies_processed |>
    dplyr::filter(!study_name %in% traits_processed$study_name & !trait %in% traits_processed$study_name) |>
    dplyr::select(data_type, study_name, trait_name) |>
    dplyr::rename(trait=trait_name)

  traits_processed <- dplyr::bind_rows(traits_processed, traits_unprocessed) |> dplyr::distinct()
  if (!"category" %in% colnames(traits_processed)) {
    traits_processed <- dplyr::mutate(traits_processed, category = NA) 
  }

  return(list(
    extracted_studies = extracted_studies,
    standardised_studies = standardised_studies,
    imputed_studies = imputed_studies,
    finemapped_studies = finemapped_studies,
    coloc_pairwise_results = coloc_pairwise_results,
    coloc_clustered_results = coloc_clustered_results,
    rare_results = compare_rare_results,
    studies_processed = studies_processed,
    traits_processed = traits_processed
  ))
}

create_study_extractions <- function(pipeline_data) {
  finemapped_studies <- pipeline_data$finemapped_studies |>
    dplyr::filter(min_p <= p_value_threshold) |>
    dplyr::select(study, unique_study_id, file, snp, chr, bp, min_p, cis_trans, ld_block, svg_file, file_with_lbfs, ignore)

  finemapped_studies$known_gene <- pipeline_data$studies_processed$gene[match(finemapped_studies$study, pipeline_data$studies_processed$study_name)]

  rare_genes_study_map <- pipeline_data$studies_processed |>
    dplyr::filter(variant_type != variant_types$common & !is.na(gene)) |>
    dplyr::select(study_name, gene) |>
    dplyr::rename(known_gene=gene)

  if (nrow(pipeline_data$rare_results) > 0) {
    rare_studies <- pipeline_data$rare_results |>
      dplyr::mutate(candidate_snp=trimws(candidate_snp)) |>
      tidyr::separate_rows(traits, min_ps, genes, files, sep=", ") |>
      dplyr::rename(unique_study_id=traits, min_p=min_ps, situated_gene=genes, file=files, snp=candidate_snp) |>
      dplyr::mutate(min_p = as.numeric(min_p)) |>
      tidyr::separate(unique_study_id, into = c("study", "ancestry", "chr", "bp"), sep = "_", remove = F) |>
      dplyr::mutate(bp = as.numeric(sub('-.*', '', bp))) |>
      dplyr::mutate(file_with_lbfs = NA, svg_file = NA, ignore = F) |>
      dplyr::mutate(
        known_gene = rare_genes_study_map$known_gene[match(study, rare_genes_study_map$study_name)],
        cis_trans = dplyr::case_when(
          is.na(known_gene) | is.na(situated_gene) ~ NA_character_,
          known_gene == situated_gene ~ "cis",
          TRUE ~ "trans"
        )
      ) |>
      dplyr::select(study, unique_study_id, file, snp, chr, bp, min_p, cis_trans, ld_block, file_with_lbfs, svg_file, ignore, known_gene, situated_gene)
  } else {
    rare_studies <- data.frame()
  }

  all_studies <- dplyr::bind_rows(finemapped_studies, rare_studies)
  return(all_studies)
}

aggregate_pipeline_metadata <- function(pipeline_data, ld_info) {
  metadata_per_ld_block <- lapply(ld_info$block, function(block) {
    extracted_per_block <- dplyr::filter(pipeline_data$extracted_studies, ld_block == block)
    standardised_per_block <- dplyr::filter(pipeline_data$standardised_studies, ld_block == block)
    imputed_per_block <- dplyr::filter(pipeline_data$imputed_studies, ld_block == block)
    finemapped_per_block <- dplyr::filter(pipeline_data$finemapped_studies, ld_block == block)

    unique_finemapped_per_block <- dplyr::filter(finemapped_per_block, grepl('_1$', unique_study_id))

    return(data.frame(ld_block = block,
                      number_extracted=nrow(extracted_per_block),
                      number_standardised=nrow(standardised_per_block),
                      mean_snps_removed_by_reference_panel=mean(standardised_per_block$snps_removed_by_reference_panel, na.rm=T),
                      number_imputed=nrow(imputed_per_block),
                      significant_snps_imputed=mean(imputed_per_block$significant_rows_imputed, na.rm=T),
                      significant_imputed_snps_filtered=mean(imputed_per_block$significant_rows_filtered, na.rm=T),
                      number_finemapped=nrow(dplyr::filter(finemapped_per_block, min_p <= p_value_threshold)),
                      finemapped_per_imputed=nrow(finemapped_per_block) / nrow(unique_finemapped_per_block),
                      num_finemap_failed=sum(finemapped_per_block$finemap_message == 'failed', na.rm=T),
                      standardised_time_taken=mean(as.difftime(standardised_per_block$time_taken), na.rm=T),
                      imputed_time_taken=mean(as.difftime(imputed_per_block$time_taken), na.rm=T),
                      finemapped_time_taken=mean(as.difftime(finemapped_per_block$time_taken), na.rm=T)
    ))
  }) |> dplyr::bind_rows()

  return(metadata_per_ld_block)
}

main()