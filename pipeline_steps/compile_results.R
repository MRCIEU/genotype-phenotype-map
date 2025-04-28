source('constants.R')

posterior_prob_threshold <- 0.5

parser <- argparser::arg_parser('Compile results from pipeline')
#INPUT
parser <- argparser::add_argument(parser, '--studies_to_process', help = 'Studies to process', type = 'character')
parser <- argparser::add_argument(parser, '--studies_processed', help = 'Current state of processed studies', type = 'character')
parser <- argparser::add_argument(parser, '--traits_processed', help = 'Current state of processed traits', type = 'character')
#OUTPUT
parser <- argparser::add_argument(parser, '--new_studies_processed_file', help = 'New studies processed file to save', type = 'character')
parser <- argparser::add_argument(parser, '--new_traits_processed_file', help = 'New traits processed file to save', type = 'character')
parser <- argparser::add_argument(parser, '--study_extractions_file', help = 'Compiled result file to save', type = 'character')
parser <- argparser::add_argument(parser, '--coloc_results_file', help = 'Compiled result file to save', type = 'character')
parser <- argparser::add_argument(parser, '--rare_results_file', help = 'Compiled result file to save', type = 'character')
parser <- argparser::add_argument(parser, '--compiled_results_metadata_file', help = 'Compiled result metadata file to save', type = 'character')

args <- argparser::parse_args(parser)

main <- function() {
  ld_blocks <- vroom::vroom('data/ld_blocks.tsv', show_col_types = F)
  ld_info <- construct_ld_block(ld_blocks$ancestry, ld_blocks$chr, ld_blocks$start, ld_blocks$stop) |>
    dplyr::filter(dir.exists(ld_block_data))

  pipeline_data <- aggregate_data_produced_by_pipeline(ld_info, args$studies_to_process, args$studies_processed, args$traits_processed)

  pipeline_data$study_extractions <- create_study_extractions(pipeline_data)
  results_metadata <- aggregate_pipeline_metadata(pipeline_data, ld_info)

  vroom::vroom_write(pipeline_data$studies_processed, args$new_studies_processed_file)
  vroom::vroom_write(pipeline_data$traits_processed, args$new_traits_processed_file)
  vroom::vroom_write(pipeline_data$raw_coloc_results, args$coloc_results_file)
  vroom::vroom_write(pipeline_data$raw_rare_results, args$rare_results_file)
  vroom::vroom_write(pipeline_data$study_extractions, args$study_extractions_file)
  vroom::vroom_write(results_metadata, args$compiled_results_metadata_file)

  # if (is.na(TEST_RUN)) {
  #   rmarkdown::render("pipeline_summary.Rmd", output_file = args$pipeline_summary)
  # } else {
  #   vroom::vroom_write(data.frame(), args$pipeline_summary_file)
  # }
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
  imputed_studies <- vroom::vroom(imputed_studies_files, show_col_types = F, col_types = imputed_column_types)

  finemapped_studies_files <- Filter(function(file) file.exists(file), glue::glue('{ld_info$ld_block_data}/finemapped_studies.tsv'))
  finemapped_studies <- vroom::vroom(finemapped_studies_files, show_col_types = F, col_types = finemapped_column_types)

  coloc_input_files <- Filter(function(file) file.exists(file), glue::glue('{ld_info$ld_block_data}/coloc_results.tsv'))
  raw_coloc_results <- vroom::vroom(coloc_input_files, delim='\t', show_col_types = F) |>
    dplyr::filter(!is.na(traits) & traits != 'None' & posterior_prob >= posterior_prob_threshold)

  compare_rare_input_files <- Filter(function(file) file.exists(file), glue::glue('{ld_info$ld_block_data}/compare_rare_results.tsv'))
  if (length(compare_rare_input_files) == 0) {
    compare_rare_results <- data.frame()
  } else {
    compare_rare_results <- vroom::vroom(compare_rare_input_files, delim='\t', show_col_types = F)
  }

  #update studies_processed.tsv.gz with studies_to_process.tsv
  gene_name_map <- vroom::vroom(glue::glue('{liftover_dir}/gene_name_map.tsv'), show_col_types=F)

  studies_to_process <- vroom::vroom(studies_to_process_file, show_col_types=F)
  gene_names <- gene_name_map$GENE_NAME[match(studies_to_process$gene, gene_name_map$ENSEMBL_ID)]
  gene_names[is.na(gene_names)] <- studies_to_process$gene[is.na(gene_names)]
  studies_to_process$gene <- gene_names

  if (file.exists(studies_processed_file)) {
    studies_processed <- vroom::vroom(studies_processed_file, show_col_types=F)
    studies_processed <- rbind(studies_processed, studies_to_process) |> dplyr::distinct()
  } else {
    studies_processed <- studies_to_process
  }

  traits_unprocessed <- studies_processed |>
    dplyr::filter(!study_name %in% traits_processed$study_name) |>
    dplyr::select(data_type, study_name, trait_name)
  traits_processed <- vroom::vroom(traits_processed_file, show_col_types=F)

  traits_processed <- rbind(traits_processed, traits_unprocessed) |> dplyr::distinct()

  return(list(
    extracted_studies = extracted_studies,
    standardised_studies = standardised_studies,
    imputed_studies = imputed_studies,
    finemapped_studies = finemapped_studies,
    raw_coloc_results = raw_coloc_results,
    raw_rare_results = compare_rare_results,
    studies_processed = studies_processed,
    traits_processed = traits_processed
  ))
}

create_study_extractions <- function(pipeline_data) {
  finemapped_studies <- pipeline_data$finemapped_studies |>
    dplyr::filter(min_p <= p_value_threshold) |>
    dplyr::select(study, unique_study_id, file, chr, bp, min_p, cis_trans, ld_block)

  finemapped_studies$known_gene <- pipeline_data$studies_processed$gene[match(finemapped_studies$study, pipeline_data$studies_processed$study_name)]

  standardised_studies <- pipeline_data$standardised_studies |>
    dplyr::filter(variant_type != variant_types$common) |>
    dplyr::mutate(unique_study_id = paste0(study, '_', chr, '_', bp))
  
  rare_studies <- pipeline_data$raw_rare_results |>
    dplyr::mutate(candidate_snp=trimws(candidate_snp)) |>
    tidyr::separate_rows(traits, min_ps, genes, sep=", ") |>
    dplyr::rename(unique_study_id=traits, min_p=min_ps, known_gene=genes) |>
    dplyr::mutate(min_p = as.numeric(min_p)) |>
    dplyr::left_join(standardised_studies, by="unique_study_id") |>
    dplyr::select(study, unique_study_id, file, chr, bp, min_p, cis_trans, ld_block, known_gene)

  all_studies <- rbind(finemapped_studies, rare_studies)
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