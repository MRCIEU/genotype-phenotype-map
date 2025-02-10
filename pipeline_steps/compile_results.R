source('constants.R')

posterior_prob_threshold <- 0.5

parser <- argparser::arg_parser('Compile results from pipeline')
#INPUT
parser <- argparser::add_argument(parser, '--studies_to_process', help = 'Studies to process', type = 'character')
parser <- argparser::add_argument(parser, '--studies_processed', help = 'Current state of processed studies', type = 'character')
#OUTPUT
parser <- argparser::add_argument(parser, '--all_study_blocks_file', help = 'Compiled result file to save', type = 'character')
parser <- argparser::add_argument(parser, '--raw_coloc_results_file', help = 'Raw coloc result files to amalgamate', type = 'character')
parser <- argparser::add_argument(parser, '--coloc_results_file', help = 'Compiled result file to save', type = 'character')
parser <- argparser::add_argument(parser, '--rare_results_file', help = 'Compiled result file to save', type = 'character')
parser <- argparser::add_argument(parser, '--compiled_results_metadata_file', help = 'Compiled result metadata file to save', type = 'character')
parser <- argparser::add_argument(parser, '--variant_annotations_file', help = 'Variant Annotations of Candidate SNPs', type = 'character')
parser <- argparser::add_argument(parser, '--pipeline_summary_file', help = 'Rendered Rmd file of output', type = 'character')

args <- argparser::parse_args(parser)

main <- function() {
  ld_blocks <- vroom::vroom('data/ld_blocks.tsv', show_col_types = F)
  ld_info <- construct_ld_block(ld_blocks$ancestry, ld_blocks$chr, ld_blocks$start, ld_blocks$stop) |>
    dplyr::filter(dir.exists(ld_block_data))

  pipeline_data <- aggregate_data_produced_by_pipeline(ld_info, args$studies_to_process, args$studies_processed)
  # cleanup_studies_with_no_extractions(pipeline_data)

  coloc_results <- compile_coloc_results(pipeline_data)
  rare_results <- compile_rare_results(pipeline_data)
  variant_annotations <- annotate_variants(pipeline_data)
  pipeline_data$all_study_blocks <- compile_study_blocks(pipeline_data)
  results_metadata <- aggregate_pipeline_metadata(pipeline_data, ld_info)

  validate_results(pipeline_data, variant_annotations, coloc_results)

  vroom::vroom_write(pipeline_data$raw_coloc_results, args$raw_coloc_results_file)
  vroom::vroom_write(rare_results, args$rare_results_file)
  vroom::vroom_write(coloc_results, args$coloc_results_file)
  vroom::vroom_write(variant_annotations, args$variant_annotations_file)
  vroom::vroom_write(pipeline_data$all_study_blocks, args$all_study_blocks_file)
  vroom::vroom_write(results_metadata, args$compiled_results_metadata_file)
  #this should always be the last thing done in the step, as we want to be able to rerun the pipeline other things fail
  vroom::vroom_write(pipeline_data$studies_processed, args$studies_processed)
  file.copy(args$studies_processed, dirname(args$coloc_results))

  if (is.na(TEST_RUN)) {
    rmarkdown::render("pipeline_summary.Rmd", output_file = args$pipeline_summary)
  } else {
    vroom::vroom_write(data.frame(), args$pipeline_summary_file)
  }
}

aggregate_data_produced_by_pipeline <- function(ld_info, studies_to_process_file, studies_processed_file) {
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
    dplyr::filter(!is.na(traits) & traits != 'None')

  compare_rare_input_files <- Filter(function(file) file.exists(file), glue::glue('{ld_info$ld_block_data}/compare_rare_results.tsv'))
  if (length(compare_rare_input_files) == 0) {
    compare_rare_results <- data.frame()
  } else {
    compare_rare_results <- vroom::vroom(compare_rare_input_files, delim='\t', show_col_types = F)
  }

  #update studies_processed.tsv with studies_to_process.tsv
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

  return(list(
    extracted_studies = extracted_studies,
    standardised_studies = standardised_studies,
    imputed_studies = imputed_studies,
    finemapped_studies = finemapped_studies,
    raw_coloc_results = raw_coloc_results,
    rare_results = compare_rare_results,
    studies_processed = studies_processed
  ))
}

compile_study_blocks <- function(pipeline_data) {
  finemapped_studies <- pipeline_data$finemapped_studies |>
    dplyr::filter(min_p <= p_value_threshold) |>
    dplyr::select(study, unique_study_id, file, chr, bp, min_p, cis_trans, ld_block)

  finemapped_studies$known_gene <- pipeline_data$studies_processed$gene[match(finemapped_studies$study, pipeline_data$studies_processed$study_name)]

  # rare_studies <- pipeline_data$standardised_studies |>
  #   dplyr::filter(variant_type != variant_types$common) |>
  #   dplyr::select(study, file, chr, bp, min_p, cis_trans, ld_block) |>
  #   dplyr::mutate(unique_study_id = NA, known_gene = if ('GENE' %in% colnames(rare_results)) rare_results$GENE else NA)

  # all_studies <- rbind(finemapped_studies, rare_studies)
  return(finemapped_studies)
}

# TODO: add all finemapped lead snps to this list in the future
annotate_variants <- function(pipeline_data) {
  coloc_candidate_snps <- unique(pipeline_data$raw_coloc_results$candidate_snp)
  variant_annotations <- vroom::vroom(glue::glue('{variant_annotation_dir}/vep_annotations_hg38.tsv.gz'), show_col_types = F) |>
    dplyr::filter(SNP %in% coloc_candidate_snps)
  
  return(variant_annotations)
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

#TODO: what should be checked here?
ingested_data_integrity_check <- function() {
  ld_blocks <- vroom::vroom('data/ld_blocks.tsv')
  ld_info <- construct_ld_block(ld_blocks$ancestry, ld_blocks$chr, ld_blocks$start, ld_blocks$stop)

  lapply(ld_info$ld_block_data, function(ld_block) {
    extracted_studies_file <- glue::glue('{ld_block}/extracted_studies.tsv')
    imputed_studies_file <- glue::glue('{ld_block}/imputed_studies.tsv')
    finemapped_studies_file <- glue::glue('{ld_block}/finemapped_studies.tsv')
    if (!file.exists(extracted_studies_file)) return()

    extracted_studies <- vroom::vroom(extracted_studies_file, show_col_types = F)
    imputed_studies <- vroom::vroom(extracted_studies_file, show_col_types = F)
    finemapped_studies <- vroom::vroom(extracted_studies_file, show_col_types = F)

  })
}

compile_rare_results <- function(pipeline_data) {
  pairwise_results <- apply(pipeline_data$rare_results, 1, function(result) {
    traits <- strsplit(result[['traits']], ', ')[[1]]
    ordered_traits <- order_trait_by_type(traits, pipeline_data$studies_processed)
    if (length(ordered_traits$unique_study_id) < 2) {
      return()
    }
    paired_results <- data.frame(t(utils::combn(ordered_traits$unique_study_id, 2))) |>
      dplyr::rename(unique_study_a = X1, unique_study_b = X2) |>
      dplyr::mutate(candidate_snp = as.character(result['candidate_snp']))

    #this figures out if each paired result is 'directed', meaning if the relationship is from earlier in the biological
    #causal pathway, as defined by ordered_data_types (ie. gene_expression -> protein -> phenotype)
    #this might have to get more complicated, as relationships between types is not always easy to define
    first_study_data_types <- ordered_traits$data_type[match(paired_results$unique_study_a, ordered_traits$unique_study_id)]
    second_study_data_types <- ordered_traits$data_type[match(paired_results$unique_study_b, ordered_traits$unique_study_id)]
    paired_results$directed <- first_study_data_types != second_study_data_types & second_study_data_types == ordered_data_types$phenotype

    return(paired_results)
  }) |>
    dplyr::bind_rows() |>
    dplyr::mutate(study_a = sub('.*_', '', unique_study_a), study_b = sub('.*_', '', unique_study_b))

  return(pairwise_results)
}

compile_coloc_results <- function(pipeline_data) {
  significant_results <- dplyr::filter(pipeline_data$raw_coloc_results, !is.na(traits) & !is.na(posterior_prob) & posterior_prob >= posterior_prob_threshold)

  # num_parallel_jobs <- 10
  # pairwise_significant_results <- parallel::mclapply(X=1:nrow(significant_results), mc.cores=num_parallel_jobs, FUN=function(i) {
  #   result <- significant_results[i, ]
  pairwise_significant_results <- apply(significant_results, 1, function(result) {
    traits <- strsplit(result[['traits']], ', ')[[1]]
    ordered_traits <- order_trait_by_type(traits, pipeline_data$studies_processed)
    if (length(ordered_traits$unique_study_id) < 2) {
      return()
    }

    paired_results <- data.frame(t(utils::combn(ordered_traits$unique_study_id, 2))) |>
      dplyr::rename(unique_study_a = X1, unique_study_b = X2) |>
      dplyr::mutate(posterior_prob = as.numeric(result['posterior_prob']),
                    candidate_snp = as.character(result['candidate_snp']),
                    posterior_explained_by_snp = as.numeric(result['posterior_explained_by_snp'])
      )

    #this figures out if each paired result is 'directed', meaning if the relationship is from earlier in the biological
    #causal pathway, as defined by ordered_data_types (ie. gene_expression -> protein -> phenotype)
    #this might have to get more complicated, as relationships between types is not always easy to define
    first_study_data_types <- ordered_traits$data_type[match(paired_results$unique_study_a, ordered_traits$unique_study_id)]
    second_study_data_types <- ordered_traits$data_type[match(paired_results$unique_study_b, ordered_traits$unique_study_id)]
    paired_results$directed <- first_study_data_types != second_study_data_types & second_study_data_types == ordered_data_types$phenotype

    return(paired_results)
  }) |> dplyr::bind_rows()

  #remove duplicate rows (of either study_a, study_b or study_b, study_a)
  cols <- c('unique_study_a','unique_study_b')
  pairwise_significant_results <- pairwise_significant_results[!duplicated(t(apply(pairwise_significant_results[cols], 1, sort))), ]

  #remove duplicates where the same 2 studies are colocalising on the same candidate SNP
  duplicate_candidate_snps <- data.frame(
    study_a=sub('_.*', '', pairwise_significant_results$unique_study_a),
    study_b=sub('_.*', '', pairwise_significant_results$unique_study_b),
    candidate_snp=pairwise_significant_results$candidate_snp
  )
  same_candidate_snp_duplicates <- duplicated(duplicate_candidate_snps)
  pairwise_significant_results <- pairwise_significant_results[!same_candidate_snp_duplicates,]

  return(pairwise_significant_results)
}

order_trait_by_type <- function(traits, studies_processed) {
  trait_studies <- sub('_.*', '', traits)
  studies <- dplyr::filter(studies_processed, study_name %in% trait_studies) |>
    dplyr::select(study_name, data_type)

  studies <- studies[order(studies$study_name), ]
  studies$unique_study_id <- traits[order(traits)]
  studies <- studies[order(match(studies$data_type, ordered_data_types)), ]
  return(studies)
}

split_string_into_vector <- function(input_string) {
  return(unlist(strsplit(input_string, '[ ]')))
}

validate_results <- function(pipeline_data, variant_annotations, coloc_results, rare_results) {
  #check that all unique ids in coloc results are in finemapped studies
  all_unique_ids <- unique(c(coloc_results$unique_study_a, coloc_results$unique_study_b))
  missing_unique_ids <- setdiff(all_unique_ids, pipeline_data$finemapped_studies$unique_study_id)
  if (length(missing_unique_ids) > 0) {
    message('Error: there are ', length(missing_unique_ids), ' unique study ids in coloc results are not in finemapped studies')
    message(paste(missing_unique_ids, collapse = ', '))
  }

  #check that all SNPs in coloc results are in the variant annotations 
  all_candidate_snps <- unique(coloc_results$candidate_snp)
  missing_candidate_snps <- setdiff(all_candidate_snps, variant_annotations$SNP)
  if (length(missing_candidate_snps) > 0) {
    message('Error: there are ', length(missing_candidate_snps), ' candidate SNPs in coloc results are not in variant annotations')
    message(paste(missing_candidate_snps, collapse = ', '))
  } 

  #check that all studies in rare results are in studies processed
  # all_unique_studies <- unique(c(rare_results$study_a, rare_results$study_b))
  # missing_studies <- setdiff(all_unique_studies, pipeline_data$studies_processed$study_name)
  # if (length(missing_studies) > 0) {
  #   message('Error: there are', length(missing_studies), 'studies in rare results are not in studies processed')
  #   stop(missing_studies)
  # }

  #check that all study ids in finemapped studies are in studies processed
  all_study_ids <- unique(pipeline_data$finemapped_studies$study)
  missing_study_ids <- setdiff(all_study_ids, pipeline_data$studies_processed$study_name)
  if (length(missing_study_ids) > 0) {
    message('Error: there are', length(missing_study_ids), 'study ids in finemapped studies are not in studies processed')
    message(paste(missing_study_ids, collapse = ', '))
  }
}

cleanup_studies_with_no_extractions <- function(pipeline_data) {
  study_dirs  <- Sys.glob(glue::glue('{extracted_study_dir}/*'))
  empty_study_dirs <- Filter(function(e) file.size(glue::glue('{e}/extracted_snps.tsv')) == 0, study_dirs)
  message('Studies with no extractions that will be cleaned up: ', length(empty_study_dirs))
  for (empty_study in empty_study_dirs) {
    system(glue::glue('rm -r {empty_study}'))
  }
}

main()