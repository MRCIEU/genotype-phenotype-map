source('constants.R')

POSTERIOR_PROB_THRESHOLD <- 0.5

parser <- argparser::arg_parser('Compile results from pipeline')
#INPUT
parser <- argparser::add_argument(parser, '--studies_to_process', help = 'Studies to process', type = 'character')
parser <- argparser::add_argument(parser, '--studies_processed', help = 'Current state of processed studies', type = 'character')
#OUTPUT
parser <- argparser::add_argument(parser, '--all_study_blocks_file', help = 'Compiled result file to save', type = 'character')
parser <- argparser::add_argument(parser, '--raw_coloc_results_file', help = 'Raw coloc result files to amalgamate', type = 'character')
parser <- argparser::add_argument(parser, '--coloc_results_file', help = 'Compiled result file to save', type = 'character')
parser <- argparser::add_argument(parser, '--compiled_results_metadata_file', help = 'Compiled result metadata file to save', type = 'character')

args <- argparser::parse_args(parser)


main <- function(args) {
  ld_blocks <- vroom::vroom('data/ld_blocks.tsv', show_col_types = F)
  ld_info <- construct_ld_block(ld_blocks$ancestry, ld_blocks$chr, ld_blocks$start, ld_blocks$stop)

  coloc_input_files <- glue::glue('{ld_info$ld_block_results}/coloc_results.tsv')
  coloc_input_files <- Filter(function(file) file.exists(file), coloc_input_files)

  raw_coloc_results <- vroom::vroom(coloc_input_files, delim='\t', show_col_types = F) |>
    dplyr::filter(!is.na(traits) & traits != 'None')

  all_studies_processed <- update_processed_study_metadata(args$studies_to_process, args$studies_processed)
  coloc_results <- compile_coloc_results(raw_coloc_results, all_studies_processed)

  all_study_blocks <- compile_entire_list_of_extracted_study_regions(all_studies_processed, ld_info)
  results_metadata <- aggregate_pipeline_metadata(ld_info)

  vroom::vroom_write(raw_coloc_results, args$raw_coloc_results_file)
  vroom::vroom_write(coloc_results, args$coloc_results_file)
  vroom::vroom_write(all_study_blocks, args$all_study_blocks_file)
  vroom::vroom_write(results_metadata, args$compiled_results_metadata_file)
  vroom::vroom_write(all_studies_processed, args$studies_processed)
  file.copy(args$studies_processed, dirname(args$coloc_results))
}

update_processed_study_metadata <- function(studies_to_process_file, studies_processed_file) {
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
  return(studies_processed)
}

compile_entire_list_of_extracted_study_regions <- function(all_studies, ld_info) {
  all_finemapped_studies <- apply(ld_info, 1, function(ld_block) {
    finemap_study <- glue::glue('{ld_block["ld_block_data"]}/finemapped_studies.tsv')
    if (!file.exists(finemap_study) || file.size(finemap_study) == 0L) return(data.frame())

    finemapped_studies <- vroom::vroom(finemap_study, delim = '\t', show_col_types = F) |>
      dplyr::select(study, unique_study_id, file, chr, bp, min_p, cis_trans) |>
      dplyr::mutate(chr = as.character(chr), bp = as.numeric(bp), min_p = as.numeric(min_p))
    finemapped_studies$ld_block <- ld_block['block']
    return(finemapped_studies)
  }) |> dplyr::bind_rows()

  all_finemapped_studies$known_gene <- all_studies$gene[match(all_finemapped_studies$study, all_studies$study_name)]
  all_finemapped_studies <- find_suspected_gene_associated_with_position(all_finemapped_studies)
  return(all_finemapped_studies)
}

find_suspected_gene_associated_with_position <- function(all_finemapped_studies) {
  coords <- GenomicFeatures::genes(TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene)
  gene_names <- as.data.frame(org.Hs.eg.db::org.Hs.egSYMBOL)

  genomic_ranges <- dplyr::select(all_finemapped_studies, chrom=chr, start=bp, end=bp, unique_study_id=unique_study_id) |>
    dplyr::mutate(chrom=glue::glue('chr{chrom}')) |>
    GenomicRanges::makeGRangesFromDataFrame(na.rm = T, keep.extra.columns = T)
  genes <- IRanges::mergeByOverlaps(coords, genomic_ranges)

  genes$gene_name <- gene_names$symbol[match(genes$gene_id, gene_names$gene_id)]
  all_finemapped_studies$suspected_gene <- genes$gene_name[match(all_finemapped_studies$unique_study_id, genes$unique_study_id)]
  return(all_finemapped_studies)
}

aggregate_pipeline_metadata <- function(ld_info) {
  metadata_per_ld_block <- apply(ld_info, 1, function(ld) {
    ld_block_data <- ld['ld_block_data']
    if (!file.exists(glue::glue('{ld_block_data}/finemapped_studies.tsv'))) {
      return(data.frame())
    }
    extracted_studies <- vroom::vroom(glue::glue('{ld_block_data}/extracted_studies.tsv'), show_col_types = F)
    standardised_studies <- vroom::vroom(glue::glue('{ld_block_data}/standardised_studies.tsv'), show_col_types = F)
    imputed_studies <- vroom::vroom(glue::glue('{ld_block_data}/imputed_studies.tsv'), show_col_types = F)
    finemapped_studies <- vroom::vroom(glue::glue('{ld_block_data}/finemapped_studies.tsv'), show_col_types = F)
    above_threshold <- nrow(dplyr::filter(finemapped_studies, min_p <= p_value_threshold))

    return(data.frame(ld_block = ld['block'],
                      extracted_regions=nrow(extracted_studies),
                      standardised_time_taken=mean(as.difftime(standardised_studies$time_taken), na.rm=T),
                      mean_snps_removed_by_reference_panel=mean(standardised_studies$snps_removed_by_reference_panel, na.rm=T),
                      studies_imputed=nrow(imputed_studies),
                      imputed_time_taken=mean(as.difftime(imputed_studies$time_taken), na.rm=T),
                      mean_snps_imputed=mean(imputed_studies$rows_imputed, na.rm=T),
                      number_finemapped=nrow(finemapped_studies),
                      number_finemapped_above_threshold=above_threshold,
                      finemapped_time_taken=mean(as.difftime(finemapped_studies$time_taken), na.rm=T),
                      finemap_failed=sum(finemapped_studies$message == 'failed', na.rm=T),
                      finemap_no_need=sum(finemapped_studies$message == 'less_than_2_cs', na.rm=T)
    ))
  }) |> dplyr::bind_rows()

  # means <- colMeans(metadata_per_ld_block[-1])
  # means$ld_block <- 'mean'
  # totals <- colSums(metadata_per_ld_block[-1])
  # totals$ld_block <- 'total'
  # metadata_per_ld_block <- dplyr::bind_rows(metadata_per_ld_block, means, totals)
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

compile_coloc_results <- function(raw_coloc_results, studies_processed) {
  significant_results <- dplyr::filter(raw_coloc_results, !is.na(traits) & !is.na(posterior_prob) & posterior_prob >= POSTERIOR_PROB_THRESHOLD)

  pairwise_significant_results <- apply(significant_results, 1, function(result) {
    traits <- strsplit(result[['traits']], ', ')[[1]]
    ordered_traits <- order_trait_by_type(traits, studies_processed)
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

main(args)
