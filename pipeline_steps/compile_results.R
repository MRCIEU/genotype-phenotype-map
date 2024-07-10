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

  all_study_regions <- compile_entire_list_of_extracted_study_regions(all_studies_processed)
  results_metadata <- aggregate_pipeline_metadata()

  vroom::vroom_write(coloc_results, args$coloc_results_file)
  vroom::vroom_write(all_study_regions, args$all_study_regions_file)
  vroom::vroom_write(results_metadata, args$compiled_results_metadata_file)
  vroom::vroom_write(all_studies_processed, args$studies_processed)
  file.copy(args$studies_processed, dirname(args$coloc_results))
}

update_processed_study_metadata <- function(studies_to_process_file, studies_processed_file) {
  studies_to_process <- vroom::vroom(studies_to_process_file, show_col_types=F)
  if (file.exists(studies_processed_file)) {
    studies_processed <- vroom::vroom(studies_processed_file, show_col_types=F)
    studies_processed <- rbind(studies_processed, studies_to_process) |> dplyr::distinct()
  } else {
    studies_processed <- studies_to_process
  }
  return(studies_processed)
}

compile_entire_list_of_extracted_study_regions <- function(all_studies) {
  ld_regions <- vroom::vroom('data/ld_regions.tsv', show_col_types = F)
  ld_info <- construct_ld_block(ld_regions$ancestry, ld_regions$chr, ld_regions$start, ld_regions$stop)

  all_finemapped_studies <- lapply(ld_info$ld_block_data, function(ld_block_dir) {
    finemap_study <- paste0(ld_block_dir, 'finemapped_studies.tsv')
    if (!file.exists(finemap_study)) return(data.frame())
    finemapped_studies <- vroom::vroom(finemap_study, show_col_types = F) |>
      dplyr::select(study, unique_study_id, file, chr, bp, cis_trans) |>
      dplyr::mutate(chr = as.character(chr), bp = as.numeric(bp))
  }) |>
    dplyr::bind_rows()

  all_finemapped_studies$known_gene <- all_studies$gene[match(all_finemapped_studies$study, all_studies$study_name)]
  all_finemapped_studies <- find_suspected_gene_associated_with_position(all_finemapped_studies)
  return(all_finemapped_studies)
}

find_suspected_gene_associated_with_position <- function(all_finemapped_studies) {
  coords <- GenomicFeatures::genes(TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene)
  gene_names <- as.data.frame(org.Hs.eg.db::org.Hs.egSYMBOL)

  genomic_ranges <- dplyr::select(all_finemapped_studies, chrom=chr, start=bp, end=bp, unique_study_id=unique_study_id) |>
    dplyr::mutate(chrom=paste0('chr', chrom)) |>
    GenomicRanges::makeGRangesFromDataFrame(na.rm = T, keep.extra.columns = T)
  genes <- IRanges::mergeByOverlaps(coords, genomic_ranges)

  genes$gene_name <- gene_names$symbol[match(genes$gene_id, gene_names$gene_id)]
  all_finemapped_studies$suspected_gene <- genes$gene_name[match(all_finemapped_studies$unique_study_id, genes$unique_study_id)]
  return(all_finemapped_studies)
}

aggregate_pipeline_metadata <- function() {
  ld_regions <- vroom::vroom('data/ld_regions.tsv', show_col_types = F)
  ld_block_dirs <- paste0(ld_block_data_dir, ld_regions$ancestry, '/', ld_regions$chr, '/', ld_regions$start, '_', ld_regions$stop, '/')

  metadata_per_ld_region <- lapply(ld_block_dirs, function(ld_block_dir) {
    if (!file.exists(paste0(ld_block_dir, 'finemapped_studies.tsv'))) {
      return(data.frame())
    }
    extracted_studies <- vroom::vroom(paste0(ld_block_dir, 'extracted_studies.tsv'), show_col_types = F)
    imputed_studies <- vroom::vroom(paste0(ld_block_dir, 'imputed_studies.tsv'), show_col_types = F)
    finemapped_studies <- vroom::vroom(paste0(ld_block_dir, 'finemapped_studies.tsv'), show_col_types = F)

    return(data.frame(ld_region = ld_block_dir,
                      extracted_regions=nrow(extracted_studies),
                      studies_imputed=nrow(imputed_studies),
                      mean_snps_imputed=mean(imputed_studies$rows_imputed),
                      number_finemapped=nrow(finemapped_studies),
                      finemap_failed=sum(finemapped_studies$message == 'failed'),
                      finemap_no_need=sum(finemapped_studies$message == 'less_than_2_cs')
    ))
  }) |> dplyr::bind_rows()
  
  means <- colMeans(metadata_per_ld_region[-1])
  means$ld_region <- 'mean'
  totals <- colSums(metadata_per_ld_region[-1])
  totals$ld_region <- 'total'
  metadata_per_ld_region <- dplyr::bind_rows(metadata_per_ld_region, means, totals)
  return(metadata_per_ld_region)
}

compile_coloc_results <- function(coloc_input_files, studies_processed) {
  significant_results <- lapply(coloc_input_files, function(coloc_file) {
    coloc <- vroom::vroom(coloc_file, show_col_types = F)
    if (is.null(coloc) || nrow(coloc) == 0) return ()
    significant_result <- dplyr::filter(coloc, !is.na(posterior_prob) & posterior_prob >= POSTERIOR_PROB_THRESHOLD)
    if (nrow(significant_result) > 0) return(significant_result)
  }) |> dplyr::bind_rows()

  pairwise_significant_results <- apply(significant_results, 1, function(result) {
    traits <- strsplit(result['traits'], ', ')[[1]]
    ordered_traits <- order_trait_by_type(traits, studies_processed)
    if (nrow(ordered_traits) < 2) {
      return()
    }

    paired_results <- data.frame(t(utils::combn(ordered_traits$unique_study_id, 2))) |>
      dplyr::rename(unique_study_a = X1, unique_study_b = X2) |>
      dplyr::mutate(posterior_prob = as.numeric(result['posterior_prob']),
                    candidate_snp = as.character(result['candidate_snp']),
                    posterior_explained_by_snp = as.numeric(result['posterior_explained_by_snp']),
                    study_a = ordered_traits$study_name[match(unique_study_a, ordered_traits$unique_study_id)],
                    study_b = ordered_traits$study_name[match(unique_study_b, ordered_traits$unique_study_id)]
      )

    #this figures out if each paired result is 'directed', meaning if the relationship is from earlier in the biological
    #causal pathway, as defined by ordered_data_types (ie. gene_expression -> protein -> phenotype)
    #this might have to get more complicated, as relationships between types is not always easy to define
    first_study_data_types <- ordered_traits$data_type[match(paired_results$unique_study_a, ordered_traits$unique_study_id)]
    second_study_data_types <- ordered_traits$data_type[match(paired_results$unique_study_b, ordered_traits$unique_study_id)]
    paired_results$directed <- first_study_data_types != second_study_data_types & second_data_types == ordered_data_types$phenotype

    return(paired_results)
  }) |> dplyr::bind_rows()
  return(pairwise_significant_results)
}

order_trait_by_type <- function(traits, studies_processed) {
  traits <- sort(traits)
  trait_studies <- unlist(lapply(traits, function(trait) strsplit(trait, '_')[[1]][1]))
  studies <- dplyr::filter(studies_processed, study_name %in% trait_studies) |>
    dplyr::select(study_name, data_type) |>
    dplyr::arrange(study_name)
  studies$unique_study_id <- traits
  studies <- studies[order(match(studies$data_type, ordered_data_types)), ]
  return(studies)
}

split_string_into_vector <- function(input_string) {
  return(unlist(strsplit(input_string, '[ ]')))
}

main(args)
