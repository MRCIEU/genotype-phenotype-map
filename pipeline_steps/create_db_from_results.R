source('constants.R')
source('database_definitions.R')
source('gwas_calculations.R')

parser <- argparser::arg_parser('Create DuckDB from pipeline results')
parser <- argparser::add_argument(parser, '--studies_db_file', help = 'Studies DB file', type = 'character')
parser <- argparser::add_argument(parser, '--associations_db_file', help = 'Associations DB file', type = 'character')
parser <- argparser::add_argument(parser, '--coloc_pairs_full_db_file', help = 'Coloc pairs DB file', type = 'character')
parser <- argparser::add_argument(parser, '--coloc_pairs_significant_db_file', help = 'Coloc pairs DB file', type = 'character')
parser <- argparser::add_argument(parser, '--ld_db_file', help = 'LD DB file', type = 'character')
parser <- argparser::add_argument(parser, '--gwas_upload_db_file', help = 'GWAS upload DB file', type = 'character')
parser <- argparser::add_argument(parser, '--completed_output_file', help = 'Completed output file', type = 'character')

args <- argparser::parse_args(parser)

main <- function() {
  if (file.exists(args$studies_db_file)) file.remove(args$studies_db_file)
  if (file.exists(args$associations_db_file)) file.remove(args$associations_db_file)
  if (file.exists(args$coloc_pairs_full_db_file)) file.remove(args$coloc_pairs_full_db_file)
  if (file.exists(args$coloc_pairs_significant_db_file)) file.remove(args$coloc_pairs_significant_db_file)
  if (file.exists(args$ld_db_file)) file.remove(args$ld_db_file)
  if (file.exists(args$gwas_upload_db_file)) file.remove(args$gwas_upload_db_file)

  latest_studies_conn <- duckdb::dbConnect(duckdb::duckdb(), glue::glue("{latest_results_dir}/studies.db"), read_only = TRUE)
  studies_conn <- duckdb::dbConnect(duckdb::duckdb(), args$studies_db_file)
  associations_conn <- duckdb::dbConnect(duckdb::duckdb(), args$associations_db_file)
  coloc_pairs_full_conn <- duckdb::dbConnect(duckdb::duckdb(), args$coloc_pairs_full_db_file)
  coloc_pairs_significant_conn <- duckdb::dbConnect(duckdb::duckdb(), args$coloc_pairs_significant_db_file)
  ld_conn <- duckdb::dbConnect(duckdb::duckdb(), args$ld_db_file)
  gwas_upload_conn <- duckdb::dbConnect(duckdb::duckdb(), args$gwas_upload_db_file)

  lapply(studies_db, \(table) {
    DBI::dbExecute(studies_conn, table$query)
    if (!is.null(table$indexes)) DBI::dbExecute(studies_conn, table$indexes)
  })

  lapply(gwas_upload_db, \(table) DBI::dbExecute(gwas_upload_conn, table$query))
  DBI::dbExecute(coloc_pairs_full_conn, coloc_pairs_full_table$query)
  DBI::dbExecute(coloc_pairs_significant_conn, coloc_pairs_significant_table$query)
  DBI::dbExecute(ld_conn, ld_table$query)

  message('Creating studies db')
  studies_db <- populate_existing_row_ids(latest_studies_conn, studies_db)
  studies_db <- load_data_for_studies_db(studies_db, studies_conn)

  all_relevant_snps <- find_relevant_snps(studies_db)

  message('Creating coloc pairs db...')
  load_data_into_coloc_pairs_db(coloc_pairs_full_conn, coloc_pairs_significant_conn, studies_db)

  message('Creating ld db...')
  load_data_into_ld_db(ld_conn, studies_db, all_relevant_snps)
  match_causal_snps_in_high_ld(ld_conn, studies_conn, studies_db)

  message('Creating associations db...')
  load_data_into_associations_db(associations_conn, studies_db, all_relevant_snps)

  DBI::dbDisconnect(latest_studies_conn, shutdown=TRUE)
  DBI::dbDisconnect(studies_conn, shutdown=TRUE)
  DBI::dbDisconnect(associations_conn, shutdown=TRUE)
  DBI::dbDisconnect(coloc_pairs_full_conn, shutdown=TRUE)
  DBI::dbDisconnect(coloc_pairs_significant_conn, shutdown=TRUE)
  DBI::dbDisconnect(ld_conn, shutdown=TRUE)
  DBI::dbDisconnect(gwas_upload_conn, shutdown=TRUE)

  vroom::vroom_write(data.frame(), args$completed_output_file)
}

populate_existing_row_ids <- function(conn, tables) {
  for (table_name in names(tables)) {
    table <- tables[[table_name]]
    if (!is.na(table$persist_id_from)) {
      existing_ids <- DBI::dbGetQuery(conn, glue::glue("SELECT id, {table$persist_id_from} FROM {table_name}"))
      tables[[table_name]]$existing_ids <- existing_ids
    }
  }
  return(tables)
}

populate_missing_row_ids <- function(table, id_column_name) {
  max_id <- suppressWarnings(max(table[[id_column_name]], na.rm = TRUE))
  if (!is.finite(max_id)) max_id <- 0

  na_idx <- which(is.na(table[[id_column_name]]))
  if (length(na_idx) > 0) {
    table[[id_column_name]][na_idx] <- seq(from = max_id + 1, length.out = length(na_idx))
  }
  return(table)
}

# Adds existing ids to the table if they exist, and assigns new sequential ids to rows where id is NA
resolve_ids_for_table <- function(table, existing_ids=NULL, join_by) {
  if (is.null(existing_ids)) return(table)

  table <- table |> dplyr::left_join(existing_ids, by=join_by)

  max_id <- suppressWarnings(max(table$id, na.rm = TRUE))
  if (!is.finite(max_id)) max_id <- 0

  na_idx <- which(is.na(table$id))
  if (length(na_idx) > 0) {
    table$id[na_idx] <- seq(from = max_id + 1, length.out = length(na_idx))
  }
  table <- table[order(table$id), ]
  return(table)
}


load_data_for_studies_db <- function(studies_db, studies_conn) {
  studies_db$study_sources$data <- vroom::vroom(file.path("data/study_sources.csv"), show_col_types = F) |>
    resolve_ids_for_table(studies_db$study_sources$existing_ids, studies_db$study_sources$persist_id_from) |>
    dplyr::select(get_table_column_names(studies_db$study_sources))

  studies_db$ld_blocks$data <- vroom::vroom(file.path("data/ld_blocks.tsv"), show_col_types = F) |>
    dplyr::mutate(ld_block=paste0(ancestry, "/", chr, "/", start, "-", stop)) |>
    resolve_ids_for_table(studies_db$ld_blocks$existing_ids, studies_db$ld_blocks$persist_id_from) |>
    dplyr::select(get_table_column_names(studies_db$ld_blocks))

  studies_db$gene_annotations$data <- vroom::vroom(file.path(variant_annotation_dir, "gene_info.tsv"), show_col_types = F) |>
    resolve_ids_for_table(studies_db$gene_annotations$existing_ids, studies_db$gene_annotations$persist_id_from) |>
    dplyr::select(get_table_column_names(studies_db$gene_annotations))

  studies_db$study_extractions$data <- vroom::vroom(file.path(current_results_dir, "study_extractions.tsv.gz"), show_col_types = F)

  sources_subset <- studies_db$study_sources$data |>
    dplyr::select(source, id) |>
    dplyr::rename(source_id=id)

  gene_subset <- studies_db$gene_annotations$data |>
    dplyr::select(ensembl_id, id) |>
    dplyr::rename(gene_id=id)
  
  # Remove the studies that don't have any study extractions
  studies_db$studies$data <- vroom::vroom(file.path(current_results_dir, "studies_processed.tsv.gz"), show_col_types = F) |>
    dplyr::left_join(sources_subset, by=c("source"="source")) |>
    dplyr::left_join(gene_subset, by=c("ensg"="ensembl_id")) |>
    dplyr::filter(study_name %in% studies_db$study_extractions$data$study & trait %in% studies_db$study_extractions$data$study) |>
    resolve_ids_for_table(studies_db$studies$existing_ids, studies_db$studies$persist_id_from)

  # Remove the traits that don't have any study extractions
  studies_db$traits$data <- vroom::vroom(file.path(current_results_dir, "traits_processed.tsv.gz"), show_col_types = F) |>
    dplyr::filter(study_name %in% studies_db$studies$data$study_name) |>
    resolve_ids_for_table(studies_db$traits$existing_ids, studies_db$traits$persist_id_from) |>
    dplyr::rename(trait_name=trait, trait=study_name, trait_category=category) |>
    dplyr::select(get_table_column_names(studies_db$traits))

  traits_subset <- studies_db$traits$data |>
    dplyr::select(trait, id) |>
    dplyr::rename(trait_id=id)
  
  studies_db$studies$data <- studies_db$studies$data |>
    dplyr::left_join(traits_subset, by=c("trait"="trait"))

  studies_db$studies$data <- studies_db$studies$data |>
    dplyr::select(get_table_column_names(studies_db$studies))
  
  gene_subset <- studies_db$gene_annotations$data |>
    dplyr::select(gene, id) |>
    dplyr::rename(gene_id=id)

  studies_db$snp_annotations$data <- vroom::vroom(file.path(variant_annotation_dir, "vep_annotations_hg38.tsv.gz"), show_col_types =  F) |>
    resolve_ids_for_table(studies_db$snp_annotations$existing_ids, studies_db$snp_annotations$persist_id_from) |>
    dplyr::left_join(gene_subset, by=c("gene"="gene")) |>
    dplyr::mutate(snp=trimws(snp)) |>
    dplyr::mutate(rsid=sub(",.*", "", rsid)) |>
    dplyr::select(get_table_column_names(studies_db$snp_annotations))

  studies_db$study_extractions$data <- studies_db$study_extractions$data |>
    format_study_extractions(studies_db)

  studies_db$coloc_groups$data <- vroom::vroom(file.path(current_results_dir, "coloc_clustered_results.tsv.gz"), show_col_types = F) |>
    format_clustered_colocs(studies_db)

  studies_db$rare_results$data <- vroom::vroom(file.path(current_results_dir, "rare_results.tsv.gz"), show_col_types = F) |>
    dplyr::rename_with(tolower) |>
    format_rare_results(studies_db)
  
  lapply(studies_db, \(table) DBI::dbAppendTable(studies_conn, table$name, table$data))
  create_wide_tables(studies_conn)

  return(studies_db)
}

format_study_extractions <- function(study_extractions, studies_db) {
  gene_subset <- studies_db$gene_annotations$data |>
    dplyr::select(gene, id) |>
    dplyr::rename(gene_id=id)

  snp_annotations_subset <- studies_db$snp_annotations$data |>
    dplyr::select(snp, display_snp, id) |>
    dplyr::rename(snp_id=id)
  
  studies_subset <- studies_db$studies$data |>
    dplyr::select(study_name, id) |>
    dplyr::rename(study_id=id)
  
  duplicate_studies <- studies_db$studies$data[duplicated(studies_db$studies$data$study_name), ]
  if (nrow(duplicate_studies) > 0) {
    vroom::vroom_write(duplicate_studies, file.path(current_results_dir, "duplicate_studies.tsv"))
    stop("ERROR: Found ", nrow(duplicate_studies), " rows with duplicate study_name values.  Saving to current results dir")
  }

  duplicate_snps <- studies_db$snp_annotations$data[duplicated(studies_db$snp_annotations$data$snp), ]
  if (nrow(duplicate_snps) > 0) {
    vroom::vroom_write(duplicate_snps, file.path(current_results_dir, "duplicate_snps.tsv"))
    stop("ERROR: Found ", nrow(duplicate_snps), " rows with duplicate display_snp values.  Saving to current results dir")
  }

  duplicate_genes <- studies_db$gene_annotations$data[duplicated(studies_db$gene_annotations$data$gene), ]
  if (nrow(duplicate_genes) > 0) {
    vroom::vroom_write(duplicate_genes, file.path(current_results_dir, "duplicate_genes.tsv"))
    stop("ERROR: Found ", nrow(duplicate_genes), " rows with duplicate gene values.  Saving to current results dir")
  }

  ld_blocks_subset <- studies_db$ld_blocks$data |>
    dplyr::select(ld_block, id) |>
    dplyr::rename(ld_block_id=id)

  study_extractions <- study_extractions |>
    dplyr::filter(ignore == F) |>
    resolve_ids_for_table(studies_db$study_extractions$existing_ids, studies_db$study_extractions$persist_id_from) |>
    dplyr::mutate(
      file=sub(ld_block_data_dir, "", file),
      svg_file=sub(data_dir, "", svg_file),
      file_with_lbfs=sub(data_dir, "", file_with_lbfs),
      snp=trimws(snp)
    ) |>
    dplyr::rename(gene=known_gene) |>
    dplyr::left_join(studies_subset, by=c("study"="study_name")) |>
    dplyr::left_join(ld_blocks_subset, by="ld_block") |>
    dplyr::left_join(gene_subset, by="gene") |>
    dplyr::left_join(snp_annotations_subset, by="snp") |>
    populate_missing_row_ids("id") |>
    dplyr::select(get_table_column_names(studies_db$study_extractions))

  duplicates <- study_extractions[duplicated(study_extractions$unique_study_id), ]
  
  if (nrow(duplicates) > 0) {
    vroom::vroom_write(duplicates, file.path(current_results_dir, "duplicate_study_extractions.tsv"))
    stop("ERROR: Found ", nrow(duplicates), " rows with duplicate unique_study_id values.  Saving to current results dir")
  }

  missing_snps <- study_extractions |>
    dplyr::filter(is.na(snp_id) | is.null(snp_id))
  
  if (nrow(missing_snps) > 0) {
    message("WARNING: Found ", nrow(missing_snps), " rows with missing snp_id values.  Removing study extractions with missing snp_id values")
    vroom::vroom_write(missing_snps, file.path(current_results_dir, "missing_snps_in_study_extractions.tsv"))
    study_extractions <- study_extractions[!is.na(study_extractions$snp_id) & !is.null(study_extractions$snp_id), ]
  }

  return(study_extractions)
}

format_clustered_colocs <- function(clustered_colocs, studies_db) {
  study_extractions_subset <- studies_db$study_extractions$data |>
    dplyr::select(id, study_id, unique_study_id, ld_block_id) |>
    dplyr::rename(study_extraction_id=id)

  snp_annotations_subset <- studies_db$snp_annotations$data |>
    dplyr::select(snp, id) |>
    dplyr::rename(snp_id=id)

  clustered_colocs <- clustered_colocs |>
    dplyr::arrange(coloc_group_id) |>
    dplyr::left_join(study_extractions_subset, by=c("unique_study_id"="unique_study_id")) |>
    dplyr::left_join(snp_annotations_subset, by="snp")
  
  missing_stuff <- clustered_colocs |>
    dplyr::filter(is.na(study_extraction_id) | is.na(snp_id))
  
  if (nrow(missing_stuff) > 0) {
    message("WARNING: Found ", nrow(missing_stuff), " rows with missing study_extraction_id or snp_id values.  Removing rows with missing data.")
    vroom::vroom_write(missing_stuff, file.path(current_results_dir, "missing_stuff_in_clustered_colocs.tsv"))
    clustered_colocs <- clustered_colocs[!is.na(clustered_colocs$study_extraction_id) & !is.na(clustered_colocs$snp_id), ]
  }

  clustered_colocs <- clustered_colocs |>
    dplyr::select(get_table_column_names(studies_db$coloc_groups)) |>
    dplyr::arrange(coloc_group_id)
  
  return(clustered_colocs)
}

format_rare_results <- function(rare_results, studies_db) {
  study_extractions_subset <- studies_db$study_extractions$data |>
    dplyr::select(id, unique_study_id, ld_block_id, study_id) |>
    dplyr::rename(study_extraction_id=id)

  snp_annotations_subset <- studies_db$snp_annotations$data |>
    dplyr::select(snp, id) |>
    dplyr::rename(snp_id=id)

  gene_annotations_subset <- studies_db$gene_annotations$data |>
    dplyr::select(gene, id) |>
    dplyr::rename(gene_id=id)
  
  rare_results <- rare_results |>
    dplyr::mutate(rare_result_group_id=1:dplyr::n()) |>
    tidyr::separate_rows(traits, genes, files, sep=", ") |>
    dplyr::rename(unique_study_id=traits, gene=genes) |>
    dplyr::mutate(candidate_snp=trimws(candidate_snp), unique_study_id=trimws(unique_study_id), gene=trimws(gene), files=trimws(files)) |>
    dplyr::left_join(study_extractions_subset, by="unique_study_id") |>
    dplyr::left_join(snp_annotations_subset, by=c("candidate_snp"="snp"), relationship="many-to-many") |>
    dplyr::left_join(gene_annotations_subset, by="gene", relationship="many-to-one")

  missing_study_extractions <- dplyr::filter(rare_results, is.na(study_extraction_id) | is.na(snp_id))
  if (nrow(missing_study_extractions) > 0) {
    message("WARNING: Found ", nrow(missing_study_extractions), " rows with missing study extractions.  Saving to current results dir")
    vroom::vroom_write(missing_study_extractions, file.path(current_results_dir, "missing_study_extractions_in_rare_results.tsv"))
    rare_results <- rare_results[!is.na(rare_results$study_extraction_id) & !is.na(rare_results$snp_id), ]
  }

  rare_results <- rare_results |>
    dplyr::select(get_table_column_names(studies_db$rare_results)) |>
    dplyr::arrange(rare_result_group_id)

  return(rare_results)
}

create_wide_tables <- function(studies_conn) {
  safe_lapply(additional_studies_tables, function(table) {
    DBI::dbExecute(studies_conn, table$query)
    DBI::dbExecute(studies_conn, table$indexes)
  })
}

load_data_into_coloc_pairs_db <- function(coloc_pairs_full_conn, coloc_pairs_significant_conn, studies_db) {
  pairwise_colocs <- data.table::fread(file.path(current_results_dir, "coloc_pairwise_results.tsv.gz"))
  pairwise_colocs <- pairwise_colocs[!is.na(PP.H4.abf) & ignore == FALSE]

  study_extractions_subset <- studies_db$study_extractions$data |>
    dplyr::select(id, unique_study_id, ld_block_id) |>
    dplyr::rename(study_extraction_id=id) |>
    data.table::as.data.table()
  
  data.table::setkey(pairwise_colocs, unique_study_a)
  data.table::setkey(study_extractions_subset, unique_study_id)
  pairwise_colocs <- study_extractions_subset[pairwise_colocs, on = "unique_study_id==unique_study_a"]
  data.table::setnames(pairwise_colocs, c("study_extraction_id", "ld_block_id"), c("study_extraction_a_id", "ld_block_id_a"))
  
  data.table::setkey(pairwise_colocs, unique_study_b)
  pairwise_colocs <- study_extractions_subset[pairwise_colocs, on = "unique_study_id==unique_study_b"]
  data.table::setnames(pairwise_colocs, c("study_extraction_id", "ld_block_id"), c("study_extraction_b_id", "ld_block_id_b"))
  data.table::setnames(pairwise_colocs, c("PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf"), c("h0", "h1", "h2", "h3", "h4"))
  
  pairwise_colocs[, ld_block_id := ld_block_id_a]
  pairwise_colocs[, c("ld_block_id_a", "ld_block_id_b") := NULL]

  missing_study_extractions <- pairwise_colocs[is.na(study_extraction_a_id) | is.na(study_extraction_b_id)]
  if (nrow(missing_study_extractions) > 0) {
    message("WARNING: Found ", nrow(missing_study_extractions), " rows with missing study extractions. Removing rows with missing data.")
    data.table::fwrite(missing_study_extractions, file.path(current_results_dir, "missing_study_extractions_in_pairwise_colocs.tsv"))
    pairwise_colocs <- pairwise_colocs[!is.na(study_extraction_a_id) & !is.na(study_extraction_b_id)]
  }

  pairwise_colocs <- pairwise_colocs |>
    dplyr::select(get_table_column_names(coloc_pairs_full_table)) |>
    dplyr::arrange(study_extraction_a_id, study_extraction_b_id)

  DBI::dbAppendTable(coloc_pairs_full_conn, coloc_pairs_full_table$name, pairwise_colocs)

  significant_colocs <- pairwise_colocs[h4 > posterior_prob_threshold]
  
  study_extraction_snp_map <- studies_db$coloc_groups$data |>
    dplyr::select(study_extraction_id, snp_id) |>
    dplyr::distinct()

  significant_colocs <- significant_colocs |>
    dplyr::left_join(study_extraction_snp_map, by=c("study_extraction_a_id"="study_extraction_id")) |>
    dplyr::filter(!is.na(snp_id)) |>
    dplyr::select(get_table_column_names(coloc_pairs_significant_table)) |>
    dplyr::arrange(snp_id, study_extraction_a_id, study_extraction_b_id)

  # DBI::dbAppendTable(coloc_pairs_significant_conn, coloc_pairs_significant_table$name, significant_colocs)

  significant_colocs_chunk_info <- split_large_dataframe_into_chunks(significant_colocs, "snp_id")
  DBI::dbExecute(coloc_pairs_significant_conn, coloc_pairs_significant_db$coloc_pairs_metadata$query)

  DBI::dbAppendTable(coloc_pairs_significant_conn, coloc_pairs_significant_db$coloc_pairs_metadata$name,
    data.frame(
      start_snp_id = significant_colocs_chunk_info$dataframe_metadata$start_snp_id,
      stop_snp_id = significant_colocs_chunk_info$dataframe_metadata$end_snp_id,
      coloc_pairs_table_name = glue::glue("coloc_pairs_{1:nrow(significant_colocs_chunk_info$dataframe_metadata)}")
    )
  )

  lapply(seq_along(significant_colocs_chunk_info$dataframe_chunks), function(i) {
    table_chunk_name <- glue::glue("coloc_pairs_{i}")
    create_table_query <- sub("table_name", table_chunk_name, coloc_pairs_significant_db$coloc_pairs$query)
    DBI::dbExecute(coloc_pairs_significant_conn, create_table_query)

    indexes_query <- gsub("table_name", table_chunk_name, coloc_pairs_significant_db$coloc_pairs$indexes)
    DBI::dbExecute(coloc_pairs_significant_conn, indexes_query)

    DBI::dbAppendTable(coloc_pairs_significant_conn, table_chunk_name, significant_colocs_chunk_info$dataframe_chunks[[i]])
    message('Added ', nrow(significant_colocs_chunk_info$dataframe_chunks[[i]]), ' rows to ', table_chunk_name)
  })
}

# Find relevant snps: all coloc SNPs, all rare association SNPs and all finemapped SNPs that colocalise with nothing (at genome wide significance)
find_relevant_snps <- function(studies_db) {
  snp_annotations_subset <- data.table::as.data.table(studies_db$snp_annotations$data)[, .(candidate_snp = snp, snp_id = id)]
  ld_blocks_subset <- data.table::as.data.table(studies_db$ld_blocks$data)[, .(ld_block, ld_block_id = id)]
  studies_subset <- data.table::as.data.table(studies_db$studies$data)[, .(variant_type, study_name, study_id = id)]

  colocalising_snps <- data.table::as.data.table(studies_db$coloc_groups$data)[, .(snp_id, ld_block_id)]
  colocalising_snps <- unique(colocalising_snps)
  colocalising_snps <- snp_annotations_subset[colocalising_snps, on = "snp_id"]
  colocalising_snps <- ld_blocks_subset[colocalising_snps, on = "ld_block_id"]
  colocalising_snps$variant_type <- variant_types$common
  colocalising_snps$study_name <- NA

  rare_results_snps <- data.table::as.data.table(studies_db$rare_results$data)[, .(snp_id, ld_block_id)]
  rare_results_snps <- unique(rare_results_snps)
  rare_results_snps <- snp_annotations_subset[rare_results_snps, on = "snp_id"]
  rare_results_snps <- ld_blocks_subset[rare_results_snps, on = "ld_block_id"]
  rare_results_snps$variant_type <- variant_types$rare_exome
  rare_results_snps$study_name <- NA

  study_extractions_snps <- data.table::as.data.table(studies_db$study_extractions$data)[, .(id, min_p, snp_id, study_id, ld_block_id)]
  study_extractions_snps <- unique(study_extractions_snps)
  study_extractions_snps <- study_extractions_snps[
    !id %in% studies_db$coloc_groups$data$study_extraction_id &
    min_p <= lowest_p_value_threshold
  ]
  study_extractions_snps <- snp_annotations_subset[study_extractions_snps, on = "snp_id"]
  study_extractions_snps <- ld_blocks_subset[study_extractions_snps, on = "ld_block_id"]
  study_extractions_snps <- studies_subset[study_extractions_snps, on = "study_id"]
  study_extractions_snps$study_id <- NULL
  study_extractions_snps$id <- NULL
  study_extractions_snps$min_p <- NULL

  message("Found ", nrow(colocalising_snps), " colocalising SNPs, ",
    nrow(rare_results_snps), " rare SNPs, ",
    nrow(study_extractions_snps), " significant non-colocalising SNPs"
  )
  relevant_snps <- dplyr::bind_rows(colocalising_snps, rare_results_snps, study_extractions_snps)
  vroom::vroom_write(relevant_snps, file.path(current_results_dir, "relevant_snps.tsv"))
  return(relevant_snps)
}

load_data_into_ld_db <- function(ld_conn, studies_db, all_relevant_snps) {
  all_relevant_snps <- all_relevant_snps[variant_type == variant_types$common & is.na(study_name)]
  relevant_ld_blocks <- unique(studies_db$study_extractions$data$ld_block)

  ld_data <- parallel::mclapply(relevant_ld_blocks, mc.cores=25, \(relevant_ld_block) {
    gc()
    tryCatch({
      relevant_snps <- data.table::as.data.table(all_relevant_snps)[ld_block == relevant_ld_block]
      generate_ld_obj(relevant_ld_block, relevant_snps$candidate_snp)
    }, error = function(e) {
      message("Error generating LD object for ", relevant_ld_block, " - ", e)
      throw(e)
    })
  }) |> data.table::rbindlist(fill = TRUE)

  variant_annotations_subset_lead <- data.table::as.data.table(studies_db$snp_annotations$data)[, .(snp, lead_flipped=flipped, lead_snp_id = id)]
  variant_annotations_subset_variant <- data.table::as.data.table(studies_db$snp_annotations$data)[, .(snp, variant_flipped=flipped, variant_snp_id = id)]
  ld_blocks_subset <- data.table::as.data.table(studies_db$ld_blocks$data)[, .(ld_block, ld_block_id = id)]

  ld_data <- variant_annotations_subset_lead[ld_data, on = c("snp" = "lead")]
  ld_data <- variant_annotations_subset_variant[ld_data, on = c("snp" = "variant")]
  ld_data <- ld_blocks_subset[ld_data, on = "ld_block"]

  ld_data$r <- ld_data$r * ifelse(ld_data$lead_flipped, -1, 1) * ifelse(ld_data$variant_flipped, -1, 1)

  ld_data <- ld_data[, .(lead_snp_id, variant_snp_id, ld_block_id, r)]
  ld_data <- ld_data[order(lead_snp_id, variant_snp_id)]

  DBI::dbAppendTable(ld_conn, ld_table$name, ld_data)
}

generate_ld_obj <- function(ld_block, snps) {
  message("Generating LD object for ", ld_block)
  ld_file <- file.path(ld_reference_panel_dir, glue::glue("{ld_block}.unphased.vcor1"))
  ld <- data.table::fread(ld_file, header = FALSE, showProgress = FALSE)
  ldvars <- data.table::fread(glue::glue("{ld_file}.vars"), sep = " ", header = FALSE, showProgress = FALSE)

  data.table::setnames(ld, ldvars$V1)
  ld[, lead := ldvars$V1]
  ind <- which(ldvars$V1 %in% snps)
  ld <- ld[ind]

  ld_long <- data.table::melt(ld, id.vars = "lead", variable.name = "variant", value.name = "r")
  ld_long <- ld_long[r^2 > 0.8 | variant %in% ld$lead]
  ld_long <- ld_long[lead != variant]
  ld_long[, ld_block := ld_block]
  
  return(ld_long)
}

match_causal_snps_in_high_ld <- function(ld_conn, studies_conn, studies_db) {
  r2_threshold <- 0.98
  causal_snps <- unique(DBI::dbGetQuery(studies_conn, "SELECT DISTINCT snp_id FROM coloc_groups"))

  formatted_snp_ids <- paste0(unique(causal_snps$snp_id), collapse = ",")
  snp_pair_ld_scores <- DBI::dbGetQuery(ld_conn, glue::glue("SELECT lead_snp_id, variant_snp_id
    FROM ld
    WHERE lead_snp_id IN ({formatted_snp_ids})
    AND variant_snp_id IN ({formatted_snp_ids})
    AND r*r > {r2_threshold}"
  ))

  # Ensure only unique, non-redundant pairs are used for the undirected graph
  ld_edges <- snp_pair_ld_scores |>
    dplyr::rowwise() |>
    dplyr::mutate(s1 = min(lead_snp_id, variant_snp_id), s2 = max(lead_snp_id, variant_snp_id)) |>
    dplyr::ungroup() |>
    dplyr::distinct(s1, s2) |>
    dplyr::select(s1, s2)

  snp_graph <- igraph::graph_from_data_frame(ld_edges, directed = FALSE)
  clusters <- igraph::components(snp_graph)
  snp_mapping_intermediate <- tibble::tibble(snp_id = as.integer(names(clusters$membership)), cluster_id = clusters$membership)

  snp_mapping <- snp_mapping_intermediate |>
    dplyr::group_by(cluster_id) |>
    dplyr::mutate(representative_snp_id = min(snp_id)) |>
    dplyr::ungroup() |>
    dplyr::select(snp_id, representative_snp_id)

  studies_db$coloc_groups$data <- studies_db$coloc_groups$data |>
    dplyr::left_join(snp_mapping, by = "snp_id") |>
    dplyr::mutate(snp_id = coalesce(representative_snp_id, snp_id)) |>
    dplyr::select(-representative_snp_id)

  # Update the DuckDB coloc_groups table
  DBI::dbExecute(studies_conn, "DROP TABLE coloc_groups")
  DBI::dbAppendTable(studies_conn, "coloc_groups", studies_db$coloc_groups$data)

  return(studies_db)
}

load_data_into_associations_db <- function(conn, studies_db, all_relevant_snps) {
  relevant_snps_per_ld_block <- split(all_relevant_snps, all_relevant_snps$ld_block)

  associations <- parallel::mclapply(names(relevant_snps_per_ld_block), mc.cores=30, \(ld_block) {
    gc()
    relevant_snps <- relevant_snps_per_ld_block[[ld_block]]
    specific_snps <- relevant_snps[!is.na(study_name)]
    general_snps <- relevant_snps[is.na(study_name)]
    message("Processing ld_block: ", ld_block, " - ", nrow(general_snps), " general SNPs, ", nrow(specific_snps), " study specific SNPs")

    associations <- extract_associations_for_ld_block(ld_block, general_snps, specific_snps, studies_db)
    return(associations)
  }) 

  associations <- associations[!sapply(associations, is.null)]
  associations <- data.table::rbindlist(associations, fill = TRUE)

  associations <- flip_alleles(associations, associations$flipped)
  associations <- associations[, .(snp_id, study_id, beta=BETA, se=SE, imputed, p, eaf=EAF)]
  associations <- unique(associations, by = c("snp_id", "study_id"))

  original_num_rows <- nrow(associations)
  associations <- associations[!is.na(beta) & !is.na(se) & !is.na(p) & !is.na(eaf)]
  message('Removed ', original_num_rows - nrow(associations), ' rows with missing values')

  associations <- associations[order(snp_id, study_id)]
  association_chunk_info <- split_large_dataframe_into_chunks(associations, "snp_id")

  DBI::dbExecute(conn, associations_db$associations_metadata$query)
  DBI::dbAppendTable(conn, associations_db$associations_metadata$name,
    data.frame(
      start_snp_id = association_chunk_info$dataframe_metadata$start_snp_id,
      stop_snp_id = association_chunk_info$dataframe_metadata$end_snp_id,
      associations_table_name = glue::glue("associations_{1:nrow(association_chunk_info$dataframe_metadata)}")
    )
  )

  lapply(seq_along(association_chunk_info$dataframe_chunks), function(i) {
    table_chunk_name <- glue::glue("associations_{i}")
    create_table_query <- sub("table_name", table_chunk_name, associations_db$associations$query)
    DBI::dbExecute(conn, create_table_query)
    DBI::dbAppendTable(conn, table_chunk_name, association_chunk_info$dataframe_chunks[[i]])
    message('Added ', nrow(association_chunk_info$dataframe_chunks[[i]]), ' rows to ', table_chunk_name)
  })
}

extract_associations_for_ld_block <- function(ld_block, general_snps, specific_snps, studies_db) {
  snp_annotations_subset <- data.table::as.data.table(studies_db$snp_annotations$data)[, .(snp, flipped, snp_id = id)]
  studies_subset <- data.table::as.data.table(studies_db$studies$data)[, .(study_name, study_id = id)]

  tryCatch({
    imputed_studies_file <- file.path(ld_block_data_dir, ld_block, "imputed_studies.tsv")
    standard_studies_file <- file.path(ld_block_data_dir, ld_block, "standardised_studies.tsv")
    if (!file.exists(imputed_studies_file)) return(NULL)
    
    imputed_studies <- data.table::fread(imputed_studies_file, showProgress = FALSE)
    imputed_studies <- imputed_studies[study %in% studies_db$studies$data$study_name]
    standardised_studies <- data.table::fread(standard_studies_file, showProgress = FALSE)
    standardised_studies <- standardised_studies[study %in% studies_db$studies$data$study_name & variant_type != variant_types$common]

    all_studies <- rbind(imputed_studies[, .(study, file)], standardised_studies[, .(study, file)])

    associations <- apply(all_studies, 1, function(study) {
      tryCatch({
        extractions <- data.table::fread(study[['file']], showProgress = FALSE, nThread = 1)
        needed_cols <- c("SNP", "BETA", "SE", "P", "EAF", "IMPUTED")
        missing_cols <- setdiff(needed_cols, names(extractions))
        for (col in missing_cols) {
          if (col == "IMPUTED") {
            extractions[, (col) := FALSE]
          } else {
            extractions[, (col) := NA]
          }
        }
        extractions <- extractions[, ..needed_cols]
        general_extractions <- extractions[SNP %in% general_snps$candidate_snp]
        general_extractions[, study_name := study[['study']]]
        data.table::setnames(general_extractions, tolower(names(general_extractions)))

        specific_extractions <- data.table::data.table()
        if (study[['study']] %in% specific_snps$study_name) {
          message('Processing study specific SNPs for ', study[['study']])
          specific_extractions <- specific_snps[study_name == study[['study']]]
          specific_extractions <- extractions[SNP %in% specific_extractions$candidate_snp]
          specific_extractions[, study_name := study[['study']]]
          data.table::setnames(specific_extractions, tolower(names(specific_extractions)))
        }
        all_extractions <- rbind(general_extractions, specific_extractions)
        return(all_extractions)
      }, error = function(e) {
        message('Error processing file: ', study[['file']], ' - ', e)
        return(NULL)
      })
    })
    associations <- associations[!sapply(associations, is.null)]
    if (length(associations) == 0) return(NULL)
    associations <- data.table::rbindlist(associations, fill = TRUE)

    associations <- snp_annotations_subset[associations, on = "snp"]
    associations <- studies_subset[associations, on = c("study_name" = "study_name")]
    associations <- associations[, .(snp_id, study_id, BETA=beta, SE=se, imputed, p, EAF=eaf, flipped)]

    message('Extracted ', nrow(associations), ' associations for ', ld_block)
    return(associations)
  }, error = function(e) {
    message('Error processing ld_block: ', ld_block, ' - ', e)
    return(NULL)
  })
}


split_large_dataframe_into_chunks <- function(dataframe, specific_id) {
  chunk_size <- 5000000
  num_chunks <- ceiling(nrow(dataframe) / chunk_size)
  dataframe_size <- nrow(dataframe)
  start_id <- paste0("start_", specific_id)
  stop_id <- paste0("end_", specific_id)

  dataframe_metadata_per_chunk <- lapply(seq_len(num_chunks), function(i) {
    start_boundary <- (i - 1) * chunk_size + 1
    stop_boundary <- min(i * chunk_size, dataframe_size)
    id_at_start_boundary <- dataframe[[specific_id]][start_boundary]
    id_at_stop_boundary <- dataframe[[specific_id]][stop_boundary]

    dataframe_chunk <- dataframe[start_boundary:stop_boundary, ]
    result_list <- list(
      dataframe_chunk=dataframe_chunk,
      chunk_number=i
    )
    result_list[[start_id]] <- id_at_start_boundary
    result_list[[stop_id]] <- id_at_stop_boundary
    return(result_list)
  })
  dataframe_chunks <- lapply(dataframe_metadata_per_chunk, function(x) x$dataframe_chunk)

  dataframe_metadata <- data.frame(
    sapply(dataframe_metadata_per_chunk, function(x) x$chunk_number),
    sapply(dataframe_metadata_per_chunk, function(x) x[[start_id]]),
    sapply(dataframe_metadata_per_chunk, function(x) x[[stop_id]])
  )
  names(dataframe_metadata) <- c("chunk_number", start_id, stop_id)

  return(list(dataframe_chunks=dataframe_chunks, dataframe_metadata=dataframe_metadata))
}


get_table_column_names <- function(table) {
  message('Populated ', table$name)
  lines <- unlist(strsplit(table$query, "\n")) |> trimws()
  table_column_names <- sub("^([A-Za-z0-9_]+).*", "\\1", lines)
  table_column_names <- table_column_names[
    table_column_names != "CREATE" &
    table_column_names != "FOREIGN" &
    table_column_names != "PRIMARY" &
    table_column_names != "REFERENCES" &
    table_column_names != "" &
    table_column_names != "(" &
    table_column_names != ")"
  ]

  return(table_column_names)
}

main()
