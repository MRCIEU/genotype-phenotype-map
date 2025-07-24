source('constants.R')
source('database_definitions.R')

parser <- argparser::arg_parser('Create DuckDB from pipeline results')
parser <- argparser::add_argument(parser, '--studies_db_file', help = 'Studies DB file', type = 'character')
parser <- argparser::add_argument(parser, '--associations_db_file', help = 'Associations DB file', type = 'character')
parser <- argparser::add_argument(parser, '--coloc_pairs_db_file', help = 'Coloc pairs DB file', type = 'character')
parser <- argparser::add_argument(parser, '--ld_db_file', help = 'LD DB file', type = 'character')
parser <- argparser::add_argument(parser, '--gwas_upload_db_file', help = 'GWAS upload DB file', type = 'character')
parser <- argparser::add_argument(parser, '--completed_output_file', help = 'Completed output file', type = 'character')

args <- argparser::parse_args(parser)

max_cores <- 40

main <- function() {
  if (file.exists(args$studies_db_file)) file.remove(args$studies_db_file)
  if (file.exists(args$associations_db_file)) file.remove(args$associations_db_file)
  if (file.exists(args$coloc_pairs_db_file)) file.remove(args$coloc_pairs_db_file)
  if (file.exists(args$ld_db_file)) file.remove(args$ld_db_file)
  if (file.exists(args$gwas_upload_db_file)) file.remove(args$gwas_upload_db_file)

  latest_studies_conn <- duckdb::dbConnect(duckdb::duckdb(), glue::glue("{latest_results_dir}/studies.db"))
  studies_conn <- duckdb::dbConnect(duckdb::duckdb(), args$studies_db_file)
  associations_conn <- duckdb::dbConnect(duckdb::duckdb(), args$associations_db_file)
  coloc_pairs_conn <- duckdb::dbConnect(duckdb::duckdb(), args$coloc_pairs_db_file)
  ld_conn <- duckdb::dbConnect(duckdb::duckdb(), args$ld_db_file)
  gwas_upload_conn <- duckdb::dbConnect(duckdb::duckdb(), args$gwas_upload_db_file)

  lapply(studies_db, \(table) DBI::dbExecute(studies_conn, table$query))
  lapply(gwas_upload_db, \(table) DBI::dbExecute(gwas_upload_conn, table$query))
  DBI::dbExecute(associations_conn, associations_table$query)
  DBI::dbExecute(coloc_pairs_conn, coloc_pairs_table$query)
  DBI::dbExecute(ld_conn, ld_table$query)

  message('Loading data for studies db')
  studies_db <- get_existing_row_ids(latest_studies_conn, studies_db)
  studies_db <- load_data_for_studies_db(studies_db)

  message('Writing data to studies db')
  lapply(studies_db, \(table) DBI::dbAppendTable(studies_conn, table$name, table$data))

  message('Creating coloc pairs db...')
  load_data_into_coloc_pairs_db(coloc_pairs_conn, studies_db)
  q()

  all_relevant_snps <- find_relevant_snps(studies_db)
  message("Found ", nrow(all_relevant_snps), " relevant SNPs")
  message('Creating ld db...')
  load_data_into_ld_db(ld_conn, studies_db, all_relevant_snps)

  message('Creating associations db...')
  load_data_into_associations_db(associations_conn, studies_db, all_relevant_snps)

  DBI::dbDisconnect(latest_studies_conn, shutdown=TRUE)
  DBI::dbDisconnect(studies_conn, shutdown=TRUE)
  DBI::dbDisconnect(associations_conn, shutdown=TRUE)
  DBI::dbDisconnect(coloc_pairs_conn, shutdown=TRUE)
  DBI::dbDisconnect(ld_conn, shutdown=TRUE)
  DBI::dbDisconnect(gwas_upload_conn, shutdown=TRUE)

  file.copy(args$studies_db_file, file.path(latest_results_dir, "studies.db"), overwrite = TRUE)
  file.copy(args$associations_db_file, file.path(latest_results_dir, "associations.db"), overwrite = TRUE)
  file.copy(args$ld_db_file, file.path(latest_results_dir, "ld.db"), overwrite = TRUE)
  file.copy(args$gwas_upload_db_file, file.path(latest_results_dir, "gwas_upload.db"), overwrite = TRUE)

  vroom::vroom_write(data.frame(), args$completed_output_file)
}

get_existing_row_ids <- function(conn, tables) {
  for (table_name in names(tables)) {
    table <- tables[[table_name]]
    if (!is.na(table$persist_id_from)) {
      existing_ids <- DBI::dbGetQuery(conn, glue::glue("SELECT id, {table$persist_id_from} FROM {table_name}"))
      tables[[table_name]]$existing_ids <- existing_ids
    }
  }
  return(tables)
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
  return(table)
}

load_data_for_studies_db <- function(studies_db) {
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
    #TODO: remove tolower
    dplyr::rename_with(tolower) |>
    resolve_ids_for_table(studies_db$snp_annotations$existing_ids, studies_db$snp_annotations$persist_id_from) |>
    dplyr::left_join(gene_subset, by=c("gene"="gene")) |>
    dplyr::mutate(snp=trimws(snp)) |>
    dplyr::select(get_table_column_names(studies_db$snp_annotations))

  studies_db$study_extractions$data <- studies_db$study_extractions$data |>
    format_study_extractions(studies_db)

  studies_db$coloc_groups$data <- vroom::vroom(file.path(current_results_dir, "coloc_clustered_results.tsv.gz"), show_col_types = F) |>
    format_clustered_colocs(studies_db)

  studies_db$rare_results$data <- vroom::vroom(file.path(current_results_dir, "rare_results.tsv.gz"), show_col_types = F) |>
    dplyr::rename_with(tolower) |>
    format_rare_results(studies_db)

  # studies_db$results_metadata$data <- vroom::vroom(file.path(current_results_dir, "results_metadata.tsv"), show_col_types = F) |>
  #   dplyr::left_join(dplyr::select(studies_db$ld_blocks$data, ld_block, id) |>
  #   dplyr::rename(ld_block_id=id), by="ld_block")

  return(studies_db)
}

format_study_extractions <- function(study_extractions, studies_db) {
  gene_subset <- studies_db$gene_annotations$data |>
    dplyr::select(gene, id) |>
    dplyr::rename(gene_id=id)

  snp_annotations_subset <- studies_db$snp_annotations$data |>
    dplyr::select(snp, id) |>
    dplyr::rename(snp_id=id)
  
  studies_subset <- studies_db$studies$data |>
    dplyr::select(study_name, id) |>
    dplyr::rename(study_id=id)

  ld_blocks_subset <- studies_db$ld_blocks$data |>
    dplyr::select(ld_block, id) |>
    dplyr::rename(ld_block_id=id)

  study_extractions <- study_extractions |>
    resolve_ids_for_table(studies_db$study_extractions$existing_ids, studies_db$study_extractions$persist_id_from) |>
    dplyr::filter(ignore == F) |>
    dplyr::mutate(
      file=sub(ld_block_data_dir, "", file),
      svg_file=sub(data_dir, "", svg_file),
      file_with_lbfs=sub(data_dir, "", file_with_lbfs)
    ) |>
    dplyr::rename(gene=known_gene) |>
    dplyr::left_join(studies_subset, by=c("study"="study_name")) |>
    dplyr::left_join(ld_blocks_subset, by="ld_block") |>
    dplyr::left_join(gene_subset, by="gene") |>
    dplyr::left_join(snp_annotations_subset, by="snp") |>
    dplyr::select(get_table_column_names(studies_db$study_extractions))

  duplicates <- study_extractions[duplicated(study_extractions$unique_study_id), ]
  
  if (nrow(duplicates) > 0) {
    message("WARNING: Found ", nrow(duplicates), " rows with duplicate unique_study_id values.  Saving to current results dir")
    vroom::vroom_write(duplicates, file.path(current_results_dir, "duplicate_study_extractions.tsv"))
    study_extractions <- study_extractions[!duplicated(study_extractions$unique_study_id), ]
  }

  missing_snps <- study_extractions |>
    dplyr::filter(is.na(snp_id) | is.null(snp_id))
  
  if (nrow(missing_snps) > 0) {
    message("WARNING: Found ", nrow(missing_snps), " rows with missing snp_id values.  Saving to current results dir")
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
    dplyr::left_join(study_extractions_subset, by=c("unique_study_id"="unique_study_id")) |>
    dplyr::left_join(snp_annotations_subset, by="snp") |>
    dplyr::select(get_table_column_names(studies_db$coloc_groups))
  
  missing_stuff <- clustered_colocs |>
    dplyr::filter(is.na(study_extraction_id) | is.na(snp_id) | is.na(ld_block_id))
  
  if (nrow(missing_stuff) > 0) {
    message("WARNING: Found ", nrow(missing_stuff), " rows with missing study_extraction_id, snp_id, or ld_block_id values.  Saving to current results dir")
    vroom::vroom_write(missing_stuff, file.path(current_results_dir, "missing_stuff_in_clustered_colocs.tsv"))
    clustered_colocs <- clustered_colocs[!is.na(clustered_colocs$study_extraction_id) & !is.na(clustered_colocs$snp_id) & !is.na(clustered_colocs$ld_block_id), ]
  }
  
  return(clustered_colocs)
}

load_data_into_coloc_pairs_db <- function(coloc_pairs_conn, studies_db) {
  pairwise_colocs <- vroom::vroom(file.path(current_results_dir, "coloc_pairwise_results.tsv.gz"), show_col_types = F)

  study_extractions_subset <- studies_db$study_extractions$data |>
    dplyr::select(id, unique_study_id, ld_block_id) |>
    dplyr::rename(study_extraction_id=id)

  pairwise_colocs <- pairwise_colocs |>
    dplyr::filter(!is.na(h4) & ignore == F) |>
    dplyr::left_join(study_extractions_subset, by=c("unique_study_a"="unique_study_id"), relationship="many-to-one") |>
    dplyr::left_join(study_extractions_subset, by=c("unique_study_b"="unique_study_id"), relationship="many-to-one") |>
    dplyr::rename(h3=PP.H3.abf, study_extraction_a_id=study_extraction_id.x, study_extraction_b_id=study_extraction_id.y, ld_block_id=ld_block_id.x)

  missing_study_extractions <- dplyr::filter(pairwise_colocs, is.na(study_extraction_a_id) | is.na(study_extraction_b_id))
  if (nrow(missing_study_extractions) > 0) {
    message("WARNING: Found ", nrow(missing_study_extractions), " rows with missing study extractions.  Saving to current results dir")
    vroom::vroom_write(missing_study_extractions, file.path(current_results_dir, "missing_study_extractions_in_pairwise_colocs.tsv"))
    pairwise_colocs <- pairwise_colocs[!is.na(pairwise_colocs$study_extraction_a_id) & !is.na(pairwise_colocs$study_extraction_b_id), ]
  }

  pairwise_colocs <- pairwise_colocs |>
    dplyr::select(get_table_column_names(coloc_pairs_table))

  DBI::dbAppendTable(coloc_pairs_conn, coloc_pairs_table$name, pairwise_colocs)
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
    dplyr::mutate(candidate_snp=trimws(candidate_snp)) |>
    tidyr::separate_rows(traits, genes, sep=", ") |>
    dplyr::rename(unique_study_id=traits, gene=genes) |>
    dplyr::left_join(study_extractions_subset, by="unique_study_id") |>
    dplyr::left_join(snp_annotations_subset, by=c("candidate_snp"="snp"), relationship="many-to-many") |>
    dplyr::left_join(gene_annotations_subset, by="gene", relationship="many-to-one") |>
    dplyr::select(get_table_column_names(studies_db$rare_results))

  missing_study_extractions <- dplyr::filter(rare_results, is.na(study_extraction_id) | is.na(snp_id))
  if (nrow(missing_study_extractions) > 0) {
    message("WARNING: Found ", nrow(missing_study_extractions), " rows with missing study extractions.  Saving to current results dir")
    vroom::vroom_write(missing_study_extractions, file.path(current_results_dir, "missing_study_extractions_in_rare_results.tsv"))
    rare_results <- rare_results[!is.na(rare_results$study_extraction_id) & !is.na(rare_results$snp_id), ]
  }

  return(rare_results)
}

# Find relevant snps: all coloc SNPs and all finemapped SNPs that colocalise with nothing (at genome wide significance)
find_relevant_snps <- function(studies_db) {
  colocalising_snps <- dplyr::select(studies_db$colocalisations$data, candidate_snp, snp_id, ld_block) |> dplyr::distinct()

  # non_colocalising_snps <- studies_db$study_extractions$data |>
  #   dplyr::filter(!unique_study_id %in% studies_db$colocalisations$data$unique_study_id & min_p < genome_wide_p_value_threshold) |>
  #   dplyr::select(snp, snp_id, ld_block) |>
  #   dplyr::rename(candidate_snp=snp) |>
  #   dplyr::distinct()
  
  # all_relevant_snps <- dplyr::bind_rows(colocalising_snps, non_colocalising_snps) |> dplyr::distinct()
  return(colocalising_snps)
}

load_data_into_ld_db <- function(ld_conn, studies_db, all_relevant_snps) {
  relevant_ld_blocks <- unique(studies_db$study_extractions$data$ld_block)

  ld_data <- parallel::mclapply(relevant_ld_blocks, mc.cores=max_cores, \(ld_block) {
    tryCatch({
      relevant_snps <- all_relevant_snps |> dplyr::filter(ld_block == ld_block)
      generate_ld_obj(ld_block, relevant_snps$candidate_snp)
    }, error = function(e) {
      message("Error generating LD object for ", ld_block, " - ", e)
      throw(e)
    })
  }) |> dplyr::bind_rows()

  variant_annotations_subset_lead <- dplyr::select(studies_db$snp_annotations$data, snp, id) |> dplyr::rename(lead_snp_id=id)
  variant_annotations_subset_variant <- dplyr::select(studies_db$snp_annotations$data, snp, id) |> dplyr::rename(variant_snp_id=id)
  ld_blocks_subset <- dplyr::select(studies_db$ld_blocks$data, ld_block, id) |> dplyr::rename(ld_block_id=id)

  ld_data <- ld_data |>
    dplyr::left_join(variant_annotations_subset_lead, by=c("lead"="snp")) |>
    dplyr::left_join(variant_annotations_subset_variant, by=c("variant"="snp")) |>
    dplyr::left_join(ld_blocks_subset, by="ld_block") |>
    dplyr::select(-lead, -variant, -ld_block)

  DBI::dbAppendTable(ld_conn, ld_table$name, ld_data)
}

generate_ld_obj <- function(ld_block, snps) {
  message("Generating LD object for ", ld_block)
  ld_file <- file.path(ld_reference_panel_dir, glue::glue("{ld_block}.unphased.vcor1"))
  ld <- vroom::vroom(ld_file, col_names = FALSE, show_col_types = F)
  ldvars <- vroom::vroom(glue::glue("{ld_file}.vars"), delim = " ", col_names = FALSE, show_col_types = F)

  names(ld) <- ldvars$X1
  ld$lead <- ldvars$X1
  ind <- which(ldvars$X1 %in% snps)
  ld <- ld[ind,]

  ldl <- tidyr::pivot_longer(ld, cols=-lead, names_to="variant", values_to="r") |>
    dplyr::filter(r^2 > 0.8 | variant %in% ld$lead) |> 
    dplyr::filter(lead != variant) |>
    dplyr::mutate(ld_block=ld_block)
  return(ldl)
}

load_data_into_associations_db <- function(conn, studies_db, all_relevant_snps) {
  message('Retrieving ', nrow(all_relevant_snps), ' SNPs for ', nrow(studies_db$studies$data), ' studies')

  snp_annotations_subset <- studies_db$snp_annotations$data |>
    dplyr::select(snp, id) |>
    dplyr::rename(snp_id=id)

  studies_subset <- studies_db$studies$data |>
    dplyr::select(study_name, id) |>
    dplyr::rename(study_id=id)

  relevant_snps_per_ld_block <- split(all_relevant_snps, all_relevant_snps$ld_block)
  associations <- parallel::mclapply(names(relevant_snps_per_ld_block), mc.cores=max_cores, \(ld_block) {
    gc()
    message("Processing ld_block: ", ld_block)
    tryCatch({
      relevant_snps <- relevant_snps_per_ld_block[[ld_block]]
      imputed_studies_file <- file.path(ld_block_data_dir, ld_block, "imputed_studies.tsv")
      if (!file.exists(imputed_studies_file)) return(NULL)
      
      imputed_studies <- vroom::vroom(imputed_studies_file, show_col_types = F) |>
        dplyr::filter(study %in% studies_db$studies$data$study_name)
      if (nrow(imputed_studies) == 0) return(NULL)

      associations <- apply(imputed_studies, 1, \(study) {
        tryCatch({
          extractions <- vroom::vroom(study[['file']],
            show_col_types = F,
            altrep = FALSE,
            col_select = dplyr::any_of(c("SNP", "BETA", "SE", "P", "EAF", "IMPUTED"))
          )
          if (!"IMPUTED" %in% names(extractions)) {
            extractions$IMPUTED <- FALSE
          }
          extractions <- extractions |>
            dplyr::filter(SNP %in% relevant_snps$candidate_snp) |>
            dplyr::mutate(study = study[['study']]) |>
            dplyr::rename_with(tolower)
          return(extractions)
        }, error = function(e) {
          message('Error processing file: ', study[['file']], ' - ', e)
          return(NULL)
        })
      })
      associations <- associations[!sapply(associations, is.null)]
      if (length(associations) == 0) return(NULL)
      associations <- do.call(rbind, associations)

      associations <- associations |> 
        dplyr::left_join(snp_annotations_subset, by="snp") |>
        dplyr::left_join(studies_subset, by=c("study"="study_name")) |>
        dplyr::select(snp_id, study_id, beta, se, imputed, p, eaf)

      message('Extracted ', nrow(associations), ' associations for ', ld_block)

      return(associations)
    }, error = function(e) {
      message('Error processing ld_block: ', ld_block, ' - ', e)
      return(NULL)
    })
  }) 

  gc()
  associations <- associations[!sapply(associations, is.null)]
  associations <- do.call(rbind, associations)

  original_num_rows <- nrow(associations)
  associations <- associations |>
    dplyr::filter(!is.na(beta) & !is.na(se) & !is.na(p) & !is.na(eaf))
  message('Removed ', original_num_rows - nrow(associations), ' rows with missing values')

  DBI::dbAppendTable(conn, associations_table$name, associations)
  num_rows <- DBI::dbGetQuery(conn, glue::glue("SELECT COUNT(*) FROM {associations_table$name}"))
  message('Added ', num_rows$count_star, ' rows to ', associations_table$name)
  # if (num_rows$count_star > 1000000) {
    # association_files <- Sys.glob(glue::glue("{args$results_dir}/*_associations.tsv.gz"))
    # file.remove(association_files)
  # }
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
