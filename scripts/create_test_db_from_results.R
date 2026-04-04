source("../pipeline_steps/constants.R")
source("../pipeline_steps/database_definitions.R")

library(dplyr)
library(duckdb)
library(data.table)
library(validate)
library(tidyr)
library(R.utils)
library(purrr)
library(furrr)
library(parallel)

# NOTE: currently run with Rscript create_test_db_from_results.R --study_ids 5020 54929

parser <- argparser::arg_parser("Create test DuckDB from pipeline results")
parser <- argparser::add_argument(parser,
  "--study_ids",
  help = "Space separated list of study ids",
  type = "character",
  nargs = Inf
)
args <- argparser::parse_args(parser)

main <- function() {
  study_ids <- split_string_into_vector(args$study_ids)

  orig_studies_db_file <- file.path(latest_results_dir, "studies.db")
  orig_ld_db_file <- file.path(latest_results_dir, "ld.db")
  orig_associations_db_file <- file.path(latest_results_dir, "associations.db")
  orig_associations_full_db_file <- file.path(latest_results_dir, "associations_full.db")
  orig_coloc_pairs_db_file <- file.path(latest_results_dir, "coloc_pairs.db")

  studies_db_file <- file.path(latest_results_dir, "studies_small.db")
  ld_db_file <- file.path(latest_results_dir, "ld_small.db")
  associations_db_file <- file.path(latest_results_dir, "associations_small.db")
  associations_full_db_file <- file.path(latest_results_dir, "associations_full_small.db")
  coloc_pairs_db_file <- file.path(latest_results_dir, "coloc_pairs_small.db")

  unlink(studies_db_file)
  unlink(ld_db_file)
  unlink(associations_db_file)
  unlink(associations_full_db_file)
  unlink(coloc_pairs_db_file)

  orig_studies_con <- duckdb::dbConnect(duckdb::duckdb(), orig_studies_db_file, read_only = TRUE)
  orig_ld_con <- duckdb::dbConnect(duckdb::duckdb(), orig_ld_db_file, read_only = TRUE)
  orig_associations_con <- duckdb::dbConnect(duckdb::duckdb(), orig_associations_db_file, read_only = TRUE)
  orig_associations_full_con <- duckdb::dbConnect(duckdb::duckdb(), orig_associations_full_db_file, read_only = TRUE)
  orig_coloc_pairs_con <- duckdb::dbConnect(duckdb::duckdb(), orig_coloc_pairs_db_file, read_only = TRUE)

  studies_con <- duckdb::dbConnect(duckdb::duckdb(), studies_db_file)
  ld_con <- duckdb::dbConnect(duckdb::duckdb(), ld_db_file)
  associations_con <- duckdb::dbConnect(duckdb::duckdb(), associations_db_file)
  associations_full_con <- duckdb::dbConnect(duckdb::duckdb(), associations_full_db_file)
  coloc_pairs_con <- duckdb::dbConnect(duckdb::duckdb(), coloc_pairs_db_file)

  # remove text of foreign key constraints from all tables
  lapply(studies_db, \(table) {
    query <- sub(",\n.\\s+FOREIGN.*", ")", table$query)
    DBI::dbExecute(studies_con, query)
    return()
  })
  safe_lapply(additional_studies_tables, function(table) {
    DBI::dbExecute(studies_con, table$query)
    # DBI::dbExecute(studies_con, table$indexes)
    return()
  })

  DBI::dbExecute(ld_con, ld_table$query)
  DBI::dbExecute(coloc_pairs_con, coloc_pairs_significant_table$query)

  study_sources <- DBI::dbGetQuery(orig_studies_con, "SELECT * FROM study_sources")
  DBI::dbAppendTable(studies_con, "study_sources", study_sources)

  ld_blocks <- DBI::dbGetQuery(orig_studies_con, "SELECT * FROM ld_blocks")
  DBI::dbAppendTable(studies_con, "ld_blocks", ld_blocks)

  coloc_groups <- DBI::dbGetQuery(
    orig_studies_con,
    sprintf(
      "SELECT coloc_group_id FROM coloc_groups WHERE study_id IN (%s)",
      paste(study_ids, collapse = ",")
    )
  )
  coloc_groups <- DBI::dbGetQuery(
    orig_studies_con,
    sprintf(
      "SELECT * FROM coloc_groups WHERE coloc_group_id IN (%s)",
      paste(coloc_groups$coloc_group_id, collapse = ",")
    )
  )

  coloc_groups_wide <- DBI::dbGetQuery(
    orig_studies_con,
    sprintf(
      "SELECT * FROM coloc_groups_wide WHERE coloc_group_id IN (%s)",
      paste(coloc_groups$coloc_group_id, collapse = ",")
    )
  )

  study_ids <- unique(c(coloc_groups$study_id, study_ids))
  DBI::dbAppendTable(studies_con, "coloc_groups", coloc_groups)
  DBI::dbAppendTable(studies_con, "coloc_groups_wide", coloc_groups_wide)


  studies <- DBI::dbGetQuery(
    orig_studies_con,
    sprintf("SELECT * FROM studies WHERE id IN (%s)", paste(study_ids, collapse = ","))
  )
  DBI::dbAppendTable(studies_con, "studies", studies)

  traits <- DBI::dbGetQuery(
    orig_studies_con,
    sprintf("SELECT * FROM traits WHERE id IN (%s)", paste(studies$trait_id, collapse = ","))
  )
  DBI::dbAppendTable(studies_con, "traits", traits)

  study_extractions <- DBI::dbGetQuery(
    orig_studies_con,
    sprintf(
      "SELECT * FROM study_extractions WHERE id IN (%s)",
      paste(coloc_groups$study_extraction_id, collapse = ",")
    )
  )
  study_extractions_wide <- DBI::dbGetQuery(
    orig_studies_con,
    sprintf(
      "SELECT * FROM study_extractions_wide WHERE id IN (%s)",
      paste(coloc_groups$study_extraction_id, collapse = ",")
    )
  )

  all_study_ids <- unique(c(traits$study_id, studies$id))
  DBI::dbAppendTable(studies_con, "study_extractions", study_extractions)
  DBI::dbAppendTable(studies_con, "study_extractions_wide", study_extractions_wide)

  rare_results_groups <- DBI::dbGetQuery(
    orig_studies_con,
    sprintf(
      "SELECT rare_result_group_id FROM rare_results WHERE study_id IN (%s)",
      paste(all_study_ids, collapse = ",")
    )
  )
  rare_results_groups <- DBI::dbGetQuery(
    orig_studies_con,
    sprintf(
      "SELECT * FROM rare_results WHERE rare_result_group_id IN (%s)",
      paste(rare_results_groups$rare_result_group_id, collapse = ",")
    )
  )
  DBI::dbAppendTable(studies_con, "rare_results", rare_results_groups)

  rare_results <- DBI::dbGetQuery(
    orig_studies_con,
    sprintf(
      "SELECT * FROM rare_results_wide WHERE rare_result_group_id IN (%s)",
      paste(rare_results_groups$rare_result_group_id, collapse = ",")
    )
  )
  DBI::dbAppendTable(studies_con, "rare_results_wide", rare_results)

  variant_ids_to_keep <- unique(c(study_extractions$variant_id, coloc_groups$variant_id, rare_results$variant_id))

  variant_annotations <- DBI::dbGetQuery(
    orig_studies_con,
    sprintf(
      "SELECT * FROM variant_annotations WHERE id IN (%s)",
      paste(variant_ids_to_keep, collapse = ",")
    )
  )
  DBI::dbAppendTable(studies_con, "variant_annotations", variant_annotations)

  genes_to_keep <- unique(na.omit(c(study_extractions$gene_id, rare_results$gene_id, studies$gene_id)))

  gene_annotations <- DBI::dbGetQuery(
    orig_studies_con,
    sprintf(
      "SELECT * FROM gene_annotations WHERE id IN (%s)",
      paste(genes_to_keep, collapse = ",")
    )
  )
  DBI::dbAppendTable(studies_con, "gene_annotations", gene_annotations)

  variant_pleiotropy <- DBI::dbGetQuery(
    orig_studies_con,
    sprintf(
      "SELECT * FROM variant_pleiotropy WHERE variant_id IN (%s)",
      paste(variant_ids_to_keep, collapse = ",")
    )
  )
  DBI::dbAppendTable(studies_con, "variant_pleiotropy", variant_pleiotropy)

  gene_pleiotropy <- DBI::dbGetQuery(
    orig_studies_con,
    sprintf(
      "SELECT * FROM gene_pleiotropy WHERE gene_id IN (%s)",
      paste(genes_to_keep, collapse = ",")
    )
  )
  DBI::dbAppendTable(studies_con, "gene_pleiotropy", gene_pleiotropy)

  coloc_pairs <- DBI::dbGetQuery(
    orig_coloc_pairs_con,
    sprintf(
      "SELECT * FROM coloc_pairs WHERE study_extraction_a_id IN (%s) AND study_extraction_b_id IN (%s)",
      paste(coloc_groups$study_extraction_id, collapse = ","),
      paste(coloc_groups$study_extraction_id, collapse = ",")
    )
  )
  DBI::dbAppendTable(coloc_pairs_con, "coloc_pairs", coloc_pairs)


  create_associations_table_query <- sub("table_name", "associations", associations_db$associations$query)
  DBI::dbExecute(associations_con, create_associations_table_query)
  associations <- DBI::dbGetQuery(
    orig_associations_con,
    sprintf(
      "SELECT * FROM associations WHERE variant_id IN (%s) AND study_id IN (%s)",
      paste(variant_ids_to_keep, collapse = ","),
      paste(studies$id, collapse = ",")
    )
  )
  DBI::dbAppendTable(associations_con, "associations", associations)

  DBI::dbExecute(associations_full_con, associations_db$associations_metadata$query)
  associations_metadata <- DBI::dbGetQuery(orig_associations_full_con, "SELECT * FROM associations_metadata")
  all_associations <- lapply(associations_metadata$associations_table_name, function(table_name) {
    create_table_query <- sub("table_name", table_name, associations_db$associations$query)
    DBI::dbExecute(associations_full_con, create_table_query)
    associations <- DBI::dbGetQuery(
      orig_associations_full_con,
      sprintf(
        "SELECT * FROM %s WHERE variant_id IN (%s) AND study_id IN (%s)",
        table_name,
        paste(variant_ids_to_keep, collapse = ","),
        paste(studies$id, collapse = ",")
      )
    )
    return(associations)
  })
  all_associations <- do.call(rbind, all_associations)
  DBI::dbAppendTable(associations_full_con, "associations_1", all_associations)
  associations_metadata <- data.frame(
    start_variant_id = 1,
    stop_variant_id = .Machine$integer.max,
    associations_table_name = "associations_1"
  )
  DBI::dbAppendTable(associations_full_con, "associations_metadata", associations_metadata)

  ld_to_keep <- DBI::dbGetQuery(
    orig_ld_con,
    sprintf(
      "SELECT * FROM ld WHERE lead_variant_id IN (%s) OR proxy_variant_id IN (%s)",
      paste(variant_ids_to_keep, collapse = ","),
      paste(variant_ids_to_keep, collapse = ",")
    )
  )
  DBI::dbAppendTable(ld_con, "ld", ld_to_keep)

  DBI::dbDisconnect(studies_con, shutdown = TRUE)
  DBI::dbDisconnect(ld_con, shutdown = TRUE)
  DBI::dbDisconnect(associations_con, shutdown = TRUE)
  DBI::dbDisconnect(associations_full_con, shutdown = TRUE)
  DBI::dbDisconnect(coloc_pairs_con, shutdown = TRUE)
  return()
}

split_string_into_vector <- function(input_string) {
  return(unlist(strsplit(input_string, "[ ]")))
}

main()
