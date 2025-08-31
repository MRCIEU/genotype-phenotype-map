source('../pipeline_steps/constants.R')
source('../pipeline_steps/database_definitions.R')

library(dplyr)
library(duckdb)
library(data.table)
library(validate)
library(tidyr)
library(R.utils)
library(purrr)
library(furrr)
library(parallel)

#NOTE: currently run with Rscript create_test_db_from_results.R --study_ids 5020

parser <- argparser::arg_parser('Create test DuckDB from pipeline results')
parser <- argparser::add_argument(parser,
                       "--study_ids",
                       help = "Space separated list of study ids",
                       type = "character",
                       nargs = Inf
)
args <- argparser::parse_args(parser)

main <- function() {
    study_ids <- split_string_into_vector(args$study_ids)

    orig_studies_db_file <- file.path(latest_results_dir, 'studies.db')
    orig_ld_db_file <- file.path(latest_results_dir, 'ld.db')
    orig_associations_db_file <- file.path(latest_results_dir, 'associations_hunmil.db')
    orig_coloc_pairs_db_file <- file.path(latest_results_dir, 'coloc_pairs.db')

    studies_db_file <- file.path(latest_results_dir, 'studies_small.db')
    ld_db_file <- file.path(latest_results_dir, 'ld_small.db')
    associations_db_file <- file.path(latest_results_dir, 'associations_small.db')
    coloc_pairs_db_file <- file.path(latest_results_dir, 'coloc_pairs_small.db')

    unlink(studies_db_file)
    unlink(ld_db_file)
    unlink(associations_db_file)
    unlink(coloc_pairs_db_file)

    orig_studies_con <- duckdb::dbConnect(duckdb::duckdb(), orig_studies_db_file, read_only=TRUE)
    orig_ld_con <- duckdb::dbConnect(duckdb::duckdb(), orig_ld_db_file, read_only=TRUE)
    orig_associations_con <- duckdb::dbConnect(duckdb::duckdb(), orig_associations_db_file, read_only=TRUE)
    orig_coloc_pairs_con <- duckdb::dbConnect(duckdb::duckdb(), orig_coloc_pairs_db_file, read_only=TRUE)

    studies_con <- duckdb::dbConnect(duckdb::duckdb(), studies_db_file)
    ld_con <- duckdb::dbConnect(duckdb::duckdb(), ld_db_file)
    associations_con <- duckdb::dbConnect(duckdb::duckdb(), associations_db_file)
    coloc_pairs_con <- duckdb::dbConnect(duckdb::duckdb(), coloc_pairs_db_file)

    #remove text of foreign key constraints from all tables
    lapply(studies_db, \(table) {
        query <- sub(",\n.\\s+FOREIGN.*", ")", table$query)
        DBI::dbExecute(studies_con, query)
    })
    DBI::dbExecute(ld_con, ld_table$query)
    DBI::dbExecute(coloc_pairs_con, coloc_pairs_significant_table$query)

    study_sources <- DBI::dbGetQuery(orig_studies_con, "SELECT * FROM study_sources")
    DBI::dbAppendTable(studies_con, "study_sources", study_sources)

    ld_blocks <- DBI::dbGetQuery(orig_studies_con, "SELECT * FROM ld_blocks")
    DBI::dbAppendTable(studies_con, "ld_blocks", ld_blocks)

    # results_metadata <- DBI::dbGetQuery(orig_studies_con, "SELECT * FROM results_metadata")
    # DBI::dbAppendTable(studies_con, "results_metadata", results_metadata)

    coloc_groups <- DBI::dbGetQuery(orig_studies_con,
        sprintf("SELECT coloc_group_id FROM coloc_groups WHERE study_id IN (%s)",
        paste(study_ids, collapse=",")
    ))
    coloc_groups <- DBI::dbGetQuery(orig_studies_con,
        sprintf("SELECT * FROM coloc_groups WHERE coloc_group_id IN (%s)",
        paste(coloc_groups$coloc_group_id, collapse=",")
    ))

    study_ids <- unique(c(coloc_groups$study_id, study_ids))
    DBI::dbAppendTable(studies_con, "coloc_groups", coloc_groups)

    studies <- DBI::dbGetQuery(orig_studies_con,
        sprintf("SELECT * FROM studies WHERE id IN (%s)", paste(study_ids, collapse=",")
    ))
    DBI::dbAppendTable(studies_con, "studies", studies)

    traits <- DBI::dbGetQuery(orig_studies_con,
        sprintf("SELECT * FROM traits WHERE id IN (%s)", paste(studies$trait_id, collapse=",")))
    DBI::dbAppendTable(studies_con, "traits", traits)

    study_extractions <- DBI::dbGetQuery(orig_studies_con,
        sprintf("SELECT * FROM study_extractions WHERE id IN (%s)",
        paste(coloc_groups$study_extraction_id, collapse=",")
    ))

    all_study_ids <- unique(c(traits$study_id, studies$id))
    DBI::dbAppendTable(studies_con, "study_extractions", study_extractions)
    rare_results <- DBI::dbGetQuery(orig_studies_con,
        sprintf("SELECT * FROM rare_results WHERE study_id IN (%s)",
        paste(all_study_ids, collapse=",")
    ))
    DBI::dbAppendTable(studies_con, "rare_results", rare_results)

    snp_ids_to_keep <- unique(c(study_extractions$snp_id, coloc_groups$snp_id, rare_results$snp_id))

    snp_annotations <- DBI::dbGetQuery(orig_studies_con,
        sprintf("SELECT * FROM snp_annotations WHERE id IN (%s)",
        paste(snp_ids_to_keep, collapse=",")
    ))
    DBI::dbAppendTable(studies_con, "snp_annotations", snp_annotations)

    genes_to_keep <- unique(na.omit(c(study_extractions$gene_id, rare_results$gene_id, studies$gene_id)))

    gene_annotations <- DBI::dbGetQuery(orig_studies_con,
        sprintf("SELECT * FROM gene_annotations WHERE id IN (%s)",
        paste(genes_to_keep, collapse=",")
    ))
    DBI::dbAppendTable(studies_con, "gene_annotations", gene_annotations)


    coloc_pairs <- DBI::dbGetQuery(orig_coloc_pairs_con,
        sprintf("SELECT * FROM coloc_pairs WHERE study_extraction_a_id IN (%s) AND study_extraction_b_id IN (%s)",
        paste(coloc_groups$study_extraction_id, collapse=","),
        paste(coloc_groups$study_extraction_id, collapse=",")
    ))
    DBI::dbAppendTable(coloc_pairs_con, "coloc_pairs", coloc_pairs)

    associations_metadata <- DBI::dbGetQuery(orig_associations_con, "SELECT * FROM associations_metadata")
    DBI::dbExecute(associations_con, associations_db$associations_metadata$query)
    DBI::dbAppendTable(associations_con, "associations_metadata", associations_metadata)

    lapply(associations_metadata$associations_table_name, function(table_name) {
        create_table_query <- sub("table_name", table_name, associations_db$associations$query)
        DBI::dbExecute(associations_con, create_table_query)
        associations <- DBI::dbGetQuery(orig_associations_con,
          sprintf("SELECT * FROM %s WHERE snp_id IN (%s) AND study_id IN (%s)",
          table_name,
          paste(snp_ids_to_keep, collapse=","),
          paste(studies$id, collapse=","))
        )
        DBI::dbAppendTable(associations_con, table_name, associations)
    })
    
    ld_to_keep <- DBI::dbGetQuery(
        orig_ld_con,
        sprintf("SELECT * FROM ld WHERE lead_snp_id IN (%s) OR variant_snp_id IN (%s)",
        paste(snp_ids_to_keep, collapse=","),
        paste(snp_ids_to_keep, collapse=",")
    ))
    DBI::dbAppendTable(ld_con, "ld", ld_to_keep)

    DBI::dbDisconnect(studies_con, shutdown=TRUE)
    DBI::dbDisconnect(ld_con, shutdown=TRUE)
    DBI::dbDisconnect(associations_con, shutdown=TRUE)
    DBI::dbDisconnect(coloc_pairs_con, shutdown=TRUE)
}

split_string_into_vector <- function(input_string) {
  return(unlist(strsplit(input_string, '[ ]')))
}

main() 