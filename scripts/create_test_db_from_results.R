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

parser <- argparser::arg_parser('Create test DuckDB from pipeline results')
parser <- argparser::add_argument(parser, '--db_file', help = 'Path to the DuckDB file to create a small test db from', type = 'character')

args <- argparser::parse_args(parser)

main <- function() {
    gpm_db_file <- args$db_file
    small_gpm_db_file <- file.path(results_dir, 'gpm_small.db')
    upload_db_file <- file.path(results_dir, 'gwas_upload_small.db')
    unlink(small_gpm_db_file)
    unlink(upload_db_file)

    gpm_con <- duckdb::dbConnect(duckdb::duckdb(), gpm_db_file)
    small_gpm_con <- duckdb::dbConnect(duckdb::duckdb(), small_gpm_db_file)
    upload_con <- duckdb::dbConnect(duckdb::duckdb(), upload_db_file)

    lapply(simple_db_tables, \(table) DBI::dbExecute(small_gpm_con, table$query))
    DBI::dbExecute(small_gpm_con, associations_table$query)
    DBI::dbExecute(upload_con, gwas_upload_table$query)

    study_sources <- DBI::dbGetQuery(gpm_con, "SELECT * FROM study_sources")
    DBI::dbAppendTable(small_gpm_con, "study_sources", study_sources)

    ld_blocks <- DBI::dbGetQuery(gpm_con, "SELECT * FROM ld_blocks")
    DBI::dbAppendTable(small_gpm_con, "ld_blocks", ld_blocks)


    coloc_group_counts <- DBI::dbGetQuery(gpm_con, "SELECT coloc_group_id, COUNT(coloc_group_id) AS n FROM colocalisations GROUP BY coloc_group_id")

    small_coloc_groups <- coloc_group_counts |>
        dplyr::filter(n > 5 & n < 10) |>
        head(100)

    small_coloc_groups <- DBI::dbGetQuery(
        gpm_con,
        sprintf("SELECT c.coloc_group_id, se.id, se.study_id, se.known_gene FROM colocalisations c
                INNER JOIN study_extractions se ON c.study_extraction_id = se.id 
                WHERE c.coloc_group_id IN (%s) AND se.known_gene IS NOT NULL",
        paste(small_coloc_groups$coloc_group_id, collapse=",")
    ))

    coloc_group_ids_to_keep <- small_coloc_groups |>
        dplyr::pull(coloc_group_id) |>
        unique()
    
    colocs_to_keep <- DBI::dbGetQuery(
        gpm_con,
        sprintf("SELECT * FROM colocalisations WHERE coloc_group_id IN (%s)",
        paste(coloc_group_ids_to_keep, collapse=",")
    ))

    study_extractions_to_keep <- DBI::dbGetQuery(
        gpm_con,
        sprintf("SELECT * FROM study_extractions WHERE id IN (%s)",
        paste(colocs_to_keep$study_extraction_id, collapse=",")
    ))
    
    studies_to_keep <- DBI::dbGetQuery(
        gpm_con,
        sprintf("SELECT * FROM studies WHERE id IN (%s)",
        paste(study_extractions_to_keep$study_id, collapse=",")
    ))

    snp_ids_to_keep <- unique(c(study_extractions_to_keep$snp_id, colocs_to_keep$snp_id))

    snp_annotations_to_keep <- DBI::dbGetQuery(
        gpm_con,
        sprintf("SELECT * FROM snp_annotations WHERE id IN (%s)",
        paste(snp_ids_to_keep, collapse=",")
    ))

    associations_to_keep <- DBI::dbGetQuery(
        gpm_con,
        sprintf("SELECT * FROM assocs WHERE snp_id IN (%s) AND study_id IN (%s)",
        paste(snp_ids_to_keep, collapse=","),
        paste(study_extractions_to_keep$study_id, collapse=",")
    ))

    ld_to_keep <- DBI::dbGetQuery(
        gpm_con,
        sprintf("SELECT * FROM ld WHERE lead_snp_id IN (%s) OR variant_snp_id IN (%s)",
        paste(snp_ids_to_keep, collapse=","),
        paste(snp_ids_to_keep, collapse=",")
    ))

    if (nrow(ld_to_keep) > 0) {
        additional_snp_annotations_to_keep <- DBI::dbGetQuery(
            gpm_con,
            sprintf("SELECT * FROM snp_annotations WHERE id IN (%s)",
            paste(ld_to_keep$variant_snp_id, collapse=",")
        ))

        snp_annotations_to_keep <- rbind(snp_annotations_to_keep, additional_snp_annotations_to_keep)
        snp_annotations_to_keep <- snp_annotations_to_keep[!duplicated(snp_annotations_to_keep$id), ]
    }

    DBI::dbAppendTable(small_gpm_con, "studies", studies_to_keep)
    DBI::dbAppendTable(small_gpm_con, "snp_annotations", snp_annotations_to_keep)
    DBI::dbAppendTable(small_gpm_con, "study_extractions", study_extractions_to_keep)
    DBI::dbAppendTable(small_gpm_con, "colocalisations", colocs_to_keep)
    DBI::dbAppendTable(small_gpm_con, "assocs", associations_to_keep)
    DBI::dbAppendTable(small_gpm_con, "ld", ld_to_keep)

    DBI::dbDisconnect(gpm_con, shutdown=TRUE)
    DBI::dbDisconnect(small_gpm_con, shutdown=TRUE)
}

main() 