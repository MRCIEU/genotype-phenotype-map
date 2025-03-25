source('../pipeline_steps/constants.R')

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
parser <- argparser::add_argument(parser, '--small_db', help = 'Create a small test db', type = 'logical', default = FALSE)
parser <- argparser::add_argument(parser, '--n_colocs', help = 'Number of colocs to sample', type = 'numeric', default = 10)
parser <- argparser::add_argument(parser, '--genes', help = 'Genes to sample', type = 'character', default = NULL)
parser <- argparser::add_argument(parser, '--studies', help = 'Studies to sample', type = 'character', default = NULL)

args <- argparser::parse_args(parser)

main <- function() {
    last_pipeline_run <- head(sort(list.files(results_dir, pattern = "^20", include.dirs=T), decreasing=T), 1)[1]
    gpm_db_file <- file.path(results_dir, 'gpm.db')

    if (args$small_db) {
        create_small_test_db(gpm_db_file)
    }
    q()
    
    # Create test database paths
    test_studies_db <- file.path(last_result_dir, 'test_studies.db')
    test_associations_db <- file.path(last_result_dir, 'test_associations.db')
    unlink(test_studies_db)
    unlink(test_associations_db)

    # Connect to source databases
    studies_con <- dbConnect(duckdb::duckdb(), file.path(last_result_dir, 'studies.db'), read_only=TRUE)
    assocs_con <- dbConnect(duckdb::duckdb(), file.path(last_result_dir, 'associations.db'), read_only=TRUE)
    
    colocs_processed <- dbGetQuery(studies_con, sprintf("
        SELECT * FROM coloc
        WHERE id IN (1,2,3,4,5,6,7,8,9,10)"
    ))
    sampled_studies <- unique(unlist(colocs_processed %>% pull(study)))
    print(sampled_studies)
    
    if (!is.null(args$genes)) {
        additional_studies_processed <- dbGetQuery(studies_con, sprintf("
            SELECT * FROM studies_processed 
            WHERE gene IN ('%s')",
            paste(args$genes, collapse="','")
        ))
        sampled_studies <- unique(c(sampled_studies, unlist(additional_studies_processed %>% pull(study_name))))
    }
    if (!is.null(args$studies)) {
        additional_studies_processed <- dbGetQuery(studies_con, sprintf("
            SELECT * FROM studies_processed 
            WHERE study_name IN ('%s')",
            paste(args$studies, collapse="','")
        ))
        sampled_studies <- unique(c(sampled_studies, unlist(additional_studies_processed %>% pull(study_name))))
    }

    # get all coloc ids that are in these studies, THEN grab all the studies that have those coloc ids
    coloc <- dbGetQuery(studies_con, sprintf("
        SELECT * FROM coloc 
        WHERE study IN ('%s')",
        paste(sampled_studies, collapse="','")
    ))
    coloc_ids <- unique(coloc$id)

    coloc <- dbGetQuery(studies_con, sprintf("
        SELECT * FROM coloc 
        WHERE id IN (%s)",
        paste(coloc_ids, collapse=",")
    ))

    sampled_studies <- unique(coloc$traits)
    print(length(sampled_studies))

    # Get related data for sampled studies
    study_extractions <- dbGetQuery(studies_con, sprintf("
        SELECT * FROM study_extractions 
        WHERE unique_study_id IN ('%s')", 
        paste(sampled_studies, collapse="','")
    ))

    studies_processed <- dbGetQuery(studies_con, sprintf("
        SELECT * FROM studies_processed 
        WHERE study_name IN ('%s')",
        paste(study_extractions$study, collapse="','")
    ))
    
    # Get unique ld_blocks for these studies
    ld_blocks <- unique(study_extractions$ld_block)
    genes <- unique(study_extractions$known_gene)
    
    # Get related metadata
    results_metadata <- dbGetQuery(studies_con, sprintf("
        SELECT * FROM results_metadata 
        WHERE ld_block IN ('%s')",
        paste(ld_blocks, collapse="','")
    ))
    
    # Get variant annotations for these regions
    variant_annotations <- dbGetQuery(studies_con, sprintf("
        SELECT * FROM variant_annotations 
        WHERE symbol IN ('%s')",
        paste(genes, collapse="','")
    ))
    
    # Get LD information
    ld <- dbGetQuery(studies_con, sprintf("
        SELECT * FROM ld 
        WHERE ld_block IN ('%s')",
        paste(ld_blocks, collapse="','")
    ))
    
    # Get associations
    associations <- dbGetQuery(assocs_con, sprintf("
        SELECT * FROM assocs 
        WHERE study IN ('%s') AND SNP IN ('%s')",
        paste(study_extractions$study, collapse="','"),
        paste(coloc$candidate_snp, collapse="','")
    ))
    
    # Create test databases
    test_studies_con <- dbConnect(duckdb::duckdb(), test_studies_db)

    dbWriteTable(test_studies_con, "study_extractions", study_extractions)
    dbWriteTable(test_studies_con, "results_metadata", results_metadata)
    dbWriteTable(test_studies_con, "studies_processed", studies_processed)
    dbWriteTable(test_studies_con, "variant_annotations", variant_annotations)
    dbWriteTable(test_studies_con, "coloc", coloc)
    dbWriteTable(test_studies_con, "ld", ld)
    
    test_assocs_con <- dbConnect(duckdb::duckdb(), test_associations_db)
    dbWriteTable(test_assocs_con, "assocs", associations)
    
    # Cleanup
    dbDisconnect(studies_con, shutdown=TRUE)
    dbDisconnect(assocs_con, shutdown=TRUE)
    dbDisconnect(test_studies_con, shutdown=TRUE)
    dbDisconnect(test_assocs_con, shutdown=TRUE)
    
    # Validate test databases
    ensure_test_dbs_are_valid(test_studies_db, test_associations_db)
}

create_small_test_db <- function(gpm_db_file) {
    small_gpm_db_file <- file.path(results_dir, 'gpm_small.db')
    file.copy(gpm_db_file, small_gpm_db_file)

    gpm_con <- duckdb::dbConnect(duckdb::duckdb(), small_gpm_db_file)

    small_coloc_groups <- DBI::dbGetQuery(gpm_con, "select count(*) from studies")
    print(small_coloc_groups)

    coloc_group_counts <- DBI::dbGetQuery(gpm_con, "SELECT coloc_group_id, COUNT(coloc_group_id) AS n FROM colocalisations GROUP BY coloc_group_id")

    small_coloc_groups <- coloc_group_counts |>
        dplyr::filter(n > 5 & n < 10) |>
        head(50)

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
        sprintf("SELECT coloc_group_id, study_extraction_id, snp_id FROM colocalisations WHERE coloc_group_id IN (%s)",
        paste(coloc_group_ids_to_keep, collapse=",")
    ))

    study_extractions_to_keep <- DBI::dbGetQuery(
        gpm_con,
        sprintf("SELECT id, study_id, snp_id FROM study_extractions WHERE id IN (%s)",
        paste(colocs_to_keep$study_extraction_id, collapse=",")
    ))
    
    studies_to_keep <- DBI::dbGetQuery(
        gpm_con,
        sprintf("SELECT id FROM studies WHERE id IN (%s)",
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

    print(nrow(associations_to_keep))
    print(nrow(study_extractions_to_keep))
    print(nrow(colocs_to_keep))
    print(nrow(snp_annotations_to_keep))
    print(nrow(studies_to_keep))
    print(head(studies_to_keep))

    DBI::dbExecute(gpm_con, "CREATE TEMP TABLE studies_to_keep AS SELECT id FROM studies WHERE FALSE")
    DBI::dbExecute(gpm_con, sprintf("INSERT INTO studies_to_keep VALUES %s",
        paste(sprintf("(%s)", studies_to_keep$id), collapse=",")
    ))
    print('studies to keep inserted')
    DBI::dbExecute(gpm_con, "DELETE FROM studies WHERE id NOT IN (SELECT id FROM studies_to_keep)")
    DBI::dbExecute(gpm_con, "DROP TABLE studies_to_keep")
    print('studies deleted')

    DBI::dbExecute(gpm_con, sprintf("DELETE FROM colocalisations WHERE coloc_group_id NOT IN (%s)",
        paste(colocs_to_keep$coloc_group_id, collapse=",")
    ))
    print('colocs deleted')

    DBI::dbExecute(gpm_con, sprintf("DELETE FROM study_extractions WHERE id NOT IN (%s)",
        paste(study_extractions_to_keep$id, collapse=",")
    ))
    print('study extractions deleted')
    DBI::dbExecute(gpm_con, sprintf("DELETE FROM snp_annotations WHERE id NOT IN (%s)",
        paste(snp_annotations_to_keep$id, collapse=",")
    ))
    print('snp annotations deleted')

    DBI::dbExecute(gpm_con, sprintf("DELETE FROM assocs WHERE snp_id NOT IN (%s) AND study_id NOT IN (%s)",
        paste(study_extractions_to_keep$snp_id, collapse=","),
        paste(study_extractions_to_keep$study_id, collapse=",")
    ))
    print('associations deleted')

    DBI::dbExecute(gpm_con, "TRUNCATE TABLE ld, rare_results")
    print('ld and rare results truncated')
    DBI::dbDisconnect(gpm_con, shutdown=TRUE)
}

ensure_test_dbs_are_valid <- function(test_studies_db, test_associations_db) {
    # Connect to test databases
    test_studies_con <- dbConnect(duckdb::duckdb(), test_studies_db, read_only=TRUE)
    test_assocs_con <- dbConnect(duckdb::duckdb(), test_associations_db, read_only=TRUE)
    
    # Get sample counts
    study_counts <- dbGetQuery(test_studies_con, "
        SELECT 
            (SELECT COUNT(*) FROM study_extractions) as study_extractions,
            (SELECT COUNT(*) FROM results_metadata) as metadata,
            (SELECT COUNT(*) FROM studies_processed) as studies,
            (SELECT COUNT(*) FROM variant_annotations) as variants,
            (SELECT COUNT(*) FROM coloc) as coloc,
            (SELECT COUNT(*) FROM ld) as ld
    ")
    
    assoc_count <- dbGetQuery(test_assocs_con, "SELECT COUNT(*) as associations FROM assocs")
    
    # Print summary
    message("Test database summary:")
    message("Studies processed: ", study_counts$studies)
    message("Study extractions: ", study_counts$study_extractions)
    message("Results metadata: ", study_counts$metadata)
    message("Variant annotations: ", study_counts$variants)
    message("Coloc results: ", study_counts$coloc)
    message("LD relationships: ", study_counts$ld)
    message("Associations: ", assoc_count$associations)
    
    # Cleanup
    dbDisconnect(test_studies_con, shutdown=TRUE)
    dbDisconnect(test_assocs_con, shutdown=TRUE)
}

main() 