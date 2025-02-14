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
parser <- argparser::add_argument(parser, '--n_colocs', help = 'Number of studies to sample', type = 'numeric', default = 10)

args <- argparser::parse_args(parser)

main <- function() {
    last_pipeline_run <- head(sort(list.files(results_dir, pattern = "^20", include.dirs=T), decreasing=T), 1)[1]
    last_result_dir <- glue::glue('{results_dir}/{last_pipeline_run}')
    
    # Create test database paths
    test_studies_db <- file.path(last_result_dir, 'test_studies.db')
    test_associations_db <- file.path(last_result_dir, 'test_associations.db')
    unlink(test_studies_db)
    unlink(test_associations_db)

    # Connect to source databases
    studies_con <- dbConnect(duckdb::duckdb(), file.path(last_result_dir, 'studies.db'), read_only=TRUE)
    assocs_con <- dbConnect(duckdb::duckdb(), file.path(last_result_dir, 'associations.db'), read_only=TRUE)
    
    # Sample random studies
    colocs_processed <- dbGetQuery(studies_con, sprintf("SELECT * FROM coloc WHERE id IN (1,2,3,4,5,6,7,8,9,10)"))
    sampled_studies <- unique(unlist(colocs_processed %>% pull(traits)))


    # Get coloc results
    coloc <- dbGetQuery(studies_con, sprintf("
        SELECT * FROM coloc 
        WHERE traits IN ('%s')",
        paste(sampled_studies, collapse="','")
    ))
    
    # Get related data for sampled studies
    all_study_blocks <- dbGetQuery(studies_con, sprintf("
        SELECT * FROM all_study_blocks 
        WHERE unique_study_id IN ('%s')", 
        paste(sampled_studies, collapse="','")
    ))

    studies_processed <- dbGetQuery(studies_con, sprintf("
        SELECT * FROM studies_processed 
        WHERE study_name IN ('%s')",
        paste(all_study_blocks$study, collapse="','")
    ))
    
    # Get unique ld_blocks for these studies
    ld_blocks <- unique(all_study_blocks$ld_block)
    
    # Get related metadata
    results_metadata <- dbGetQuery(studies_con, sprintf("
        SELECT * FROM results_metadata 
        WHERE ld_block IN ('%s')",
        paste(ld_blocks, collapse="','")
    ))
    
    # Get variant annotations for these regions
    variant_annotations <- dbGetQuery(studies_con, sprintf("
        SELECT * FROM variant_annotations 
        WHERE ld_block IN ('%s')",
        paste(ld_blocks, collapse="','")
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
        WHERE study IN ('%s')",
        paste(all_study_blocks$study, collapse="','")
    ))
    
    # Create test databases
    test_studies_con <- dbConnect(duckdb::duckdb(), test_studies_db)

    dbWriteTable(test_studies_con, "all_study_blocks", all_study_blocks)
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

ensure_test_dbs_are_valid <- function(test_studies_db, test_associations_db) {
    # Connect to test databases
    test_studies_con <- dbConnect(duckdb::duckdb(), test_studies_db, read_only=TRUE)
    test_assocs_con <- dbConnect(duckdb::duckdb(), test_associations_db, read_only=TRUE)
    
    # Get sample counts
    study_counts <- dbGetQuery(test_studies_con, "
        SELECT 
            (SELECT COUNT(*) FROM all_study_blocks) as study_blocks,
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
    message("Study blocks: ", study_counts$study_blocks)
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