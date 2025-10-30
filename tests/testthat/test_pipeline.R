library(testthat)

setwd('../../')
source('pipeline_steps/constants.R')

test_that("Identify studies to process", {
  output <- system(glue::glue("Rscript pipeline_steps/identify_studies_to_process.R"),
    wait = TRUE,
    intern = TRUE,
    ignore.stdout = FALSE,
    ignore.stderr = FALSE
  )
  has_errors <- grepl("error", output, ignore.case = TRUE)
  expect_false(any(has_errors), info = "Identify studies to process should not contain errors")

  studies_to_process_file <- glue::glue('{pipeline_metadata_dir}/studies_to_process.tsv')
  expect_true(file.exists(studies_to_process_file), info = "Studies to process file should exist")
  studies_to_process <- vroom::vroom(studies_to_process_file, show_col_types = FALSE)
  expect_equal(nrow(studies_to_process), 11, info = "Studies to process file should have 11 studies")
})

# test_that("Pipeline execution and file validation", {
#   output <- system(glue::glue("snakemake --profile ./"),
#     wait = TRUE,
#     intern = TRUE,
#     ignore.stdout = FALSE,
#     ignore.stderr = FALSE
#   )
#   has_errors <- grepl("error", output, ignore.case = TRUE)
#   expect_false(any(has_errors), info = "Pipeline execution should not contain errors")
  
#   test_that("Main result files are created", {
#     expected_tsv_files <- c(
#       file.path(results_dir, "current/studies_processed.tsv.gz"),
#       file.path(results_dir, "current/traits_processed.tsv.gz"),
#       file.path(results_dir, "current/study_extractions.tsv.gz"),
#       file.path(results_dir, "current/coloc_clustered_results.tsv.gz"),
#       file.path(results_dir, "current/coloc_pairwise_results.tsv.gz"),
#       file.path(results_dir, "current/rare_results.tsv.gz")
#     )

#     expected_db_files <- c(
#       file.path(results_dir, "current/studies.db"),
#       file.path(results_dir, "current/associations.db"),
#       file.path(results_dir, "current/coloc_pairs_full.db"),
#       file.path(results_dir, "current/coloc_pairs.db"),
#       file.path(results_dir, "current/ld.db"),
#       file.path(results_dir, "current/gwas_upload.db")
#     )

#     for (file_path in expected_tsv_files) {
#       expect_true(file.exists(file_path), info = glue("File should exist: {file_path}"))
      
#       file_size <- file.size(file_path)
#       expect_gt(file_size, 0, info = glue("File should not be empty: {file_path}"))
#     }
#   })

#   # test_that("traits_processed and studies_processed are valid", {
#   #   traits_file <- file.path(results_dir, "current/traits_processed.tsv.gz")
#   #   expect_true(file.exists(traits_file), info = glue("File should exist: {traits_file}"))
#   #   expect_gt(file.size(traits_file), 0, info = glue("File should not be empty: {traits_file}"))

#   #   studies_file <- file.path(results_dir, "current/studies_processed.tsv.gz")
#   #   expect_true(file.exists(studies_file), info = glue("File should exist: {studies_file}"))
#   #   expect_gt(file.size(studies_file), 0, info = glue("File should not be empty: {studies_file}"))
    
#   #   studies <- vroom::vroom(studies_file, show_col_types = FALSE)
#   #   traits <- vroom::vroom(traits_file, show_col_types = FALSE)

#   #   expect_true(nrow(studies) == nrow(traits) && nrow(studies) > 0, info = "Number of studies and traits should be equal")
#   #   expect_true(all(studies$study_name %in% traits$study_name), info = "All studies should be in traits_processed")

#   #   joined_studies_traits <- dplyr::inner_join(studies, traits, by = "study_name")
#   #   expect_true(nrow(joined_studies_traits) == nrow(studies) && nrow(joined_studies_traits) > 0, info = "All studies should have an associated trait")

#   #   sapply(joined_studies_traits$extracted_location, function(extracted_location) {
#   #     expect_true(dir.exists(extracted_location), info = glue("Extracted location should exist: {extracted_location}"))
#   #     extracted_dirs <- list.dirs(extracted_location, full.names = TRUE)
#   #   })
    
#   #   required_cols <- c("study_name", "data_type", "variant_type", "ancestry", 
#   #                     "sample_size", "category", "cis_trans")
#   #   for (col in required_cols) {
#   #     expect_true(col %in% colnames(studies), info = glue("Column {col} should exist in studies_processed"))
#   #   }
    
#   #   # Check data types are valid
#   #   expect_true(all(studies$data_type %in% c("phenotype", "gene_expression", "protein", 
#   #                                           "methylation", "metabolome", "cell_trait", 
#   #                                           "plasma_protein", "splice_variant", "transcript")),
#   #               info = "All data_type values should be valid")
    
#   #   expect_true(all(studies$variant_type %in% c("common", "rare_exome", "rare_wgs")),
#   #               info = "All variant_type values should be valid")
    
#   #   expect_true(all(studies$ancestry %in% c("EUR", "EAS", "AFR", "SAS")),
#   #               info = "All ancestry values should be valid")
    
#   #   # Check sample sizes are reasonable
#   #   expect_true(all(studies$sample_size > 0, na.rm = TRUE),
#   #               info = "All sample sizes should be positive")
    
#   #   expect_true(nrow(studies) > 0, info = "Should have at least one study")
#   # })

#   # test_that("study_extractions are valid", {
#   #   extractions_file <- file.path(results_dir, "current/study_extractions.tsv.gz")
#   #   expect_true(file.exists(extractions_file), info = glue("File should exist: {extractions_file}"))
#   #   expect_gt(file.size(extractions_file), 0, info = glue("File should not be empty: {extractions_file}"))

#   #   studies <- vroom::vroom(file.path(results_dir, "current/studies_processed.tsv.gz"), show_col_types = FALSE)
#   #   extractions <- vroom::vroom(extractions_file, show_col_types = FALSE)
#   #   expect_true(nrow(extractions) == nrow(studies), info = "Number of study extractions should be equal to number of studies")
#   #   expect_true(all(extractions$study %in% studies$study_name), info = "All studies should be in study_extractions")

#   #   sapply(extractions$file, function(file_path) {
#   #     expect_true(file.exists(file_path), info = glue("Extracted file should exist: {file_path}"))
#   #   })
    
#   #   required_cols <- c("study_name", "data_type", "variant_type", "ancestry", 
#   #                     "sample_size", "category", "cis_trans")
#   # })
  
#   # test_that("Study directories contain expected files", {
#   #   study_dirs <- list.dirs(file.path(data_dir, "study"), recursive = FALSE)
#   #   expect_true(length(study_dirs) > 0, info = "Should have at least one study directory")
#   #   for (study_dir in study_dirs) {
#   #     expected_files <- c("extracted_snps.tsv", "standardised_snps.tsv", "imputed_snps.tsv", "finemapped_snps.tsv")
#   #     for (file_name in expected_files) {
#   #       file_path <- file.path(study_dir, file_name)
#   #       expect_true(file.exists(file_path), info = glue("File should exist: {file_path}"))
#   #     }
#   #   }
#   # })
  
#   # test_that("LD block data is properly structured", {
#   #   ld_blocks <- vroom::vroom('../pipeline_steps/data/ld_blocks.tsv')
#   #   ld_info <- construct_ld_block(ld_blocks$ancestry, ld_blocks$chr, ld_blocks$start, ld_blocks$stop)
#   #   ld_info <- dplyr::filter(ld_info, dir.exists(ld_block_data))

#   #   for (ld_block in ld_info$ld_block) {
#   #     ld_block_dir <- ld_info$ld_block_data[ld_info$ld_block == ld_block]
#   #     expected_files <- c("extracted_studies.tsv",
#   #       "standardised_studies.tsv",
#   #       "imputed_studies.tsv",
#   #       "finemapped_studies.tsv",
#   #       "coloc_pairwise_results.tsv.gz",
#   #       "coloc_clustered_results.tsv.gz",
#   #       "rare_results.tsv.gz"
#   #     )
#   #     for (file_name in expected_files) {
#   #       file_path <- file.path(ld_block_dir, file_name)
#   #       expect_true(file.exists(file_path), info = glue("File should exist: {file_path}"))
#   #       results <- vroom::vroom(file_path, show_col_types = FALSE)
#   #       expect_true(nrow(results) > 0, info = glue("File should have rows: {file_path}"))
#   #     }
#   #   }
#   # })
  
#   # # Test 5: Check that study extractions exist
  
#   # # Test 6: Check coloc results structure
#   # test_that("Coloc results are valid", {
#   #   coloc_file <- file.path(results_dir, "latest/coloc_results.tsv.gz")
    
#   #   if (file.exists(coloc_file)) {
#   #     coloc <- vroom::vroom(coloc_file, show_col_types = FALSE)
      
#   #     # Check required columns for coloc results
#   #     expected_cols <- c("unique_study_a", "unique_study_b", "ld_block", "h4", "h3", "h2", "h1", "h0")
#   #     for (col in expected_cols) {
#   #       if (col %in% colnames(coloc)) {
#   #         # Check that probability values are between 0 and 1
#   #         if (col %in% c("h4", "h3", "h2", "h1", "h0")) {
#   #           expect_true(all(coloc[[col]] >= 0 & coloc[[col]] <= 1, na.rm = TRUE),
#   #                       info = glue("Column {col} should contain probabilities between 0 and 1"))
#   #         }
#   #       }
#   #     }
      
#   #     if (nrow(coloc) > 0) {
#   #       # Check that probabilities sum to approximately 1
#   #       prob_cols <- intersect(c("h4", "h3", "h2", "h1", "h0"), colnames(coloc))
#   #       if (length(prob_cols) > 0) {
#   #         prob_sums <- rowSums(coloc[, prob_cols], na.rm = TRUE)
#   #         expect_true(all(abs(prob_sums - 1) < 0.01, na.rm = TRUE),
#   #                     info = "Coloc probabilities should sum to approximately 1")
#   #       }
#   #     }
#   #   }
#   # })
  
#   # # Test 7: Check that LD block directories contain expected files
#   # test_that("LD block directories contain expected files", {
#   #   ld_block_dirs <- list.dirs(file.path(data_dir, "ld_blocks"), recursive = FALSE)
    
#   #   if (length(ld_block_dirs) > 0) {
#   #     # Check a sample of LD block directories
#   #     sample_dirs <- head(ld_block_dirs, 3)
      
#   #     for (ld_dir in sample_dirs) {
#   #       expected_files <- c(
#   #         "imputed_studies.tsv",
#   #         "finemapped_studies.tsv"
#   #       )
        
#   #       for (file_name in expected_files) {
#   #         file_path <- file.path(ld_dir, file_name)
#   #         if (file.exists(file_path)) {
#   #           # Check file is not empty
#   #           file_size <- file.size(file_path)
#   #           expect_gt(file_size, 0, 
#   #                     info = glue("LD block file should not be empty: {file_path}"))
#   #         }
#   #       }
#   #     }
#   #   }
#   # })
  
#   # # Test 8: Check data quality metrics
#   # test_that("Data quality metrics are reasonable", {
#   #   studies_file <- file.path(results_dir, "latest/studies_processed.tsv.gz")
    
#   #   if (file.exists(studies_file)) {
#   #     studies <- vroom::vroom(studies_file, show_col_types = FALSE)
      
#   #     # Check sample size distribution
#   #     if ("sample_size" %in% colnames(studies)) {
#   #       sample_sizes <- studies$sample_size[!is.na(studies$sample_size)]
#   #       if (length(sample_sizes) > 0) {
#   #         expect_true(median(sample_sizes) > 100, 
#   #                     info = "Median sample size should be reasonable (>100)")
#   #         expect_true(max(sample_sizes) < 10000000, 
#   #                     info = "Maximum sample size should be reasonable (<10M)")
#   #       }
#   #     }
      
#   #     # Check data type distribution
#   #     if ("data_type" %in% colnames(studies)) {
#   #       data_type_counts <- table(studies$data_type)
#   #       expect_true(length(data_type_counts) > 0, 
#   #                   info = "Should have studies from multiple data types")
#   #     }
#   #   }
#   # })
  
#   # # Test 9: Check file timestamps are recent
#   # test_that("Result files have recent timestamps", {
#   #   result_files <- c(
#   #     file.path(results_dir, "latest/studies_processed.tsv.gz"),
#   #     file.path(results_dir, "latest/traits_processed.tsv.gz"),
#   #     file.path(results_dir, "latest/study_extractions.tsv.gz")
#   #   )
    
#   #   for (file_path in result_files) {
#   #     if (file.exists(file_path)) {
#   #       file_time <- file.info(file_path)$mtime
#   #       time_diff <- as.numeric(Sys.time() - file_time, units = "hours")
#   #       expect_lt(time_diff, 2, 
#   #                 info = glue("File {file_path} should be recently created (<2 hours old)"))
#   #     }
#   #   }
#   # })
  
#   # # Test 10: Check pipeline logs for errors
#   # test_that("Pipeline logs do not contain critical errors", {
#   #   log_file <- file.path(data_dir, "pipeline_metadata/logs/snakemake.log")
    
#   #   if (file.exists(log_file)) {
#   #     log_content <- readLines(log_file, warn = FALSE)
      
#   #     error_pattern <- "rror"
      
#   #     matches <- grep(error_pattern, log_content, ignore.case = TRUE)
#   #     if (length(matches) > 0) {
#   #       critical_errors <- log_content[matches]
#   #     }
      
#   #     # Allow some expected errors but flag unexpected ones
#   #     if (length(critical_errors) > 0) {
#   #       cat("Found potential errors in pipeline log:\n")
#   #       for (error in head(critical_errors, 5)) {
#   #         cat(glue("  {error}\n"))
#   #       }
#   #       # Don't fail the test for now, just warn
#   #       cat("Warning: Critical errors found in pipeline log\n")
#   #     }
#   #   }
#   # })
  
#   cat(glue("Pipeline execution time: {round(execution_time, 2)} seconds\n"))
#   cat(glue("Test results directory: {results_dir}\n"))
# })