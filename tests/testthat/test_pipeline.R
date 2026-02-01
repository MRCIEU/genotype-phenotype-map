library(testthat)

setwd('../../')
source('pipeline_steps/constants.R')

#If you want to add more data to this test pipeline
#smr --beqtl-summary /some/besd_file  --extract-probe probe.list  --query 1 --make-besd --out subset_of_besd_file
#where probe.list has the format:

total_studies <- 8
ld_blocks <- vroom::vroom('tests/data/ld_blocks.tsv', show_col_types = FALSE)
ld_info <- construct_ld_block(ld_blocks$ancestry, ld_blocks$chr, ld_blocks$start, ld_blocks$stop)
coloc_block <- 'EUR/1/1170341-1730405'
rare_result_block <- 'EUR/1/1730405-3355587'

test_that("Identify studies to process", {
  output <- system(glue::glue("Rscript pipeline_steps/identify_studies_to_process.R"),
    wait = TRUE,
    intern = TRUE,
    ignore.stdout = FALSE,
    ignore.stderr = FALSE
  )
  has_errors <- grepl("error", output, ignore.case = TRUE)
  print(has_errors)
  expect_false(any(has_errors), info = "Identify studies to process should not contain errors")

  studies_to_process_file <- glue::glue('{pipeline_metadata_dir}/studies_to_process.tsv')
  expect_true(file.exists(studies_to_process_file), info = "Studies to process file should exist")
  studies_to_process <- vroom::vroom(studies_to_process_file, show_col_types = FALSE)
  expect_equal(nrow(studies_to_process), total_studies, info = "Studies to process file should have {total_studies} studies")
})

test_that("Pipeline execution and file validation", {
  snakemake_log <- glue::glue('{data_dir}pipeline_metadata/logs/snakemake.log')
  system(glue::glue("snakemake --profile ./ > {snakemake_log} 2>&1"), wait = TRUE)
  output <- readLines(snakemake_log)
  has_errors <- grepl("error", output, ignore.case = TRUE)
  if (any(has_errors)) {
    error_file <- glue::glue('{data_dir}pipeline_metadata/logs/snakemake_error.log')
    writeLines(output, error_file)
    print(glue::glue("Errors written to {error_file}"))
    stop()
  }
  expect_false(any(has_errors), info = "Pipeline execution should not contain errors")

  expected_static_web_files <- list(
    opengwas_ids_file=glue::glue('{static_web_dir}opengwas_ids.json'),
    phenotype_id_map_file=glue::glue('{static_web_dir}phenotype_id_map.json'),
    robots_txt_file=glue::glue('{static_web_dir}robots.txt'),
    sitemap_xml_file=glue::glue('{static_web_dir}sitemap.xml'),
    pipeline_summary_file=glue::glue('{static_web_dir}pipeline_summary.html')
  )

  expected_tsv_files <- list(
    studies_file=file.path(current_results_dir, "studies_processed.tsv.gz"),
    traits_file=file.path(current_results_dir, "traits_processed.tsv.gz"),
    study_extractions_file=file.path(current_results_dir, "study_extractions.tsv.gz"),
    coloc_clustered_results_file=file.path(current_results_dir, "coloc_clustered_results.tsv.gz"),
    coloc_pairwise_results_file=file.path(current_results_dir, "coloc_pairwise_results.tsv.gz"),
    rare_results_file=file.path(current_results_dir, "rare_results.tsv.gz")
  )

  expected_db_files <- list(
    studies_db_file=file.path(current_results_dir, "studies.db"),
    associations_db_file=file.path(current_results_dir, "associations.db"),
    coloc_pairs_full_db_file=file.path(current_results_dir, "coloc_pairs_full.db"),
    coloc_pairs_significant_db_file=file.path(current_results_dir, "coloc_pairs.db"),
    ld_db_file=file.path(current_results_dir, "ld.db"),
    gwas_upload_db_file=file.path(current_results_dir, "gwas_upload.db")
  )
  
  test_that("Main result files are created", {
    for (file_path in expected_tsv_files) {
      expect_true(file.exists(file_path), info = glue::glue("File should exist: {file_path}"))
      
      file_size <- file.size(file_path)
      expect_true(file_size > 0, info = glue::glue("File should not be empty: {file_path}"))
    }
  })

  test_that("traits_processed and studies_processed are valid", {
    expect_true(file.exists(expected_tsv_files$traits_file), info = glue::glue("File should exist: {expected_tsv_files$traits_file}"))
    expect_true(file.size(expected_tsv_files$traits_file) > 0, info = glue::glue("File should not be empty: {expected_tsv_files$traits_file}"))

    expect_true(file.exists(expected_tsv_files$studies_file), info = glue::glue("File should exist: {expected_tsv_files$studies_file}"))
    expect_true(file.size(expected_tsv_files$studies_file) > 0, info = glue::glue("File should not be empty: {expected_tsv_files$studies_file}"))
    
    studies <- vroom::vroom(expected_tsv_files$studies_file, show_col_types = FALSE)
    traits <- vroom::vroom(expected_tsv_files$traits_file, show_col_types = FALSE)

    expect_equal(nrow(studies), nrow(traits), info = "Number of studies and traits should be equal")
    expect_true(nrow(studies) > 0, info = "Should have at least one study")
    expect_true(all(studies$study_name %in% traits$study_name), info = "All studies should be in traits_processed")
    expect_equal(nrow(studies), total_studies, info = "Should have {total_studies} studies")
    expect_equal(nrow(traits), total_studies, info = "Should have {total_studies} traits")

    joined_studies_traits <- dplyr::inner_join(studies, traits, by = "study_name")
    expect_true(nrow(joined_studies_traits) == nrow(studies) && nrow(joined_studies_traits) > 0, info = "All studies should have an associated trait")

    sapply(joined_studies_traits$extracted_location, function(extracted_location) {
      expect_true(dir.exists(extracted_location), info = glue::glue("Extracted location should exist: {extracted_location}"))
      extracted_dirs <- list.dirs(extracted_location, full.names = TRUE)
    })
    
    
    expect_true(all(studies$data_type %in% names(data_types)), info = "All data_type values should be valid")
    expect_true(all(studies$variant_type %in% names(variant_types)), info = "All variant_type values should be valid")
    expect_true(all(studies$ancestry %in% names(ancestry_map)), info = "All ancestry values should be valid")
    expect_true(all(studies$sample_size > 0, na.rm = TRUE), info = "All sample sizes should be positive")
    expect_true(nrow(studies) > 0, info = "Should have at least one study")
  })

  test_that("study_extractions are valid", {
    expect_true(file.exists(expected_tsv_files$study_extractions_file), info = glue::glue("File should exist: {expected_tsv_files$study_extractions_file}"))
    expect_true(file.size(expected_tsv_files$study_extractions_file) > 0, info = glue::glue("File should not be empty: {expected_tsv_files$study_extractions_file}"))

    studies <- vroom::vroom(expected_tsv_files$studies_file, show_col_types = FALSE)
    extractions <- vroom::vroom(expected_tsv_files$study_extractions_file, show_col_types = FALSE)
    expect_true(all(extractions$study %in% studies$study_name), info = "All studies should be in study_extractions")

    common_extractions <- dplyr::left_join(extractions, studies, by = c("study"="study_name")) |>
      dplyr::filter(variant_type == variant_types$common)
    lapply(seq_len(nrow(common_extractions)), function(i) {
      extraction <- common_extractions[i, ]
      error_msg <- glue::glue("Min p value should be less than {lowest_p_value_threshold}: {extraction$min_p}")
      expect_true(file.exists(extraction$file), info = glue::glue("Extracted file should exist: {extraction$file}"))
      expect_true(file.size(extraction$file) > 0, info = glue::glue("Extracted file should not be empty: {extraction$file}"))
      expect_true(file.exists(extraction$svg_file), info = glue::glue("Extracted file should exist: {extraction$svg_file}"))
      expect_true(file.size(extraction$svg_file) > 0, info = glue::glue("Extracted file should not be empty: {extraction$svg_file}"))
      expect_true(file.exists(extraction$file_with_lbfs), info = glue::glue("Extracted file should exist: {extraction$file_with_lbfs}"))
      expect_true(file.size(extraction$file_with_lbfs) > 0, info = glue::glue("Extracted file should not be empty: {extraction$file_with_lbfs}"))
    })
    
  })
  
  test_that("LD Block directories contain expected files", {
    coloc_block_dir <- glue::glue('{ld_block_data_dir}{coloc_block}')
    expect_true(dir.exists(coloc_block_dir), info = glue::glue("LD block directory should exist: {coloc_block_dir}"))
    expected_common_files <- c("extracted_studies.tsv", "standardised_studies.tsv", "imputed_studies.tsv", "finemapped_studies.tsv", "coloc_pairwise_results.tsv.gz", "coloc_clustered_results.tsv.gz")

    for (file_name in expected_common_files) {
      file_path <- file.path(coloc_block_dir, file_name)
      expect_true(file.exists(file_path), info = glue::glue("File should exist: {file_path}"))
      expect_true(file.size(file_path) > 0, info = glue::glue("File should not be empty: {file_path}"))
    }

    rare_block_dir <- glue::glue('{ld_block_data_dir}{rare_result_block}')
    expect_true(dir.exists(rare_block_dir), info = glue::glue("LD block directory should exist: {rare_block_dir}"))
    expected_rare_files <- c("extracted_studies.tsv", "standardised_studies.tsv", "compare_rare_results.tsv", "compare_rare_cached_studies.tsv")
    for (file_name in expected_rare_files) {
      file_path <- file.path(rare_block_dir, file_name)
      expect_true(file.exists(file_path), info = glue::glue("File should exist: {file_path}"))
      expect_true(file.size(file_path) > 0, info = glue::glue("File should not be empty: {file_path}"))
    }
  })
  
  test_that("Study directory data is properly structured", {
    studies <- vroom::vroom(expected_tsv_files$studies_file, show_col_types = FALSE)
    for (i in seq_len(nrow(studies))) {
      study <- studies[i, ]
      study_dir <- study$extracted_location
      expected_dirs <- c("extracted",
        "standardised",
        "imputed",
        "finemapped"
      )
      if (study$data_type == data_types$phenotype) {
        expected_dirs <- c(expected_dirs, "svgs")
      }
      for (dir_name in expected_dirs) {
        dir_path <- file.path(study_dir, dir_name)
        expect_true(dir.exists(dir_path), info = glue::glue("Directory should exist"))
      }
      extracted_snps <- vroom::vroom(file.path(study_dir, "extracted_snps.tsv"), show_col_types = FALSE)
      for (j in seq_len(nrow(extracted_snps))) {
        file <- extracted_snps$file[j]
        expect_true(file.exists(file), info = glue::glue("File should exist: {file}"))
        expect_true(file.size(file) > 0, info = glue::glue("File should not be empty: {file}"))
      }
      files_in_extracted_dir <- list.files(file.path(study_dir, "extracted"), full.names = TRUE)
      expect_true(nrow(extracted_snps) == length(files_in_extracted_dir), info = glue::glue("Extracted SNPs file should not be empty"))
    }
  })
  
  test_that("Coloc results are valid", {
    clustered_results_file <- file.path(ld_block_data_dir, coloc_block, "coloc_clustered_results.tsv.gz")
    pairwise_results_file <- file.path(ld_block_data_dir, coloc_block, "coloc_pairwise_results.tsv.gz")
    expect_true(file.exists(clustered_results_file), info = glue::glue("File should exist: {clustered_results_file}"))
    expect_true(file.exists(pairwise_results_file), info = glue::glue("File should exist: {pairwise_results_file}"))
    
    clustered_results <- vroom::vroom(clustered_results_file, show_col_types = FALSE)
    pairwise_results <- vroom::vroom(pairwise_results_file, show_col_types = FALSE)
    
    expect_true(nrow(clustered_results) > 0, info = "Clustered results should not be empty")
  })

  test_that("Rare results are valid", {
    expect_true(file.exists(expected_tsv_files$rare_results_file), info = glue::glue("File should exist: {expected_tsv_files$rare_results_file}"))
    rare_results <- vroom::vroom(expected_tsv_files$rare_results_file, show_col_types = FALSE)
    expect_true(nrow(rare_results) > 0, info = "Rare results should not be empty")
  })

  test_that("Static web files are ready", {
    expect_true(file.exists(expected_static_web_files$pipeline_summary_file), info = glue::glue("File should exist: {expected_static_web_files$pipeline_summary_file}"))
    expect_true(file.size(expected_static_web_files$pipeline_summary_file) > 0, info = glue::glue("File should not be empty: {expected_static_web_files$pipeline_summary_file}"))
    expect_true(file.exists(expected_static_web_files$opengwas_ids_file), info = glue::glue("File should exist: {expected_static_web_files$opengwas_ids_file}"))
    expect_true(file.size(expected_static_web_files$opengwas_ids_file) > 0, info = glue::glue("File should not be empty: {expected_static_web_files$opengwas_ids_file}"))
  })

  test_that("DB files are valid", {
    source('pipeline_steps/database_definitions.R')
    expect_true(file.exists(expected_db_files$studies_db_file), info = glue::glue("File should exist: {expected_db_files$studies_db_file}"))
    expect_true(file.size(expected_db_files$studies_db_file) > 0, info = glue::glue("File should not be empty: {expected_db_files$studies_db_file}"))
    studies_conn <- duckdb::dbConnect(duckdb::duckdb(), expected_db_files$studies_db_file)
    for (table in studies_db) {
      table_data <- DBI::dbGetQuery(studies_conn, glue::glue("SELECT count(*) as count FROM {table$name}"))
      expect_true(table_data$count > 0, info = glue::glue("Table should not be empty: {table$name}"))
    }
    for (table in additional_studies_tables) {
      table_data <- DBI::dbGetQuery(studies_conn, glue::glue("SELECT count(*) as count FROM {table$name}"))
      expect_true(table_data$count > 0, info = glue::glue("Table should not be empty: {table$name}"))
    }
    ld_conn <- duckdb::dbConnect(duckdb::duckdb(), expected_db_files$ld_db_file)
    table_data <- DBI::dbGetQuery(ld_conn, glue::glue("SELECT count(*) as count FROM {ld_table$name}"))
    expect_true(table_data$count > 0, info = glue::glue("Table should not be empty: {ld_table$name}"))

    coloc_pairs_full_conn <- duckdb::dbConnect(duckdb::duckdb(), expected_db_files$coloc_pairs_full_db_file)
    table_data <- DBI::dbGetQuery(coloc_pairs_full_conn, glue::glue("SELECT count(*) as count FROM {coloc_pairs_full_table$name}"))
    expect_true(table_data$count > 0, info = glue::glue("Table should not be empty: {coloc_pairs_full_table$name}"))

    coloc_pairs_significant_conn <- duckdb::dbConnect(duckdb::duckdb(), expected_db_files$coloc_pairs_significant_db_file)
    table_data <- DBI::dbGetQuery(coloc_pairs_significant_conn, glue::glue("SELECT count(*) as count FROM {coloc_pairs_significant_table$name}"))
    expect_true(table_data$count > 0, info = glue::glue("Table should not be empty: {coloc_pairs_significant_table$name}"))

    gwas_upload_conn <- duckdb::dbConnect(duckdb::duckdb(), expected_db_files$gwas_upload_db_file)
    for (table in gwas_upload_db) {
      table_data <- DBI::dbGetQuery(gwas_upload_conn, glue::glue("SELECT count(*) as count FROM {table$name}"))
      expect_true(table_data$count == 0, info = glue::glue("Table should be empty"))
    }
  })
})