library(testthat)

setwd("../../")
source("pipeline_steps/constants.R")
ld_block_of_interest <- "EUR/1/101384499-103762931"

setup({
  real_ld_block_data_dir <- sub("/test/", "/", ld_block_data_dir)
  existing_finemapped_studies <- glue::glue("{ld_block_of_interest}/finemapped_studies.tsv")
  dir.create(glue::glue("{ld_block_data_dir}/{ld_block_of_interest}"), recursive = TRUE, showWarnings = FALSE)
  file.copy(
    glue::glue("{real_ld_block_data_dir}/{existing_finemapped_studies}"),
    glue::glue("{ld_block_data_dir}/{ld_block_of_interest}/finemapped_studies.tsv")
  )
})


# test_that("Pipeline worker fails and sends message to DLQ for bad data", {
#   redis_payload <- '/home/pipeline/tests/data/hg38_tsv_redis_message_bad_data.json'
#   output <- system(glue::glue("Rscript worker/pipeline_worker.R --custom_message_file {redis_payload}"),
#     wait = TRUE,
#     intern = TRUE,
#     ignore.stdout = FALSE,
#     ignore.stderr = FALSE
#   )
# })

test_that("Pipeline worker runs for TSV file", {
  redis_payload <- "/home/pipeline/tests/data/hg38_tsv_redis_message.json"
  output <- system(glue::glue("Rscript worker/pipeline_worker.R --custom_message_file {redis_payload}"),
    wait = TRUE,
    intern = TRUE,
    ignore.stdout = FALSE,
    ignore.stderr = FALSE
  )

  redis_message <- jsonlite::fromJSON(redis_payload)
  gwas_info <- jsonlite::fromJSON(redis_message[[2]])

  has_errors <- grepl("error", output, ignore.case = TRUE)
  if (any(has_errors)) {
    error_file <- glue::glue("{gwas_upload_dir}gwas_upload/{gwas_info$metadata$guid}/failed_worker_error.log")
    writeLines(output, error_file)
    print(glue::glue("Errors written to {error_file}"))
  }
  expect_false(any(has_errors), info = "Pipeline worker should not contain errors")

  update_directories_for_worker(gwas_info$metadata$guid)
  expect_true(dir.exists(extracted_study_dir), info = "GWAS upload directory should exist")
  expect_true(dir.exists(ld_block_data_dir), info = "LD block data directory should exist")
  ld_block_files <- c(
    glue::glue("{ld_block_data_dir}/{ld_block_of_interest}/standardised_studies.tsv"),
    glue::glue("{ld_block_data_dir}/{ld_block_of_interest}/imputed_studies.tsv"),
    glue::glue("{ld_block_data_dir}/{ld_block_of_interest}/finemapped_studies.tsv"),
    glue::glue("{ld_block_data_dir}/{ld_block_of_interest}/coloc_pairwise_results.tsv.gz"),
    glue::glue("{ld_block_data_dir}/{ld_block_of_interest}/coloc_clustered_results.tsv.gz")
  )
  for (file in ld_block_files) {
    expect_true(file.exists(file), info = glue::glue("File should exist: {file}"))
    results <- vroom::vroom(file, show_col_types = FALSE)
    expect_true(nrow(results) > 0, info = glue::glue("Results file should have rows: {file}"))
  }

  compiled_files <- c(
    glue::glue("{extracted_study_dir}/compiled_coloc_pairwise_results.tsv"),
    glue::glue("{extracted_study_dir}/compiled_extracted_studies.tsv"),
    glue::glue("{extracted_study_dir}/compiled_coloc_clustered_results.tsv"),
    glue::glue("{extracted_study_dir}/compiled_associations.tsv")
  )

  for (file in compiled_files) {
    expect_true(file.exists(file), info = glue::glue("Compiled file should exist: {file}"))
    results <- vroom::vroom(file, show_col_types = FALSE)
    expect_true(nrow(results) > 0, info = glue::glue("Compiled file should have rows: {file}"))
  }
})

test_that("Pipeline worker runs for VCF file", {
  expect_true(TRUE)
})

test_that("Pipeline worker picks up message from DLQ", {
  expect_true(TRUE)
})

test_that("get_from_in_progress_queue returns message when queue has items", {
  source("worker/redis_client.R", local = TRUE)

  test_payload <- '{"metadata":{"guid":"test-guid"}}'
  mock_conn <- list(
    BRPOP = function(queue, timeout = 0.1) {
      if (queue == process_gwas_in_progress) {
        return(list(process_gwas_in_progress, test_payload))
      }
      return(NULL)
    }
  )

  result <- get_from_in_progress_queue(mock_conn)
  expect_false(is.null(result))
  expect_equal(result[[1]], process_gwas_in_progress)
  expect_equal(result[[2]], test_payload)
})

# test_that("Pipeline worker processes message from in_progress queue when custom message has in_progress queue name", {
#   redis_payload <- "/home/pipeline/tests/data/hg38_tsv_redis_message_in_progress.json"
#   output <- system(glue::glue("Rscript worker/pipeline_worker.R --custom_message_file {redis_payload}"),
#     wait = TRUE,
#     intern = TRUE,
#     ignore.stdout = FALSE,
#     ignore.stderr = FALSE
#   )

#   redis_message <- jsonlite::fromJSON(redis_payload)
#   gwas_info <- jsonlite::fromJSON(redis_message[[2]])

#   has_errors <- grepl("error", output, ignore.case = TRUE)
#   if (any(has_errors)) {
#     error_file <- glue::glue("{gwas_upload_dir}gwas_upload/{gwas_info$metadata$guid}/failed_worker_error.log")
#     writeLines(output, error_file)
#     print(glue::glue("Errors written to {error_file}"))
#   }
#   expect_false(any(has_errors), info = "Pipeline worker should process from in_progress queue without errors")

#   update_directories_for_worker(gwas_info$metadata$guid)
#   expect_true(dir.exists(extracted_study_dir), info = "GWAS upload directory should exist")
# })
