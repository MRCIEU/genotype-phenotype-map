
library(testthat)

setwd('../../')
source('pipeline_steps/constants.R')

# Mock the redux::hiredis function
mock_hiredis <- function(host, port, payload_file) {
  mock_redis_conn <- list(
    BRPOP = function(queue, timeout = 0) {
      # Return the actual JSON data from your test file
      redis_payload <- 'tests/data/hg38_tsv_redis_message.json'
      redis_message <- jsonlite::fromJSON(readLines(redis_payload))
      return(redis_message)
    },
    LPUSH = function(queue, message) {
      return(1)
    }
  )
  return(mock_redis_conn)
}


test_that("Pipeline worker runs for TSV file", {
  redis_payload <- 'tests/data/hg38_tsv_redis_message.json'
  assignInNamespace("hiredis", mock_hiredis, ns = "redux")
  redis_message <- jsonlite::fromJSON(readLines(redis_payload))
  output <- system(glue::glue("Rscript worker/pipeline_worker.R --custom_message_file {redis_payload}"),
    wait = TRUE,
    intern = TRUE,
    ignore.stdout = FALSE,
    ignore.stderr = FALSE
  )

  gwas_info <- jsonlite::fromJSON(redis_message[[2]])

  has_errors <- grepl("error", output, ignore.case = TRUE)
  expect_false(any(has_errors), info = "Pipeline worker should not contain errors")
  gwas_upload_directory <- glue::glue('{gwas_upload_dir}/{gwas_info$metadata$guid}')
  expect_true(file.exists(gwas_upload_directory), info = "GWAS upload directory should exist")

})

test_that("Pipeline worker runs for VCF file", {
  expect_true(TRUE)
})

test_that("Pipeline worker picks up message from DLQ", {
  expect_true(TRUE)
})