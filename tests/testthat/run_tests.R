library(testthat)

TEST_DIR <- "./tests/testthat"
OUTPUT_FILE_PATH <- "./tests/testing_complete.txt"


#Test setup: set env vars if not already set
Sys.setenv("TEST_RUN" = "test")
Sys.setenv("TIMESTAMP" = "test")
Sys.setenv("DATA_DIR" = "/local-scratch/projects/genotype-phenotype-map/test/data/")
Sys.setenv("RESULTS_DIR" = "/local-scratch/projects/genotype-phenotype-map/test/results/")
Sys.setenv("BACKUP_DIR" = "/local-scratch/projects/genotype-phenotype-map/test/backup/")
Sys.setenv("GWAS_UPLOAD_DIR" = "/local-scratch/projects/genotype-phenotype-map/test/data/gwas_upload/")

source('pipeline_steps/constants.R')

#Cleanup previous test run
# system(glue::glue('rm -r {data_dir}pipeline_metadata/studies_to_process.tsv'), ignore.stdout = TRUE, ignore.stderr = TRUE)
# system(glue::glue('rm -r {data_dir}study/*'), ignore.stdout = TRUE, ignore.stderr = TRUE)
# system(glue::glue('rm -r {data_dir}ld_blocks/*/*'), ignore.stdout = TRUE, ignore.stderr = TRUE)
# system(glue::glue('rm -r {gwas_upload_dir}*'), ignore.stdout = TRUE, ignore.stderr = TRUE)
# system(glue::glue('rm -r {data_dir}pipeline_metadata/updated_ld_blocks_to_colocalise.tsv'), ignore.stdout = TRUE, ignore.stderr = TRUE)
# system(glue::glue('rm -r {results_dir}latest/studies_processed.tsv.gz'), ignore.stdout = TRUE, ignore.stderr = TRUE)
# system(glue::glue('rm -r {results_dir}/*'), ignore.stdout = TRUE, ignore.stderr = TRUE)
dir.create(current_results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(latest_results_dir, showWarnings = FALSE, recursive = TRUE)

message(paste("Starting tests in:", normalizePath(TEST_DIR)))
tryCatch({
  testthat::test_dir(TEST_DIR, reporter = "progress", stop_on_failure = TRUE)

  branch_name <- trimws(system("git rev-parse --abbrev-ref HEAD", intern = TRUE)[1])
  status_message <- paste("SUCCESS: All tests passed on branch:", branch_name)
  writeLines(status_message, con = OUTPUT_FILE_PATH)
  message("\n✅ SUCCESS: All tests passed.")
  message(status_message)
}, error = function(e) {
  message("\n❌ FAILURE: Some tests failed or had errors.")
})