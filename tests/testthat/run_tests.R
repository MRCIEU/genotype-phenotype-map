library(testthat)
library(argparser)

parser <- argparser::arg_parser("Run tests")
parser <- argparser::add_argument(
  parser,
  "--only",
  help = "Run only one test file: 'worker' or 'pipeline'",
  type = "character",
  default = NA
)
args <- argparser::parse_args(parser)

only_test <- if (!is.na(args$only) && nchar(trimws(args$only)) > 0) {
  match.arg(trimws(args$only), c("worker", "pipeline"))
} else {
  NULL
}

TEST_DIR <- "./tests/testthat"
OUTPUT_FILE_PATH <- "./tests/testing_complete.txt"


# Test setup: set env vars if not already set
Sys.setenv("TEST_RUN" = "test")
Sys.setenv("TIMESTAMP" = "test")
Sys.setenv("DATA_DIR" = "/local-scratch/projects/genotype-phenotype-map/test/data/")
Sys.setenv("RESULTS_DIR" = "/local-scratch/projects/genotype-phenotype-map/test/results/")
Sys.setenv("BACKUP_DIR" = "/local-scratch/projects/genotype-phenotype-map/test/backup/")

source("pipeline_steps/constants.R")

# Cleanup previous test run
system(
  glue::glue("rm -r {data_dir}pipeline_metadata/studies_to_process.tsv"),
  ignore.stdout = TRUE,
  ignore.stderr = TRUE
)
system(glue::glue("rm -r {data_dir}study/*"), ignore.stdout = TRUE, ignore.stderr = TRUE)
system(glue::glue("rm -r {data_dir}ld_blocks/*/*"), ignore.stdout = TRUE, ignore.stderr = TRUE)
system(
  glue::glue("rm -r {data_dir}pipeline_metadata/updated_ld_blocks_to_colocalise.tsv"),
  ignore.stdout = TRUE,
  ignore.stderr = TRUE
)
system(glue::glue("rm -r {results_dir}latest/studies_processed.tsv.gz"), ignore.stdout = TRUE, ignore.stderr = TRUE)
system(glue::glue("rm -r {results_dir}/*"), ignore.stdout = TRUE, ignore.stderr = TRUE)
system(glue::glue("rm -r {gwas_upload_dir}gwas_upload/*"), ignore.stdout = TRUE, ignore.stderr = TRUE)
system(glue::glue("rm -r {data_dir}ld_blocks/gwas_upload/*"), ignore.stdout = TRUE, ignore.stderr = TRUE)
dir.create(current_results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(latest_results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(results_analysis_dir, showWarnings = FALSE, recursive = TRUE)

message(paste("Starting tests in:", normalizePath(TEST_DIR)))
if (!is.null(only_test)) {
  message(paste("Running only:", only_test))
}
tryCatch(
  {
    if (!is.null(only_test)) {
      test_file <- file.path(TEST_DIR, paste0("test_", only_test, ".R"))
      testthat::test_file(test_file, reporter = "progress", stop_on_failure = TRUE)
    } else {
      testthat::test_dir(TEST_DIR, reporter = "progress", stop_on_failure = TRUE)
    }

    branch_name <- trimws(system("git rev-parse --abbrev-ref HEAD", intern = TRUE)[1])
    status_message <- paste("SUCCESS: All tests passed on branch:", branch_name)
    writeLines(status_message, con = OUTPUT_FILE_PATH)
    message("\n✅ SUCCESS: All tests passed.")
    message(status_message)
  },
  error = function(e) {
    message("\n❌ FAILURE: Some tests failed or had errors.")
    return()
  }
)
