library(testthat)
library(argparser)

parser <- argparser::arg_parser("Run tests")
parser <- argparser::add_argument(
  parser,
  "--only",
  help = "Run only one test suite: 'worker', 'pipeline' or 'unit'",
  type = "character",
  default = NA
)
parser <- argparser::add_argument(
  parser,
  "--dont_delete",
  flag = TRUE,
  help = "Don't delete the test data"
)
args <- argparser::parse_args(parser)

only_test <- if (!is.na(args$only) && nchar(trimws(args$only)) > 0) {
  match.arg(trimws(args$only), c("worker", "pipeline", "unit"))
} else {
  NULL
}

TEST_DIR <- "./tests/testthat"
UNIT_TEST_DIR <- "/local-scratch/projects/genotype-phenotype-map/test/unit/"
OUTPUT_FILE_PATH <- "./tests/testing_complete.txt"


# Test setup: set env vars if not already set
Sys.setenv("TEST_RUN" = "test")
Sys.setenv("DATA_DIR" = "/local-scratch/projects/genotype-phenotype-map/test/e2e/data/")
Sys.setenv("RESULTS_DIR" = "/local-scratch/projects/genotype-phenotype-map/test/e2e/results/")
Sys.setenv("BACKUP_DIR" = "/local-scratch/projects/genotype-phenotype-map/test/e2e/backup/")

source("pipeline_steps/constants.R")

if (!args$dont_delete) {
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
}

message(paste("Starting tests in:", normalizePath(TEST_DIR)))
test_names <- if (!is.null(only_test)) only_test else c("pipeline", "worker", "unit")

if ("unit" %in% test_names && !args$dont_delete) {
  # Cleanup previous unit test run, leaving the permanent fixtures and reference panel links in place
  unit_dirs_to_clean <- paste0(UNIT_TEST_DIR, c("data/study", "data/ld_blocks", "data/pipeline_metadata", "results"))
  for (unit_dir in unit_dirs_to_clean) {
    unlink(list.files(unit_dir, full.names = TRUE), recursive = TRUE)
  }
}
message(paste("Running in parallel:", paste(test_names, collapse = ", ")))

log_dir <- glue::glue("{data_dir}pipeline_metadata/logs")
dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)

tryCatch(
  {
    # Each test file runs in its own forked process, with output written to its own log file
    test_jobs <- lapply(test_names, function(test_name) {
      log_file <- glue::glue("{log_dir}/test_{test_name}.log")
      job <- parallel::mcparallel(
        {
          options(cli.dynamic = FALSE)
          log_con <- file(log_file, open = "wt")
          sink(log_con)
          sink(log_con, type = "message")
          if (test_name == "unit") {
            testthat::test_dir(file.path(TEST_DIR, "unit"), reporter = "progress", stop_on_failure = TRUE)
          } else {
            test_file <- file.path(TEST_DIR, paste0("test_", test_name, ".R"))
            testthat::test_file(test_file, reporter = "progress", stop_on_failure = TRUE)
          }
        },
        name = test_name
      )
      return(job)
    })
    test_results <- parallel::mccollect(test_jobs)

    failed_tests <- c()
    for (test_name in test_names) {
      log_file <- glue::glue("{log_dir}/test_{test_name}.log")
      message(glue::glue("\n===== {test_name} ({log_file}) ====="))
      message(paste(readLines(log_file), collapse = "\n"))

      result <- test_results[[test_name]]
      if (is.null(result) || inherits(result, "try-error")) {
        failed_tests <- c(failed_tests, test_name)
      }
    }
    if (length(failed_tests) > 0) {
      stop(paste("Failed tests:", paste(failed_tests, collapse = ", ")))
    }

    branch_name <- trimws(system("git rev-parse --abbrev-ref HEAD", intern = TRUE)[1])
    status_message <- paste("SUCCESS: All tests passed on branch:", branch_name)
    writeLines(status_message, con = OUTPUT_FILE_PATH)
    message("\n✅ SUCCESS: All tests passed.")
    message(status_message)
  },
  error = function(e) {
    message("\n❌ FAILURE: Some tests failed or had errors.")
    writeLines(e$message, con = OUTPUT_FILE_PATH)
    return()
  }
)
