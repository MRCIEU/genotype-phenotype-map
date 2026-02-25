#!/usr/bin/env Rscript
# Lint all R files and show results
# Usage: Rscript scripts/lint.R

library(lintr)

# Get all R files recursively, excluding certain directories
r_files <- list.files(
  path = ".",
  pattern = "\\.R$",
  recursive = TRUE,
  full.names = TRUE
)

# Exclude certain directories
exclude_dirs <- c(
  "docs",
  "data_formatting",
  "pre_steps/vep_annotation",
  "tests/testthat"
)

# Filter out excluded directories
r_files <- r_files[!grepl(paste(exclude_dirs, collapse = "|"), r_files)]

if (length(r_files) == 0) {
  message("No R files found to lint.")
  quit(status = 0)
}

message(paste("Linting", length(r_files), "R files..."))

# Run lintr on all files
lint_results <- list()

for (file in r_files) {
  file_lints <- lintr::lint(file)
  if (length(file_lints) > 0) {
    lint_results <- c(lint_results, file_lints)
  }
}

if (length(lint_results) > 0) {
  print(lint_results)
  message("\n")
  message(paste("Found", length(lint_results), "linting issue(s)."))
} else {
  message("No linting issues found!")
}
