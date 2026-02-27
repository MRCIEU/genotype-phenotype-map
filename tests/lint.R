#!/usr/bin/env Rscript
# Lint check script that exits with error code if issues are found
# Usage: Rscript scripts/lint_check.R

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
  "docs"
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
  # Print all linting issues
  print(lint_results)

  # Count issues by type
  issue_count <- length(lint_results)
  message("\n")
  message(paste("Found", issue_count, "linting issue(s)."))
  message("Run 'make lint' to see details.")

  quit(status = 1)
} else {
  message("No linting issues found!")
  quit(status = 0)
}
