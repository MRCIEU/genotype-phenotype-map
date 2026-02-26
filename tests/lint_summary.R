#!/usr/bin/env Rscript
# Analyze and summarize linting issues
# Usage: Rscript scripts/lint_summary.R

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
  message("No R files found.")
  quit(status = 0)
}

message("Analyzing linting issues...")

# Run lintr on all files
lint_results <- list()

for (file in r_files) {
  file_lints <- lintr::lint(file)
  if (length(file_lints) > 0) {
    lint_results <- c(lint_results, file_lints)
  }
}

if (length(lint_results) == 0) {
  message("✓ No linting issues found!")
  quit(status = 0)
}

# Categorize issues
issue_types <- sapply(lint_results, function(x) class(x$linter)[1])
issue_messages <- sapply(lint_results, function(x) x$message)

# Count by type
issue_counts <- table(issue_types)
issue_counts_sorted <- sort(issue_counts, decreasing = TRUE)

message("\n=== Linting Issues Summary ===\n")
message(paste("Total issues:", length(lint_results)))
message(paste("Files checked:", length(r_files)))
message("\nIssues by type:")
print(issue_counts_sorted)

# Show most common messages
message("\n=== Most Common Issue Messages ===\n")
common_messages <- head(sort(table(issue_messages), decreasing = TRUE), 10)
print(common_messages)

# Files with most issues
file_issue_counts <- table(sapply(lint_results, function(x) x$filename))
file_issue_counts_sorted <- head(sort(file_issue_counts, decreasing = TRUE), 10)

message("\n=== Files with Most Issues ===\n")
print(file_issue_counts_sorted)

message("\n=== Tips ===\n")
message("1. Run 'make lint-fix' to auto-fix formatting issues")
message("2. Common manual fixes needed:")
message("   - Change '=' to '<-' for assignments")
message("   - Break long lines (120 char limit)")
message("   - Remove trailing whitespace")
message("   - Fix spacing around operators")
