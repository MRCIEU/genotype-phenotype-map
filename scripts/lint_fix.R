#!/usr/bin/env Rscript
# Auto-fix linting issues using styler
# Usage: Rscript scripts/lint_fix.R
#
# Step 1: styler — tidyverse formatting (spacing, indentation, braces)
# Step 2: Trailing whitespace / tab cleanup
# Step 3: Report remaining lint issues

library(styler)
library(lintr)

exclude_dirs <- c(
  "docs",
  "data_formatting",
  "pre_steps/vep_annotation"
)

# ── Step 1: styler ────────────────────────────────────────────────────────────
message("\nStep 1: Formatting with styler (tidyverse style)...")
styler::style_dir(
  path = ".",
  filetype = "r",
  recursive = TRUE,
  exclude_dirs = exclude_dirs,
  exclude_files = c(
    "scripts/lint.R",
    "scripts/lint_check.R",
    "scripts/lint_fix.R",
    "scripts/lint_summary.R"
  ),
  transformers = styler::tidyverse_style(
    indent_by = 2,
    scope = "tokens",
    strict = TRUE
  )
)
message("✓ styler complete!")

# ── Step 2: Trailing whitespace / tab cleanup ─────────────────────────────────
message("\nStep 2: Fixing trailing whitespace and tabs...")

r_files <- list.files(
  path = ".",
  pattern = "\\.[rR]$",
  recursive = TRUE,
  full.names = TRUE
)

r_files <- r_files[!grepl(paste(exclude_dirs, collapse = "|"), r_files)]
r_files <- r_files[!grepl(
  "scripts/lint\\.R$|scripts/lint_check\\.R$|scripts/lint_fix\\.R$|scripts/lint_summary\\.R$",
  r_files
)]

additional_fixes <- 0
for (file in r_files) {
  tryCatch(
    {
      content <- readLines(file, warn = FALSE)
      original_content <- content

      content <- gsub("[ \t]+$", "", content, perl = TRUE)

      content <- sapply(content, function(line) {
        leading_tabs <- regmatches(line, regexpr("^\\t*", line))
        if (nchar(leading_tabs) > 0) {
          spaces <- paste(rep(" ", nchar(leading_tabs) * 2), collapse = "")
          gsub("^\\t+", spaces, line, perl = TRUE)
        } else {
          line
        }
      }, USE.NAMES = FALSE)

      if (!identical(content, original_content)) {
        writeLines(content, file)
        additional_fixes <- additional_fixes + 1
      }
    },
    error = function(e) {
      message(paste("  Warning: Could not process", file, "-", e$message))
    }
  )
}

if (additional_fixes > 0) {
  message(paste("✓ Fixed whitespace in", additional_fixes, "files"))
} else {
  message("✓ No whitespace fixes needed")
}

# ── Step 3: Report remaining lint issues ──────────────────────────────────────
message("\nStep 3: Checking remaining lint issues...")

lint_results <- list()
for (file in r_files) {
  file_lints <- lintr::lint(file)
  if (length(file_lints) > 0) {
    lint_results <- c(lint_results, file_lints)
  }
}

if (length(lint_results) > 0) {
  issue_types <- sapply(lint_results, function(x) x$linter)
  issue_counts <- table(issue_types)

  message("\nRemaining issues by type:")
  print(sort(issue_counts, decreasing = TRUE))

  message(paste("\nTotal remaining issues:", length(lint_results)))
  message("Run 'make lint' to see details.")
} else {
  message("\n✓ No linting issues found!")
}
