#!/usr/bin/env Rscript
# Auto-fix linting issues using styler
# Usage: Rscript scripts/lint_format.R
#
# Step 1: styler — tidyverse formatting (spacing, indentation, braces)
# Step 2: Trailing whitespace / tab cleanup
# Step 3: Report remaining lint issues

library(styler)
library(lintr)

exclude_dirs <- c(
  "docs"
)

# ── Step 1: styler ────────────────────────────────────────────────────────────
message("\nStep 1: Formatting with styler (tidyverse style)...")
styler::style_dir(
  path = ".",
  filetype = "r",
  recursive = TRUE,
  exclude_dirs = exclude_dirs,
  transformers = styler::tidyverse_style(
    indent_by = 2,
    scope = "spaces",
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
          return(gsub("^\\t+", spaces, line, perl = TRUE))
        } else {
          return(line)
        }
      }, USE.NAMES = FALSE)

      if (!identical(content, original_content)) {
        writeLines(content, file)
        additional_fixes <- additional_fixes + 1
      }
    },
    error = function(e) {
      message(paste("  Warning: Could not process", file, "-", e$message))
      return()
    }
  )
}

if (additional_fixes > 0) {
  message(paste("✓ Fixed whitespace in", additional_fixes, "files"))
} else {
  message("✓ No whitespace fixes needed")
}
