source("../pipeline_steps/constants.R")

parser <- argparser::arg_parser(paste(
  "Repair finemapped_studies.tsv lead SNPs affected by the LBF_P underflow tie-break bug,",
  "using the diff produced by find_finemap_lead_snp_inconsistencies.R.",
  "Only rewrites snp, bp, min_p, and cis_trans metadata; the finemapped GWAS files",
  "themselves are never modified"
))
parser <- argparser::add_argument(
  parser,
  "--diff_file",
  help = "Diff TSV produced by find_finemap_lead_snp_inconsistencies.R",
  type = "character",
  default = "finemap_lead_snp_inconsistencies.tsv"
)
parser <- argparser::add_argument(
  parser,
  "--apply",
  help = "Apply the repairs (pass --apply TRUE). Without this the script only reports intended changes",
  type = "logical",
  default = FALSE
)
args <- argparser::parse_args(parser)

main <- function() {
  diff <- vroom::vroom(args$diff_file, show_col_types = FALSE)
  required_columns <- c(
    "source", "study", "unique_study_id", "ld_block", "old_snp", "old_bp",
    "old_min_p", "old_cis_trans", "new_snp", "new_bp", "new_min_p",
    "new_cis_trans", "lead_snp_changed"
  )
  if (!all(required_columns %in% names(diff))) {
    stop(glue::glue("Missing columns in {args$diff_file}: ",
      paste(setdiff(required_columns, names(diff)), collapse = ", ")
    ))
  }

  changed <- diff |>
    dplyr::filter(lead_snp_changed)

  skipped <- changed |>
    dplyr::filter(!grepl("^pipeline$|^upload:", source))
  if (nrow(skipped) > 0) {
    message(glue::glue(
      "Skipping {nrow(skipped)} row(s) from standalone study_dir copies; ",
      "repair the canonical per-block finemapped_studies.tsv trees instead"
    ))
    message(paste(unique(skipped$source), collapse = ", "))
  }

  changed <- changed |>
    dplyr::filter(grepl("^pipeline$|^upload:", source)) |>
    dplyr::mutate(metadata_file = dplyr::case_when(
      source == "pipeline" ~ file.path(data_dir, "ld_blocks", ld_block, "finemapped_studies.tsv"),
      grepl("^upload:", source) ~ file.path(
        data_dir, "ld_blocks", "gwas_upload", sub("^upload:", "", source), ld_block,
        "finemapped_studies.tsv"
      ),
      TRUE ~ NA_character_
    )) |>
    dplyr::filter(!is.na(metadata_file) & file.exists(metadata_file))

  metadata_files <- unique(changed$metadata_file)
  message(glue::glue(
    "{if (args$apply) 'Repairing' else 'Dry run: would repair'} ",
    "{nrow(changed)} row(s) across {length(metadata_files)} finemapped_studies.tsv file(s)"
  ))

  n_repaired <- 0
  n_skipped <- 0
  for (metadata_file in metadata_files) {
    rows <- dplyr::filter(changed, metadata_file == !!metadata_file)
    fm <- vroom::vroom(metadata_file, col_types = finemapped_column_types, show_col_types = FALSE)

    repaired_rows <- 0
    skipped_rows <- character(0)
    for (i in seq_len(nrow(rows))) {
      row <- rows[i, ]
      m <- match(row$unique_study_id, fm$unique_study_id)
      if (length(m) != 1 || is.na(m)) {
        skipped_rows <- c(skipped_rows, glue::glue("{row$unique_study_id}: not found"))
        next
      }
      metadata_matches_scan <- as.character(fm$snp[m]) == as.character(row$old_snp) &&
        as.numeric(fm$bp[m]) == as.numeric(row$old_bp) &&
        as.numeric(fm$min_p[m]) == as.numeric(row$old_min_p)
      if (!metadata_matches_scan) {
        skipped_rows <- c(
          skipped_rows,
          glue::glue("{row$unique_study_id}: metadata changed since scan, rerun the scan")
        )
        next
      }

      fm$snp[m] <- as.character(row$new_snp)
      fm$bp[m] <- as.numeric(row$new_bp)
      fm$min_p[m] <- as.numeric(row$new_min_p)
      if ("cis_trans" %in% names(fm) && !is.na(row$new_cis_trans)) {
        fm$cis_trans[m] <- as.character(row$new_cis_trans)
      }
      repaired_rows <- repaired_rows + 1
    }

    if (args$apply) {
      vroom::vroom_write(fm, metadata_file)
    }
    n_repaired <- n_repaired + repaired_rows
    n_skipped <- n_skipped + length(skipped_rows)
    message(glue::glue(
      "{metadata_file}: {repaired_rows} row(s) {if (args$apply) 'repaired' else 'to repair'}",
      ", {length(skipped_rows)} skipped"
    ))
    for (skip_reason in skipped_rows) {
      message(glue::glue("  skipped: {skip_reason}"))
    }
  }

  message(glue::glue(
    "{n_repaired} row(s) {if (args$apply) 'repaired' else 'to repair'}, {n_skipped} skipped"
  ))
  if (!args$apply) {
    message("Dry run only; rerun with --apply TRUE to write the repairs")
  }
  message(paste(
    "After applying, delete coloc_complete/clustering_complete sentinels and the coloc/",
    "clustering outputs of the affected blocks, then rerun coloc, clustering, compile,",
    "and the DB/static web steps"
  ))
  return(invisible(changed))
}

invisible(main())
