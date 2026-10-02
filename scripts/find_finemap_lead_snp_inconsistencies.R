source("../pipeline_steps/constants.R")
source("../pipeline_steps/gwas_calculations.R")

parser <- argparser::arg_parser(
  "Find finemapped lead SNPs affected by the LBF_P underflow tie-break bug"
)
parser <- argparser::add_argument(
  parser,
  "--study_dir",
  help = paste(
    "Optional: standalone GWAS upload copy directory containing compiled_extracted_studies.tsv.",
    "Scans that file instead of the per-block finemapped_studies.tsv trees"
  ),
  type = "character",
  default = NA
)
parser <- argparser::add_argument(
  parser,
  "--min_p_filter",
  help = paste(
    "Only re-verify rows whose stored min_p is below this value. Lead SNP choice can only",
    "differ from the corrected logic when multiple p-values underflowed to the same double,",
    "which requires min_p in the subnormal/zero range"
  ),
  type = "numeric",
  default = 1e-300
)
parser <- argparser::add_argument(
  parser,
  "--diff_file",
  help = "Output TSV of all re-verified rows with old/new lead SNP, bp, min_p, and cis_trans",
  type = "character",
  default = "finemap_lead_snp_inconsistencies.tsv"
)
parser <- argparser::add_argument(
  parser,
  "--affected_blocks_file",
  help = "Output TSV of affected study/blocks (study_name, ld_block, variant_type, n_changed_rows)",
  type = "character",
  default = "finemap_lead_snp_affected_blocks.tsv"
)
parser <- argparser::add_argument(
  parser,
  "--mc_cores",
  help = "Cores for parallel re-verification of finemapped files",
  type = "numeric",
  default = 8
)
args <- argparser::parse_args(parser)

main <- function() {
  if (!is.na(args$study_dir)) {
    candidates <- load_candidates_from_compiled(args$study_dir)
  } else {
    candidates <- load_candidates_from_ld_blocks()
  }

  message(glue::glue("Re-verifying {nrow(candidates)} candidate row(s) with {args$mc_cores} core(s)"))
  recomputed <- parallel::mclapply(seq_len(nrow(candidates)), function(i) {
    data.table::setDTthreads(1)
    return(recompute_lead_snp(candidates$file[i]))
  }, mc.cores = args$mc_cores)

  candidates$new_snp <- vapply(recomputed, function(r) r$snp, character(1))
  candidates$new_bp <- vapply(recomputed, function(r) r$bp, numeric(1))
  candidates$new_min_p <- vapply(recomputed, function(r) r$min_p, numeric(1))
  candidates$file_readable <- vapply(recomputed, function(r) r$readable, logical(1))

  candidates <- candidates |>
    dplyr::mutate(
      lead_snp_changed = file_readable & !is.na(new_snp) & (is.na(old_snp) | old_snp != new_snp),
      new_cis_trans = dplyr::if_else(
        is.na(orig_cis_trans) | orig_cis_trans != cis_trans$cis_only | is.na(new_bp),
        old_cis_trans,
        dplyr::if_else(abs(orig_bp - new_bp) < 1000000, cis_trans$cis_only, cis_trans$trans_only)
      )
    ) |>
    dplyr::select(
      source, study, unique_study_id, ld_block, variant_type, file,
      old_snp, old_bp, old_min_p, old_cis_trans,
      new_snp, new_bp, new_min_p, new_cis_trans,
      orig_bp, lead_snp_changed, file_readable
    )

  vroom::vroom_write(candidates, args$diff_file)

  affected_blocks <- candidates |>
    dplyr::filter(lead_snp_changed) |>
    dplyr::group_by(study, ld_block, variant_type) |>
    dplyr::summarise(n_changed_rows = dplyr::n(), .groups = "drop") |>
    dplyr::rename(study_name = study) |>
    dplyr::arrange(study_name, ld_block)
  vroom::vroom_write(affected_blocks, args$affected_blocks_file)

  message(glue::glue(
    "Re-verified {nrow(candidates)} row(s): {sum(candidates$lead_snp_changed)} lead SNP(s) changed ",
    "across {nrow(affected_blocks)} study/block(s), {dplyr::n_distinct(affected_blocks$study_name)} stud",
    "y(ies)"
  ))
  message(glue::glue("Unreadable finemapped file(s): {sum(!candidates$file_readable)}"))
  message(glue::glue("Wrote re-verified rows to {args$diff_file}"))
  message(glue::glue("Wrote affected study/blocks to {args$affected_blocks_file}"))
  message(paste(
    "Note: lead SNPs of unfinemapped rows (finemap_message != \"success\") cannot be re-verified",
    "from the finemapped files, which do not store the original P column"
  ))
  return(invisible(candidates))
}

load_candidates_from_ld_blocks <- function() {
  ld_blocks <- vroom::vroom("../pipeline_steps/data/ld_blocks.tsv", show_col_types = FALSE)
  ld_info <- construct_ld_block(ld_blocks$ancestry, ld_blocks$chr, ld_blocks$start, ld_blocks$stop) |>
    dplyr::filter(dir.exists(ld_block_data))

  pipeline_files <- file.path(ld_info$ld_block_data, "finemapped_studies.tsv")
  pipeline_files <- pipeline_files[file.exists(pipeline_files)]
  metadata_files <- data.frame(file = pipeline_files, source = "pipeline")

  upload_root <- glue::glue("{data_dir}ld_blocks/gwas_upload")
  if (dir.exists(upload_root)) {
    upload_metadata <- lapply(list.dirs(upload_root, recursive = FALSE), function(guid) {
      files <- list.files(guid, pattern = "^finemapped_studies\\.tsv$", recursive = TRUE, full.names = TRUE)
      if (length(files) == 0) {
        return(data.frame(file = character(0), source = character(0)))
      }
      return(data.frame(file = files, source = paste0("upload:", basename(guid))))
    })
    metadata_files <- dplyr::bind_rows(metadata_files, upload_metadata)
  }

  message(glue::glue("Scanning {nrow(metadata_files)} finemapped_studies.tsv file(s)"))
  all_candidates <- lapply(seq_len(nrow(metadata_files)), function(i) {
    return(load_candidates_from_metadata_file(metadata_files$file[i], metadata_files$source[i]))
  })
  candidates <- dplyr::bind_rows(all_candidates)
  return(candidates)
}

load_candidates_from_metadata_file <- function(metadata_file, source) {
  required_columns <- c(
    "study", "unique_study_id", "ld_block", "variant_type", "file", "snp", "bp",
    "min_p", "cis_trans", "finemap_message"
  )
  fm <- tryCatch(
    vroom::vroom(metadata_file, col_types = finemapped_column_types, show_col_types = FALSE),
    error = function(e) NULL
  )
  if (is.null(fm) || !all(required_columns %in% names(fm))) {
    message(glue::glue("Skipping unreadable or old-format file: {metadata_file}"))
    return(data.frame())
  }

  candidates <- fm |>
    dplyr::filter(
      finemap_message == "success",
      is.na(min_p) | min_p < args$min_p_filter
    ) |>
    dplyr::transmute(
      study, unique_study_id, ld_block, variant_type,
      file, old_snp = snp, old_bp = bp, old_min_p = min_p, old_cis_trans = cis_trans
    )
  if (nrow(candidates) == 0) {
    return(candidates)
  }

  candidates$source <- source
  candidates$file <- dplyr::if_else(
    grepl("^gwas_upload|^study", candidates$file),
    file.path(data_dir, candidates$file),
    candidates$file
  )

  imputed_file <- file.path(dirname(metadata_file), "imputed_studies.tsv")
  if (file.exists(imputed_file)) {
    imputed <- vroom::vroom(imputed_file, show_col_types = FALSE) |>
      dplyr::transmute(study, orig_bp = bp, orig_cis_trans = cis_trans)
    candidates <- dplyr::left_join(candidates, imputed, by = "study")
  } else {
    candidates$orig_bp <- NA_real_
    candidates$orig_cis_trans <- NA_character_
  }

  return(candidates)
}

load_candidates_from_compiled <- function(study_dir) {
  compiled_file <- file.path(study_dir, "compiled_extracted_studies.tsv")
  if (!file.exists(compiled_file)) {
    stop(glue::glue("File not found: {compiled_file}"))
  }

  compiled <- vroom::vroom(compiled_file, show_col_types = FALSE)
  required_columns <- c("study", "unique_study_id", "ld_block", "file", "snp", "bp", "min_p")
  if (!all(required_columns %in% names(compiled))) {
    stop(glue::glue("Missing columns in {compiled_file}: ",
      paste(setdiff(required_columns, names(compiled)), collapse = ", ")
    ))
  }

  candidates <- compiled |>
    dplyr::filter(is.na(min_p) | min_p < args$min_p_filter) |>
    dplyr::transmute(
      source = paste0("study_dir:", normalizePath(study_dir)),
      study, unique_study_id, ld_block,
      variant_type = NA_character_,
      file = dplyr::if_else(
        grepl("^gwas_upload/", file),
        file.path(study_dir, sub("^gwas_upload/[^/]+/", "", file)),
        file
      ),
      old_snp = snp, old_bp = bp, old_min_p = min_p,
      old_cis_trans = NA_character_, orig_bp = NA_real_, orig_cis_trans = NA_character_
    )
  return(candidates)
}

recompute_lead_snp <- function(file) {
  result <- list(snp = NA_character_, bp = NA_real_, min_p = NA_real_, readable = FALSE)
  if (!file.exists(file)) {
    return(result)
  }

  gwas <- tryCatch(
    data.table::fread(
      file,
      sep = "\t",
      select = c("SNP", "BP", "LBF", "SE"),
      colClasses = c(SNP = "character", BP = "numeric", LBF = "numeric", SE = "numeric"),
      showProgress = FALSE,
      nThread = 1,
      verbose = FALSE
    ),
    error = function(e) NULL
  )
  if (is.null(gwas) || nrow(gwas) == 0) {
    return(result)
  }

  log_p <- convert_lbf_to_log_p_value(gwas$LBF, gwas$SE)
  important_row <- which.min(log_p)
  if (length(important_row) != 1 || is.na(important_row) || important_row < 1) {
    return(result)
  }

  result$snp <- gwas$SNP[important_row]
  result$bp <- gwas$BP[important_row]
  result$min_p <- exp(min(log_p, na.rm = TRUE))
  result$readable <- TRUE
  return(result)
}

invisible(main())
