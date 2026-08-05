source("constants.R")
source("common_extraction_functions.R")
source("gwas_calculations.R")

parser <- argparser::arg_parser("Standardise GWAS for pipeline")
parser <- argparser::add_argument(
  parser,
  "--ld_block",
  help = "LD block that the ",
  type = "character"
)
parser <- argparser::add_argument(
  parser,
  "--completed_output_file",
  help = "Completed output file",
  type = "character"
)
parser <- argparser::add_argument(
  parser,
  "--block_list",
  help = "CSV of studies to exclude from standardisation (columns: id_pattern, cis_trans)",
  type = "character",
  default = NA
)
parser <- argparser::add_argument(
  parser,
  "--worker_guid",
  help = "Worker GUID (if invoked by worker)",
  type = "character",
  default = NA
)
args <- argparser::parse_args(parser)

standardisation_skip_rule_version <- glue::glue(
  "min_dense={minimum_extraction_size_for_dense_coverage},min_sparse={minimum_extraction_size_for_sparse_coverage}"
)

main <- function() {
  if (!is.na(args$worker_guid)) {
    update_directories_for_worker(args$worker_guid)
  }
  paths <- ld_block_file_paths(args$ld_block, block_list = args$block_list)

  if (!file.exists(paths$extracted_studies)) {
    message(glue::glue("{args$ld_block}: No extracted studies; writing empty standardised results."))
    vroom::vroom_write(empty_standardised_studies(), paths$standardised_studies)
    if (!file.exists(paths$standardised_skipped)) {
      vroom::vroom_write(empty_standardised_skipped(), paths$standardised_skipped)
    }
    vroom::vroom_write(data.frame(), args$completed_output_file)
    return()
  }

  extracted_studies <- vroom::vroom(paths$extracted_studies, show_col_types = F)

  block_list <- NULL
  if (!is.null(args$block_list) && !is.na(args$block_list) && file.exists(args$block_list)) {
    block_list <- vroom::vroom(args$block_list, show_col_types = FALSE)
  }

  if (!is.null(block_list) && nrow(block_list) > 0 && nrow(extracted_studies) > 0) {
    blocked <- is_study_blocked(block_list, extracted_studies$study, extracted_studies$cis_trans)
    if (sum(blocked) > 0) {
      extracted_studies <- extracted_studies |> dplyr::filter(!blocked)
      message(glue::glue("{args$ld_block}: Excluded {sum(blocked)} blocked studies from standardisation"))
    }
  }

  standardised_studies_file <- paths$standardised_studies
  if (file.exists(standardised_studies_file)) {
    existing_standardised_studies <- vroom::vroom(standardised_studies_file,
      show_col_types = F,
      col_types = standardised_column_types
    )
  } else {
    existing_standardised_studies <- empty_standardised_studies()
  }

  existing_skipped <- load_standardised_skipped(paths$standardised_skipped)
  active_skipped_inputs <- existing_skipped$input_file[
    existing_skipped$rule_version == standardisation_skip_rule_version
  ]

  expected_standardised_files <- sub("extracted", "standardised", extracted_studies$file)
  already_done <- expected_standardised_files %in% existing_standardised_studies$file |
    extracted_studies$file %in% active_skipped_inputs
  todo_idx <- which(!already_done)

  if (length(todo_idx) == 0L) {
    message(glue::glue(
      "{args$ld_block}: All {nrow(extracted_studies)} studies already standardised or skipped; skipping."
    ))
    if (!file.exists(standardised_studies_file)) {
      vroom::vroom_write(existing_standardised_studies, standardised_studies_file)
    }
    if (!file.exists(paths$standardised_skipped)) {
      vroom::vroom_write(existing_skipped, paths$standardised_skipped)
    }
    vroom::vroom_write(data.frame(), args$completed_output_file)
    return()
  }

  message(glue::glue(
    "{args$ld_block}: Standardising {length(todo_idx)} / {nrow(extracted_studies)} studies"
  ))
  ld_matrix_info <- vroom::vroom(paths$ld_matrix_tsv, show_col_types = F)
  studies_to_standardise <- extracted_studies[todo_idx, , drop = FALSE]

  results <- lapply(seq_len(nrow(studies_to_standardise)), function(i) {
    study <- studies_to_standardise[i, , drop = FALSE]
    start_time <- Sys.time()
    result <- perform_standardisation(study, ld_matrix_info)
    n_variants <- nrow(result$gwas)

    if (n_variants < minimum_extraction_size_for_dense_coverage &&
        study[["variant_type"]] == variant_types$common &&
        study[["coverage"]] == coverage_types$dense
    ) {
      return(list(
        status = "skipped",
        skipped = make_standardised_skipped_row(
          study = study,
          reason = "below_min_dense",
          n_variants = n_variants
        )
      ))
    }

    if (
      n_variants < minimum_extraction_size_for_sparse_coverage &&
        study[["variant_type"]] == variant_types$common &&
        study[["coverage"]] == coverage_types$sparse
    ) {
      return(list(
        status = "skipped",
        skipped = make_standardised_skipped_row(
          study = study,
          reason = "below_min_sparse",
          n_variants = n_variants
        )
      ))
    }

    vroom::vroom_write(result$gwas, result$study$file)

    result$study$time_taken <- hms::as_hms(difftime(Sys.time(), start_time))
    return(list(status = "ok", study = result$study))
  })

  new_standardised_studies <- dplyr::bind_rows(
    lapply(results[vapply(results, function(r) identical(r$status, "ok"), logical(1))], `[[`, "study")
  )
  new_skipped <- dplyr::bind_rows(
    lapply(results[vapply(results, function(r) identical(r$status, "skipped"), logical(1))], `[[`, "skipped")
  )

  if (nrow(new_standardised_studies) > 0) {
    new_standardised_studies <- type.convert(new_standardised_studies, as.is = TRUE)
    new_standardised_studies$chr <- as.character(new_standardised_studies$chr)
  }

  if (nrow(new_standardised_studies) > 0) {
    standardised_studies <- dplyr::bind_rows(existing_standardised_studies, new_standardised_studies) |>
      dplyr::distinct()

    vroom::vroom_write(standardised_studies, standardised_studies_file)
  } else if (nrow(existing_standardised_studies) > 0) {
    vroom::vroom_write(existing_standardised_studies, standardised_studies_file)
  } else {
    vroom::vroom_write(empty_standardised_studies(), standardised_studies_file)
  }

  if (nrow(new_skipped) > 0) {
    skipped_studies <- dplyr::bind_rows(existing_skipped, new_skipped) |>
      dplyr::distinct(input_file, rule_version, .keep_all = TRUE)
    vroom::vroom_write(skipped_studies, paths$standardised_skipped)
  } else if (!file.exists(paths$standardised_skipped)) {
    vroom::vroom_write(existing_skipped, paths$standardised_skipped)
  }

  vroom::vroom_write(data.frame(), args$completed_output_file)
  return()
}

perform_standardisation <- function(study, ld_matrix_info) {
  # study may be a 1-row data.frame (preferred) or a named vector from apply()
  get_field <- function(name) {
    if (is.data.frame(study)) {
      return(study[[name]][[1]])
    }
    return(study[[name]])
  }

  standardised_file <- sub("extracted", "standardised", get_field("file"))
  gwas <- vroom::vroom(get_field("file"), show_col_types = F, col_types = vroom::cols(
    EA = vroom::col_character(),
    OA = vroom::col_character()
  ))
  is_rare_study <- get_field("variant_type") != variant_types$common

  response <- standardise_alleles(gwas) |>
    standardise_extracted_gwas(ld_matrix_info, is_rare_study)

  response$gwas <- gwas_health_check(response$gwas) |>
    filter_gwas(is_rare_study)

  study_out <- if (is.data.frame(study)) {
    as.list(study[1, , drop = TRUE])
  } else {
    as.list(study)
  }
  study_out[["ld_block"]] <- args$ld_block
  study_out[["file"]] <- standardised_file
  study_out[["eaf_from_reference_panel"]] <- response$eaf_from_reference_panel
  study_out[["snps_removed_by_reference_panel"]] <- response$snps_removed_by_reference_panel
  study_out[["reference_build"]] <- NULL

  return(list(gwas = response$gwas, study = study_out))
}

empty_standardised_studies <- function() {
  return(data.frame(
    study = character(),
    file = character(),
    ancestry = character(),
    chr = character(),
    bp = numeric(),
    p_value_threshold = numeric(),
    category = character(),
    sample_size = numeric(),
    cis_trans = character(),
    eaf_from_reference_panel = logical(),
    snps_removed_by_reference_panel = numeric(),
    time_taken = character(),
    variant_type = character(),
    coverage = character()
  ))
}

empty_standardised_skipped <- function() {
  return(data.frame(
    study = character(),
    input_file = character(),
    reason = character(),
    n_variants = integer(),
    coverage = character(),
    variant_type = character(),
    ld_block = character(),
    rule_version = character(),
    stringsAsFactors = FALSE
  ))
}

load_standardised_skipped <- function(skipped_file) {
  if (!file.exists(skipped_file)) {
    return(empty_standardised_skipped())
  }
  skipped <- vroom::vroom(skipped_file, show_col_types = FALSE)
  expected_cols <- names(empty_standardised_skipped())
  missing_cols <- setdiff(expected_cols, names(skipped))
  for (col in missing_cols) {
    skipped[[col]] <- NA
  }
  return(skipped[, expected_cols, drop = FALSE])
}

make_standardised_skipped_row <- function(study, reason, n_variants) {
  get_field <- function(name) {
    if (is.data.frame(study)) {
      return(study[[name]][[1]])
    }
    return(study[[name]])
  }
  return(data.frame(
    study = get_field("study"),
    input_file = get_field("file"),
    reason = reason,
    n_variants = as.integer(n_variants),
    coverage = get_field("coverage"),
    variant_type = get_field("variant_type"),
    ld_block = args$ld_block,
    rule_version = as.character(standardisation_skip_rule_version),
    stringsAsFactors = FALSE
  ))
}

standardise_extracted_gwas <- function(gwas, ld_matrix_info, is_rare_study = F) {
  eaf_from_reference_panel <- FALSE
  original_gwas_size <- nrow(gwas)
  gwas <- dplyr::distinct(gwas, CHR, BP, EA, OA, .keep_all = TRUE)

  if (!"Z" %in% colnames(gwas)) {
    gwas <- dplyr::mutate(gwas, SE = replace(SE, SE == 0, 0.00001), Z = BETA / SE)
  }

  if (!"P" %in% colnames(gwas) && "LP" %in% colnames(gwas)) {
    gwas$LP <- as.numeric(gwas$LP)
    gwas <- dplyr::mutate(gwas, P = 10^(-LP)) |>
      dplyr::select(-LP)
  }

  if (is_rare_study) {
    return(list(
      gwas = gwas,
      eaf_from_reference_panel = FALSE,
      snps_removed_by_reference_panel = 0
    ))
  }
  gwas <- dplyr::filter(gwas, SNP %in% ld_matrix_info$SNP)
  ld_matrix_info <- dplyr::filter(ld_matrix_info, SNP %in% gwas$SNP)

  if (all(is.na(gwas$EAF))) {
    gwas <- dplyr::select(gwas, -EAF) |>
      dplyr::left_join(ld_matrix_info |> dplyr::select(SNP, EAF), by = "SNP")
    eaf_from_reference_panel <- TRUE
  }

  columns_to_coerce <- c("EAF") # Add BETA and SE if needed
  gwas <- tidyr::drop_na(gwas, dplyr::all_of(columns_to_coerce)) |>
    dplyr::arrange(match(SNP, ld_matrix_info$SNP))

  return(list(
    gwas = gwas,
    eaf_from_reference_panel = eaf_from_reference_panel,
    snps_removed_by_reference_panel = original_gwas_size - nrow(gwas)
  ))
}

standardise_alleles <- function(gwas) {
  gwas$EA <- toupper(gwas$EA)
  gwas$OA <- toupper(gwas$OA)
  columns_to_coerce <- c("EAF", "BETA", "SE")
  gwas <- dplyr::mutate(gwas, dplyr::across(dplyr::all_of(columns_to_coerce), as.numeric))

  to_flip <- (gwas$EA > gwas$OA) & (!gwas$EA %in% c("D", "I"))
  gwas <- flip_alleles(gwas, to_flip)

  gwas$SNP <- format_unique_snp_string(gwas$CHR, gwas$BP, gwas$EA, gwas$OA)
  return(gwas)
}

gwas_health_check <- function(gwas) {
  if (any(gwas$P < 0 | gwas$P > 1, na.rm = T)) {
    stop("GWAS has some P values outside accepted range.  Please fix GWAS or remove it from pipeline")
  }
  if (any(as.numeric(gwas$SE) < 0, na.rm = T)) {
    stop("GWAS has some SE values outside accepted range.  Please fix GWAS or remove it from pipeline")
  }
  return(gwas)
}

filter_gwas <- function(gwas, is_rare_study = F) {
  if (is_rare_study) {
    return(gwas)
  }
  gwas <- dplyr::filter(
    gwas,
    (EAF < 0.995 & EAF > 0.005) &
      !is.na(CHR) & !is.na(BP) & !is.na(EA) & !is.na(OA) &
      !is.na(BETA) & !is.na(SE)
  )

  return(gwas)
}

invisible(main())
