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
  "--worker_guid",
  help = "Worker GUID (if invoked by worker)",
  type = "character",
  default = NA
)
args <- argparser::parse_args(parser)

main <- function() {
  if (!is.na(args$worker_guid)) {
    update_directories_for_worker(args$worker_guid)
  }
  ld_info <- ld_block_dirs(args$ld_block)
  ld_matrix_info <- vroom::vroom(glue::glue("{ld_info$ld_reference_panel_prefix}.tsv"), show_col_types = F)

  extracted_studies_file <- glue::glue("{ld_info$ld_block_data}/extracted_studies.tsv")
  extracted_studies <- vroom::vroom(extracted_studies_file, show_col_types = F)

  standardised_studies_file <- glue::glue("{ld_info$ld_block_data}/standardised_studies.tsv")
  if (file.exists(standardised_studies_file)) {
    existing_standardised_studies <- vroom::vroom(standardised_studies_file,
      show_col_types = F,
      col_types = standardised_column_types
    )
  } else {
    existing_standardised_studies <- empty_standardised_studies()
  }

  if (nrow(extracted_studies) > 0) {
    standardised_studies <- apply(extracted_studies, 1, function(study) {
      start_time <- Sys.time()
      standardised_file <- sub("extracted", "standardised", study[["file"]])

      if (standardised_file %in% existing_standardised_studies$file) {
        return()
      }

      result <- perform_standardisation(study, ld_matrix_info)

      if (nrow(result$gwas) < minimum_extraction_size &&
          study[["variant_type"]] == variant_types$common &&
          study[["coverage"]] == coverage_types$dense
      ) {
        return()
      }
      # TODO: Remove this once we want to ingest UK Biobank PPP trans-only studies
      if (!is.na(study[["cis_trans"]]) && study[["cis_trans"]] == cis_trans$trans_only) {
        return()
      }
      vroom::vroom_write(result$gwas, result$study$file)

      result$study$time_taken <- hms::as_hms(difftime(Sys.time(), start_time))
      return(result$study)
    }) |>
      dplyr::bind_rows() |>
      type.convert(as.is = T)

    if (nrow(standardised_studies) > 0) {
      standardised_studies$chr <- as.character(standardised_studies$chr)
    }
  }

  if (nrow(standardised_studies) > 0) {
    standardised_studies <- dplyr::bind_rows(existing_standardised_studies, standardised_studies) |>
      dplyr::distinct()

    vroom::vroom_write(standardised_studies, standardised_studies_file)
  } else {
    vroom::vroom_write(empty_standardised_studies(), standardised_studies_file)
  }

  vroom::vroom_write(data.frame(), args$completed_output_file)
  return()
}

perform_standardisation <- function(study, ld_matrix_info) {
  standardised_file <- sub("extracted", "standardised", study[["file"]])
  gwas <- vroom::vroom(study[["file"]], show_col_types = F, col_types = vroom::cols(
    EA = vroom::col_character(),
    OA = vroom::col_character()
  ))
  is_rare_study <- study[["variant_type"]] != variant_types$common

  response <- standardise_alleles(gwas) |>
    standardise_extracted_gwas(ld_matrix_info, is_rare_study)

  response$gwas <- gwas_health_check(response$gwas) |>
    filter_gwas(is_rare_study)

  study["ld_block"] <- args$ld_block
  study["file"] <- standardised_file
  study["eaf_from_reference_panel"] <- response$eaf_from_reference_panel
  study["snps_removed_by_reference_panel"] <- response$snps_removed_by_reference_panel
  study <- study[-match("reference_build", names(study))]

  return(list(gwas = response$gwas, study = as.list(study)))
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

main()
