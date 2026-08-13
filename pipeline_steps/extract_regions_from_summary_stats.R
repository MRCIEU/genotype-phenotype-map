source("constants.R")
source("common_extraction_functions.R")

parser <- argparser::arg_parser("Extract regions from a delimited file (csv, tsv)")
parser <- argparser::add_argument(
  parser,
  "--worker_guid",
  help = "Worker GUID (worker pipeline path)",
  type = "character",
  default = NA
)
parser <- argparser::add_argument(
  parser,
  "--extracted_study_location",
  help = "Study extraction directory (main pipeline path)",
  type = "character",
  default = NA
)
parser <- argparser::add_argument(
  parser,
  "--extracted_output_file",
  help = "Extracted output file (main pipeline path)",
  type = "character",
  default = NA
)

args <- argparser::parse_args(parser)
is_worker <- !is.na(args$worker_guid)

main <- function() {
  if (is_worker) {
    update_directories_for_worker(args$worker_guid)
    study_metadata <- jsonlite::fromJSON(glue::glue("{extracted_study_dir}/study_metadata.json"))
    if (is.null(study_metadata$study_name) || is.na(study_metadata$study_name)) {
      study_metadata$study_name <- ifelse(
        !is.null(study_metadata$name) && !is.na(study_metadata$name),
        study_metadata$name,
        study_metadata$guid
      )
    }
    extracted_output_file <- glue::glue("{extracted_study_dir}/extracted_snps.tsv")
  } else {
    if (is.na(args$extracted_study_location) || is.na(args$extracted_output_file)) {
      stop(paste(
        "Error: either --worker_guid or both --extracted_study_location and",
        "--extracted_output_file are required"
      ))
    }

    study <- vroom::vroom(pipeline_metadata_file_paths()$studies_to_process, show_col_types = F) |>
      dplyr::filter(extracted_location == args$extracted_study_location)
    if (nrow(study) != 1) stop("Error: cant find study to process")

    extracted_study_dir <<- args$extracted_study_location
    study_metadata <- build_study_metadata_from_studies_to_process(study)
    extracted_output_file <- args$extracted_output_file
  }

  dir.create(glue::glue("{extracted_study_dir}/gwas"), showWarnings = F, recursive = T)
  dir.create(glue::glue("{extracted_study_dir}/extracted"), showWarnings = F, recursive = T)
  dir.create(glue::glue("{extracted_study_dir}/standardised"), showWarnings = F, recursive = T)
  dir.create(glue::glue("{extracted_study_dir}/imputed"), showWarnings = F, recursive = T)
  dir.create(glue::glue("{extracted_study_dir}/finemapped"), showWarnings = F, recursive = T)
  # Full-trait SVGs are main-pipeline only (worker uses finemap/upload UI paths instead)
  if (!is_worker) {
    dir.create(glue::glue("{extracted_study_dir}/svgs"), showWarnings = F, recursive = T)
    study_metadata$extracted_location <- extracted_study_dir
  }

  p_value_threshold <- ifelse(
    is.na(study_metadata$p_value_threshold) || is.null(study_metadata$p_value_threshold),
    lowest_p_value_threshold,
    as.numeric(study_metadata$p_value_threshold)
  )

  ld_blocks <- vroom::vroom("data/ld_blocks.tsv", show_col_types = F) |>
    dplyr::filter(ancestry == study_metadata$ancestry)

  extract_summary_stats_regions(
    study_metadata,
    ld_blocks,
    p_value_threshold,
    extracted_output_file,
    create_full_svgs = !is_worker
  )
  return(invisible(NULL))
}

build_study_metadata_from_studies_to_process <- function(study) {
  column_names <- list()
  if (!is.null(study$column_names) && !is.na(study$column_names) && nchar(study$column_names) > 0) {
    column_names <- jsonlite::fromJSON(study$column_names)
  }

  file_type <- study$file_type
  if (is.null(file_type) || is.na(file_type) || file_type == "") {
    file_type <- extraction_file_types$csv
  }

  return(list(
    study_name = study$study_name,
    file_location = study$study_location,
    ancestry = study$ancestry,
    sample_size = study$sample_size,
    category = study$category,
    file_type = file_type,
    reference_build = study$reference_build,
    p_value_threshold = study$p_value_threshold,
    column_names = column_names
  ))
}

extract_summary_stats_regions <- function(
  study_metadata,
  ld_blocks,
  p_value_threshold,
  extracted_output_file,
  create_full_svgs = FALSE
) {
  if (study_metadata$file_type == extraction_file_types$vcf) {
    stop("VCF support not yet implemented")
  } else if (study_metadata$file_type == extraction_file_types$csv) {
    gwas <- vroom::vroom(study_metadata$file_location, show_col_types = F)
    study_metadata$column_names <- Filter(
      \(column) !is.null(column) && !is.na(column),
      study_metadata$column_names
    )
    gwas <- change_column_names(gwas, study_metadata$column_names)
    gwas <- standardise_columns(gwas)

    start_time <- Sys.time()
    if (study_metadata$reference_build != reference_builds$GRCh38) {
      gwas <- convert_dataframe_reference_build(gwas, study_metadata$reference_build)
    }
    time_taken <- diff_time_taken(start_time)
    message(glue::glue("Time taken to convert reference build: {time_taken}"))

    if (create_full_svgs) {
      source("svg_helpers.R")
      gwas_for_svg <- dplyr::mutate(gwas, LP = -log10(P))
      create_svgs_from_gwas(study_metadata, gwas_for_svg)
    }

    extracted_regions <- split_into_regions(gwas, ld_blocks, study_metadata, p_value_threshold)
    vroom::vroom_write(extracted_regions, extracted_output_file)
    return(invisible(NULL))
  } else {
    stop(paste(c(
      "Error: file type", study_metadata$file_type, "not recocognised.",
      "File types must be one of:", extraction_file_types
    ), collapse = " "))
  }
}

invisible(main())
