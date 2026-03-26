source("constants.R")
source("common_extraction_functions.R")

parser <- argparser::arg_parser("Extract regions from a delimited file (csv, tsv)")
parser <- argparser::add_argument(
  parser,
  "--worker_guid",
  help = "Worker GUID",
  type = "character",
  default = NA
)

args <- argparser::parse_args(parser)

main <- function() {
  if (is.na(args$worker_guid)) stop("Error: worker_guid is required for summary stats extraction")
  update_directories_for_worker(args$worker_guid)

  dir.create(glue::glue("{extracted_study_dir}/gwas"), showWarnings = F, recursive = T)
  dir.create(glue::glue("{extracted_study_dir}/extracted"), showWarnings = F, recursive = T)
  dir.create(glue::glue("{extracted_study_dir}/standardised"), showWarnings = F, recursive = T)
  dir.create(glue::glue("{extracted_study_dir}/imputed"), showWarnings = F, recursive = T)
  dir.create(glue::glue("{extracted_study_dir}/finemapped"), showWarnings = F, recursive = T)

  study_metadata <- jsonlite::fromJSON(glue::glue("{extracted_study_dir}/study_metadata.json"))
  p_value_threshold <- ifelse(
    is.na(study_metadata$p_value_threshold) || is.null(study_metadata$p_value_threshold),
    lowest_p_value_threshold,
    study_metadata$p_value_threshold
  )

  ld_blocks <- vroom::vroom("data/ld_blocks.tsv", show_col_types = F) |>
    dplyr::filter(ancestry == study_metadata$ancestry)

  clumped_hits_file <- glue::glue("{extracted_study_dir}/clumped_snps.tsv")

  if (study_metadata$file_type == extraction_file_types$vcf) {
    stop("VCF support not yet implemented")
    # TODO: Add VCF support
    # if (study$reference_build == reference_builds$GRCh37) {
    # vcf_file <- convert_reference_build(study, extraction_file_types$vcf, vcf_file)
    # }
  } else if (study_metadata$file_type == extraction_file_types$csv) {
    gwas <- vroom::vroom(study_metadata$file_location, show_col_types = F)
    study_metadata$column_names <- Filter(\(column) !is.null(column) && !is.na(column), study_metadata$column_names)
    gwas <- change_column_names(gwas, study_metadata$column_names)
    gwas <- standardise_columns(gwas)

    start_time <- Sys.time()
    if (study_metadata$reference_build != reference_builds$GRCh38) {
      gwas <- convert_dataframe_reference_build(gwas, study_metadata$reference_build)
    }
    time_taken <- diff_time_taken(start_time)
    message(glue::glue("Time taken to convert reference build: {time_taken}"))

    extracted_regions <- split_into_regions(gwas, ld_blocks, study_metadata, p_value_threshold)
    extracted_output_file <- glue::glue("{extracted_study_dir}/extracted_snps.tsv")
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
