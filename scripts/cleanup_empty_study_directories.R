source("../pipeline_steps/constants.R")
source("../pipeline_steps/study_directory_helpers.R")

parser <- argparser::arg_parser("Remove study directories with no extracted regions")
parser <- argparser::add_argument(
  parser,
  "--dry_run",
  help = "List directories that would be removed without deleting them",
  type = "logical",
  default = FALSE
)
args <- argparser::parse_args(parser)

main <- function() {
  cleanup_empty_study_directories(dry_run = args$dry_run)
  return(invisible(NULL))
}

study_has_no_extractions <- function(study_dir) {
  extracted_snps_file <- glue::glue("{study_dir}/extracted_snps.tsv")
  if (!file.exists(extracted_snps_file)) {
    return(FALSE)
  }
  if (file.info(extracted_snps_file)$size == 0) {
    return(TRUE)
  }
  extracted_snps <- data.table::fread(
    extracted_snps_file,
    nThread = 1,
    verbose = FALSE,
    showProgress = FALSE
  )
  return(nrow(extracted_snps) == 0)
}

cleanup_empty_study_directories <- function(dry_run = FALSE) {
  study_dirs <- Sys.glob(glue::glue("{extracted_study_dir}*"))
  study_dirs <- study_dirs[file.info(study_dirs)$isdir]
  empty_study_dirs <- Filter(study_has_no_extractions, study_dirs)

  message("Studies with no extractions that will be cleaned up: ", length(empty_study_dirs))
  if (length(empty_study_dirs) == 0) {
    return(invisible(empty_study_dirs))
  }

  for (study_dir in empty_study_dirs) {
    if (dry_run) {
      message("Would remove: ", study_dir)
    } else {
      message("Removing: ", study_dir)
      unlink(study_dir, recursive = TRUE)
    }
  }

  return(invisible(empty_study_dirs))
}


invisible(main())
