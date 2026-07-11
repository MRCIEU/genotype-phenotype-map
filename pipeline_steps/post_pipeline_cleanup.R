source("pipeline_steps/constants.R")
source("pipeline_steps/study_directory_helpers.R")

parser <- argparser::arg_parser("Post pipeline cleanup")
# INPUT
parser <- argparser::add_argument(
  parser,
  "--current_results_dir",
  help = "Current results directory",
  type = "character"
)
parser <- argparser::add_argument(
  parser,
  "--pipeline_summary_file",
  help = "Rendered Rmd file of output",
  type = "character"
)
args <- argparser::parse_args(parser)

main <- function() {
  ld_blocks <- vroom::vroom("pipeline_steps/data/ld_blocks.tsv", show_col_types = F)
  ld_info <- construct_ld_block(ld_blocks$ancestry, ld_blocks$chr, ld_blocks$start, ld_blocks$stop) |>
    dplyr::filter(dir.exists(ld_block_data))

  cleanup_studies_with_no_extractions()

  # only copy the studies_processed.tsv.gz file to the results directory once everything else was successful
  files <- Sys.glob(glue::glue("{args$current_results_dir}/*.tsv.gz"))
  for (file in files) {
    file.copy(file, latest_results_dir, overwrite = T)
  }

  if (is.na(TEST_RUN)) {
    rmarkdown::render("pipeline_steps/pipeline_summary.Rmd",
      output_file = args$pipeline_summary_file,
      params = list(compare_results = F)
    )
  } else {
    vroom::vroom_write(data.frame(), args$pipeline_summary_file)
  }
  return()
}

cleanup_studies_with_no_extractions <- function() {
  cleanup_empty_study_directories(dry_run = FALSE)
  return(invisible(NULL))
}

main()
