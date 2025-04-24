source('pipeline_steps/constants.R')

parser <- argparser::arg_parser('Post pipeline cleanup')
#INPUT
parser <- argparser::add_argument(parser, '--studies_processed', help = 'Current state of processed studies', type = 'character')
parser <- argparser::add_argument(parser, '--pipeline_summary_file', help = 'Rendered Rmd file of output', type = 'character')
args <- argparser::parse_args(parser)

main <- function() {
  ld_blocks <- vroom::vroom('pipeline_steps/data/ld_blocks.tsv', show_col_types = F)
  ld_info <- construct_ld_block(ld_blocks$ancestry, ld_blocks$chr, ld_blocks$start, ld_blocks$stop) |>
    dplyr::filter(dir.exists(ld_block_data))

  cleanup_studies_with_no_extractions()

  #only copy the studies_processed.tsv.gz file to the results directory once everything else was successful
  file.copy(args$studies_processed, results_dir)

  # if (is.na(TEST_RUN)) {
  #   rmarkdown::render("pipeline_steps/pipeline_summary.Rmd", output_file = args$pipeline_summary_file)
  # } else {
  #   vroom::vroom_write(data.frame(), args$pipeline_summary_file)
  # }
}

cleanup_studies_with_no_extractions <- function() {
  study_dirs  <- Sys.glob(glue::glue('{extracted_study_dir}/*'))
  empty_study_dirs <- Filter(function(e) file.size(glue::glue('{e}/extracted_snps.tsv')) == 0, study_dirs)
  message('Studies with no extractions that will be cleaned up: ', length(empty_study_dirs))
  for (empty_study in empty_study_dirs) {
    system(glue::glue('rm -r {empty_study}'))
  }
}

main()
