source('pipeline_steps/constants.R')

posterior_prob_threshold <- 0.5

parser <- argparser::arg_parser('Compile results from pipeline')
#INPUT
parser <- argparser::add_argument(parser, '--studies_to_process', help = 'Studies to process', type = 'character')
parser <- argparser::add_argument(parser, '--studies_processed', help = 'Current state of processed studies', type = 'character')
#OUTPUT

args <- argparser::parse_args(parser)

main <- function() {
  pipeline_data <- aggregate_data_produced_by_pipeline(args$studies_to_process, args$studies_processed)
  cleanup_studies_with_no_extractions()

  #this should always be the last thing done in the step, as we want to be able to rerun the pipeline other things fail
  vroom::vroom_write(pipeline_data$studies_processed, args$studies_processed)
}

aggregate_data_produced_by_pipeline <- function(studies_to_process_file, studies_processed_file) {
  #update studies_processed.tsv with studies_to_process.tsv
  gene_name_map <- vroom::vroom(glue::glue('{liftover_dir}/gene_name_map.tsv'), show_col_types=F)
  studies_to_process <- vroom::vroom(studies_to_process_file, show_col_types=F)

  if (file.exists(studies_processed_file)) {
    studies_processed <- vroom::vroom(studies_processed_file, show_col_types=F)
    studies_processed <- rbind(studies_processed, studies_to_process) |> dplyr::distinct()
  } else {
    studies_processed <- studies_to_process
  }

  studies_processed$gene <- sub('\\..*', '', studies_processed$gene)
  gene_names <- gene_name_map$GENE_NAME[match(studies_processed$gene, gene_name_map$ENSEMBL_ID)]
  gene_names[is.na(gene_names)] <- studies_processed$gene[is.na(gene_names)]
  studies_processed$gene <- gene_names

  return(list(
    studies_processed = studies_processed
  ))
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
