source('constants.R')
source('svg_helpers.R')

parser <- argparser::arg_parser('Create static web files')
parser <- argparser::add_argument(parser, '--studies_db_file', help = 'Studies DB file', type = 'character')
parser <- argparser::add_argument(parser, '--pipeline_summary_file', help = 'Pipeline summary file', type = 'character')
parser <- argparser::add_argument(parser, '--opengwas_ids_file', help = 'OpenGWAS IDs file', type = 'character')
parser <- argparser::add_argument(parser, '--svg_files_ready_file', help = 'SVG files ready file', type = 'character')
args <- argparser::parse_args(parser)

main <- function() {
  dir.create(dirname(args$opengwas_ids_file), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(args$pipeline_summary_file), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(args$svg_files_ready_file), recursive = TRUE, showWarnings = FALSE)
  
  create_opengwas_map(args$studies_db_file, args$opengwas_ids_file)
  prepare_svg_files_for_use()

  markdown::render('pipeline_summary.Rmd', output_file = args$pipeline_summary_file)
  vroom::vroom_write(data.frame(), args$svg_files_ready_file)
}

create_opengwas_map <- function(studies_db_file, opengwas_ids_file) {
  studies_conn <- duckdb::dbConnect(duckdb::duckdb(), studies_db_file)

  opengwas_ids <- DBI::dbGetQuery(studies_conn, "SELECT study_name FROM studies
    JOIN study_sources ON studies.source_id = study_sources.id
    WHERE source IN ('ebi_catalog', 'ukb')
    AND data_type = 'phenotype' AND variant_type = 'common'")
  jsonlite::write_json(opengwas_ids$study_name, opengwas_ids_file, auto_unbox = TRUE, pretty = TRUE)
}

main()