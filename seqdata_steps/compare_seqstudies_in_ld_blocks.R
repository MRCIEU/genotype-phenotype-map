source('constants.R')

# In the absence of a complex LD structure surrounding rare variant associations, compare WES and WGS studies by variant identifer alone
# for all studies with a significant hit within each LD block

parser <- argparser::arg_parser('Compare rare variant sequencing studies per LD block')
parser <- argparser::add_argument(parser, '--ld_block', help = 'LD block', type = 'character')
parser <- argparser::add_argument(parser, '--completed_output_file', help = 'Result file to save', type = 'character')
args <- argparser::parse_args(parser)

main <- function() {
    #args <- list(ld_block="EUR/1/95684841-97259500")

    ld_info <- ld_block_dirs(args$ld_block)
    block <- vroom::vroom(glue::glue('{pipeline_metadata_dir}updated_ld_blocks_to_colocalise.tsv'), show_col_types=F) |>
        dplyr::filter(data_dir == ld_info$ld_block_data)

    #standardised_file <- "example_studiesinblock.txt"
    #standardised_studies <- vroom::vroom(standardised_file, delim = '\t', show_col_types = F)
    #colnames(standardised_studies) <- "file"
    #standardised_studies$study <- lapply(stringr::str_split(standardised_studies$file, pattern = "/"), function(x){x[7]}) |> unlist()
    
    standardised_file <- glue::glue('{ld_info$ld_block_data}/standardised_studies.tsv')
    
    if (file.exists(standardised_file)) {
    standardised_studies <- vroom::vroom(standardised_file, delim = '\t', show_col_types = F) |>
        dplyr::filter(varaint_type == "rare")
    }

    if (!file.exists(standardised_file) || nrow(block) == 0 || nrow(standardised_studies) == 0) {
        message(glue::glue('Nothing to process for LD region {ld_info$ld_block_data}, skipping.'))
        vroom::vroom_write(data.frame(), args$completed_output_file)
        return()
    }






