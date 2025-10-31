source('../pipeline_steps/constants.R')
source('../pipeline_steps/gwas_calculations.R')
options(scipen = 999)

parser <- argparser::arg_parser('Compare output of steps in the region')
parser <- argparser::add_argument(parser, '--ld_block', help = 'LD block that the ', type = 'character')
parser <- argparser::add_argument(parser, '--study', help = 'Name of study', type = 'character')
args <- argparser::parse_args(parser)

main <- function() {
  ld_info <- ld_block_dirs(args$ld_block)
  extracted_study <- vroom::vroom(glue::glue('{ld_info$ld_block_data}/extracted_studies.tsv'), show_col_types=F) |>
    dplyr::filter(study == args$study)
  print(extracted_study$file)
  extracted_gwas <- vroom::vroom(extracted_study$file, show_col_types=F) |>
    dplyr::mutate(SNP = RSID)

  standardised_study <- vroom::vroom(glue::glue('{ld_info$ld_block_data}/standardised_studies.tsv'), show_col_types=F) |>
    dplyr::filter(study == args$study)
  standardised_gwas <- vroom::vroom(standardised_study$file, show_col_types=F)

  imputed_study <- vroom::vroom(glue::glue('{ld_info$ld_block_data}/imputed_studies.tsv'), show_col_types=F) |>
    dplyr::filter(study == args$study)
  imputed_gwas <- vroom::vroom(imputed_study$file, show_col_types=F)

  finemapped_studies <- vroom::vroom(glue::glue('{ld_info$ld_block_data}/finemapped_studies.tsv'), show_col_types=F) |>
    dplyr::filter(study == args$study)

  finemapped_gwases <- lapply(finemapped_studies$file, function(file) {
    vroom::vroom(file, show_col_types=F)
  })

  manhattan(extracted_gwas, glue::glue('{extracted_study$study}_extracted.png'))
  manhattan(standardised_gwas, glue::glue('{standardised_study$study}_standardised.png'))
  manhattan(imputed_gwas, glue::glue('{imputed_study$study}_imputed.png'))
  apply(finemapped_studies, 1, function(study) {
    print(study["file"])
    gwas <- vroom::vroom(study["file"], show_col_types=F) |>
      dplyr::mutate(P = 2 * pnorm(-abs(convert_lbf_to_abs_z(LBF, SE))))
    png_filename <- gsub('/', '_', glue::glue('{study["unique_study_id"]}_finemapped.png'))
    manhattan(gwas, png_filename)
  })
}

manhattan <- function(gwas, manhattan_filename) {
  print(manhattan_filename)
  print(glue::glue('rows: {nrow(gwas)}'))
  manhattan_columns <- c("SNP", "CHR", "BP", "P")

  gwas$P[gwas$P == 0] <- .Machine$double.xmin
  x_range <- c(min(gwas$BP), max(gwas$BP))

  grDevices::png(manhattan_filename, width = 750, height = 250)
  qqman::manhattan(gwas, xlim = x_range, main = "Manhattan plot of GWAS")
  grDevices::dev.off()
}

main()