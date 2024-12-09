source('../pipeline_steps/constants.R')
library(parallel)

ld_blocks <- vroom::vroom('../pipeline_steps/data/ld_blocks.tsv')
ld_info <- construct_ld_block(ld_blocks$ancestry, ld_blocks$chr, ld_blocks$start, ld_blocks$stop)
ld_info <- dplyr::filter(ld_info, dir.exists(ld_block_data))

# TODO: update this method accordingly with how you want to change the data in already ingested studies
update_method <- function(studies_file) {
  if (!file.exists(studies_file)) {
    print(paste('FILE MISSING:', studies_file))
    return()
  }
  #TODO: implement update here

}

cores_to_use <- 128
update_data_in_ld_blocks <- parallel::mclapply(X=ld_info$ld_block_data, mc.cores=cores_to_use, FUN=function(ld_block) {
  print(glue::glue('block: {ld_block}\n'))
  extracted_studies_file <- glue::glue('{ld_block}/extracted_studies.tsv')
  update_method(extracted_studies_file)

  standardised_studies_file <- glue::glue('{ld_block}/standardised_studies.tsv')
  update_method(standardised_studies_file)

  imputed_studies_file <- glue::glue('{ld_block}/imputed_studies.tsv')
  update_method(imputed_studies_file)

  finemapped_studies_file <- glue::glue('{ld_block}/finemapped_studies.tsv')
  update_method(finemapped_studies_file)
})
