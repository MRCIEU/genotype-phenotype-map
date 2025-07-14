source('../pipeline_steps/constants.R')
library(parallel)

ld_blocks <- vroom::vroom('../pipeline_steps/data/ld_blocks.tsv')
ld_info <- construct_ld_block(ld_blocks$ancestry, ld_blocks$chr, ld_blocks$start, ld_blocks$stop)
ld_info <- dplyr::filter(ld_info, dir.exists(ld_block_data))

main <- function() {
  # update_study_dirs()
  update_ld_blocks()
  # update_studies_processed()
}

# TODO: update this method accordingly with how you want to change the data in already ingested studies
update_method <- function(studies_file, type) {
  if (!file.exists(studies_file)) {
    print(paste('FILE MISSING:', studies_file))
    return()
  }

  godmc_studies <- vroom::vroom(studies_file, show_col_types = F) |>
    dplyr::filter(grepl('godmc-methylation', study))

  print(paste('updating', type, 'studies for', nrow(godmc_studies), ' godmc studies'))

  lapply(seq_len(nrow(godmc_studies)), function(i) {
    current_study <- godmc_studies[i, , drop = F]

    if (type == 'imputed') {
      gwas <- vroom::vroom(current_study$file, show_col_types = F, delim = '\t')
      if (nrow(gwas) == 0 || !'BETA' %in% colnames(gwas)) {
        print(paste('DELETING: either no rows or no beta or se in file: ', current_study$file))
        vroom::vroom(studies_file, show_col_types = F) |>
          dplyr::filter(study != current_study$study) |>
          vroom::vroom_write(studies_file)
        file.remove(current_study$file)
        return()
      }
      if (!'IMPUTED' %in% colnames(gwas)) {
        gwas$IMPUTED <- ifelse(gwas$BETA == 0 & gwas$SE == 1, T, F)
        vroom::vroom_write(gwas, current_study$file)
      }
    } else if (type == 'finemapped') {
      gwas <- vroom::vroom(current_study$file_with_lbfs, show_col_types = F, delim = '\t')
      if (nrow(gwas) == 0 || !'BETA' %in% colnames(gwas)) {
        print(paste('DELETING: either no rows or no beta or se in file: ', current_study$file_with_lbfs))
        vroom::vroom(studies_file, show_col_types = F) |>
          dplyr::filter(study != current_study$study) |>
          vroom::vroom_write(studies_file)
        file.remove(current_study$file_with_lbfs)
        return()
      }
      if (!'IMPUTED' %in% colnames(gwas)) {
        gwas$IMPUTED <- ifelse(gwas$BETA == 0 & gwas$SE == 1, T, F)
        vroom::vroom_write(gwas, current_study$file_with_lbfs)
      }
    }
    return(current_study)
  })

}

update_study_dirs <- function() {
  print("starting update study dirs")
  ukb_d_dirs <- list.dirs(path = extracted_study_dir, recursive = FALSE, full.names = TRUE)

  # cores_to_use <- 64
  # update_data_in_study_dirs <- parallel::mclapply((X=ukb_d_dirs), mc.cores=cores_to_use, FUN=function(study_dir) {
  result <- lapply(tail(ukb_d_dirs, 1), function(study_dir) {
    print(paste('changing:', study_dir))
    study_id <- basename(study_dir)
    extracted_snps_file <- paste0(study_dir, '/extracted_snps.tsv')
    if (!file.exists(extracted_snps_file )) {
      print(paste('EXTRACTION MISSING: ', extracted_snps_file))
      return()
    }
    # ...
  })
}

update_ld_blocks <- function() {
  blocks <- ld_info$ld_block_data
  # cores_to_use <- 40
  # update_data_in_ld_blocks <- parallel::mclapply(X=blocks, mc.cores=cores_to_use, FUN=function(ld_block) {
  update_data_in_ld_blocks <- lapply(blocks, function(ld_block) {
    print(glue::glue('block: {ld_block}\n'))
    # extracted_studies_file <- glue::glue('{ld_block}/extracted_studies.tsv')
    # update_method(extracted_studies_file, type='extracted')

    # standardised_studies_file <- glue::glue('{ld_block}/standardised_studies.tsv')
    # update_method(standardised_studies_file, type='standardised')

    imputed_studies_file <- glue::glue('{ld_block}/imputed_studies.tsv')
    update_method(imputed_studies_file, type='imputed')

    finemapped_studies_file <- glue::glue('{ld_block}/finemapped_studies.tsv')
    update_method(finemapped_studies_file, type='finemapped')
  })
}

update_studies_processed <- function() {
  print("starting stupdate_studies_processed")
  studies_to_process_file <- glue::glue('{pipeline_metadata_dir}/studies_to_process.tsv') 
  if (file.exists(studies_to_process_file)) {
    studies_to_process <- vroom::vroom(studies_to_process_file, show_col_types = F)
    #...
    # vroom::vroom_write(studies_to_process, studies_to_process_file)
  }

  studies_processed_file <- glue::glue('{results_dir}/studies_processed.tsv') 
  processed_studies <- vroom::vroom(studies_processed_file, show_col_types = F)
  processed_studies$variant_type <- variant_types$common
  # ...
  vroom::vroom_write(processed_studies, studies_processed_file)

}

main()