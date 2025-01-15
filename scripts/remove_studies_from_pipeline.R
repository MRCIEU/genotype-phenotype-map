source('../pipeline_steps/constants.R')

# studies_to_remove <- vroom::vroom(glue::glue('{results_dir}/studies_processed.tsv')) |> 
studies_to_remove <- vroom::vroom(glue::glue('{pipeline_metadata_dir}/studies_to_process.tsv')) |> 
  dplyr::filter(data_type != 'phenotype')
studies_to_remove <- studies_to_remove$study_name
length(studies_to_remove)

ld_blocks <- vroom::vroom('../pipeline_steps/data/ld_blocks.tsv')
ld_info <- construct_ld_block(ld_blocks$ancestry, ld_blocks$chr, ld_blocks$start, ld_blocks$stop)
ld_info <- dplyr::filter(ld_info, dir.exists(ld_block_data))

safe_lapply <- function(X, FUN, ...) {
  lapply(X, function(x) {
    tryCatch(
      FUN(x, ...),
      error = function(e) {
        stop("An error occurred: ", conditionMessage(e))
      }
    )
  })
}

main <- function() {
  delete_data_in_ld_blocks <- safe_lapply(ld_info$ld_block_data, function(ld_block) {
    print(glue::glue('block: {ld_block}'))
    extracted_studies_file <- glue::glue('{ld_block}/extracted_studies.tsv')
    if (!file.exists(extracted_studies_file) || file.size(extracted_studies_file) == 0) {
      print(paste('extracted STUDIES FILE MISSING:', extracted_studies_file))
      return()
    }

    extracted_studies <- vroom::vroom(extracted_studies_file, show_col_types = F)
    entries <- nrow(extracted_studies)
    extracted_studies <- dplyr::filter(extracted_studies, !study %in% studies_to_remove)
    if (!entries - nrow(extracted_studies) == 0) {
      print(paste('removed', entries - nrow(extracted_studies), 'rows from extracted_studies'))
      vroom::vroom_write(extracted_studies, extracted_studies_file)
    }

    standardised_studies_file <- glue::glue('{ld_block}/standardised_studies.tsv')
    if (!file.exists(standardised_studies_file) || file.size(standardised_studies_file) == 0) {
      print(paste('standardised STUDIES FILE MISSING:', standardised_studies_file))
      return()
    }
    standardised_studies <- vroom::vroom(standardised_studies_file, col_types = standardised_column_types, show_col_types = F)
    entries <- nrow(standardised_studies)
    standardised_studies <- dplyr::filter(standardised_studies, !study %in% studies_to_remove)
    if (!entries - nrow(standardised_studies) == 0) {
      vroom::vroom_write(standardised_studies, standardised_studies_file)
      print(paste('removed', entries - nrow(standardised_studies), 'rows from standardised_studies'))
    }

    imputed_studies_file <- paste0(ld_block, '/imputed_studies.tsv')
    if (!file.exists(imputed_studies_file) || file.size(imputed_studies_file) == 0) {
      print(paste('IMPUTED STUDIES FILE MISSING:', imputed_studies_file))
      return()
    }

    imputed_studies <- vroom::vroom(imputed_studies_file, col_types = imputed_column_types, show_col_types = F)
    entries <- nrow(imputed_studies)
    imputed_studies <- dplyr::filter(imputed_studies, !study %in% studies_to_remove)
    if (!entries - nrow(imputed_studies) == 0) {
      print(paste('removed', entries - nrow(imputed_studies), 'rows from imputed_studies'))
      vroom::vroom_write(imputed_studies, imputed_studies_file)
    }

    finemapped_studies_file <- paste0(ld_block, '/finemapped_studies.tsv')
    if (!file.exists(finemapped_studies_file) || file.size(finemapped_studies_file) == 0) {
      print(paste('FINEMAPPED STUDIES FILE MISSING:', finemapped_studies_file))
      return()
    }

    finemapped_studies <- vroom::vroom(finemapped_studies_file, col_types = finemapped_column_types, show_col_types = F)
    entries <- nrow(finemapped_studies)
    finemapped_studies <- dplyr::filter(finemapped_studies, !study %in% studies_to_remove)
    if (!entries - nrow(finemapped_studies) == 0) {
      print(paste('removed', entries - nrow(finemapped_studies), 'rows from finemapped_studies'))
      vroom::vroom_write(finemapped_studies, finemapped_studies_file)
    }
  })

  delete_study_directories <- lapply(studies_to_remove, function(study) {
    extracted_study_dir <- glue::glue('{extracted_study_dir}{study}')
    print(paste('deleting:', study))
    unlink(extracted_study_dir, recursive = T)
  })

  #and then delete them from the results
  studies_processed_file <- glue::glue('{results_dir}/studies_processed.tsv')
  studies_processed <- vroom::vroom(studies_processed_file, show_col_types=F)
  entries <- nrow(studies_processed)
  studies_processed <- dplyr::filter(studies_processed, !study_name %in% studies_to_remove)
  print(paste('removed', entries - nrow(studies_processed), 'rows from studies_processed file'))
  vroom::vroom_write(studies_processed, studies_processed_file)

}


main()