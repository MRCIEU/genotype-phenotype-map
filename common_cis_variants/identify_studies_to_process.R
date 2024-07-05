source('constants.R')
DEFAULT_P_VALUE_THRESHOLD <- 5e-8

gwas_list <- vroom::vroom('data/gwas_list.csv', show_col_types=F)
if(!dir.exists(pipeline_metadata_dir)) dir.create(pipeline_metadata_dir)

main <- function() {
  opengwas_entries <- dplyr::filter(gwas_list, database == databases$opengwas)
  opengwas_studies_to_process <- calculate_opengwas_studies_to_process(opengwas_entries)

  #other state calculated here, then we can dplyr::bind_rows()
  studies_to_process <- dplyr::bind_rows(opengwas_studies_to_process)

  message(paste('Found', nrow(studies_to_process), 'new studies to process'))
  vroom::vroom_write(studies_to_process, paste0(pipeline_metadata_dir, 'studies_to_process.tsv'))
}


calculate_opengwas_studies_to_process <- function(entries) {
  expanded_directories <- apply(entries, 1, function(entry) {
    file_regex <- paste0(entry[['data_location']], '/', entry[['id_pattern']])
    all_directories <- Sys.glob(file_regex)
    return(data.frame(
      data_type = entry[['data_type']],
      directory = all_directories,
      ancestry = entry[['ancestry']],
      script = entry[['script']]
    ))
  }) |> dplyr::bind_rows()

  processing_information <- apply(expanded_directories, 1, function(entry) {
    directory <- entry[['directory']]
    study <- basename(directory)
    data_study_dir <- paste0(data_dir, 'study/', study, '/')
    studies_processed_file <- paste0(paste0(results_dir, 'studies_processed.tsv'))

    if (file.exists(studies_processed_file)) {
      studies_processed <- vroom::vroom(studies_processed_file, delim='\t', show_col_types=F)

      already_processed <- dplyr::filter(studies_processed, study_name == study)
      if (nrow(already_processed) > 0 & already_processed$p_value_threshold <= DEFAULT_P_VALUE_THRESHOLD) {
        return(data.frame())
      }
    }

    study_metadata <- jsonlite::fromJSON(paste0(directory, '/', study, '.json'))
    ancestry <- study_metadata$population
    category <- study_metadata$category
    if (is.null(category)) category <- NA
    if (is.null(ancestry) || ancestry != ancestry_map[[entry[['ancestry']]]]) {
      return(data.frame())
    }

    return(data.frame(
      data_type = entry[['data_type']],
      study_name = study,
      trait = study_metadata$trait,
      ancestry = reverse_ancestry_map[[ancestry]],
      sample_size = study_metadata$sample_size,
      category = category,
      study_location = directory,
      extracted_location = data_study_dir,
      p_value_threshold = format(DEFAULT_P_VALUE_THRESHOLD, scientific=FALSE),
      associated_gene = NA,
      script = entry[['script']]
    ))
  }) |> dplyr::bind_rows()
}

main()
