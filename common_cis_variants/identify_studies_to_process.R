source("constants.R")
gwas_list <- vroom::vroom("data/gwas_list.csv", show_col_types=F)
if(!dir.exists(pipeline_metadata_dir)) dir.create(pipeline_metadata_dir)
extraction_p_value <- 5e-8

main <- function() {
  opengwas_entries <- dplyr::filter(gwas_list, database == databases$opengwas)
  opengwas_data_current_state <- calculate_state_opengwas_data(opengwas_entries)

  #other state calculated here, then we can dplyr::bind_rows()
  studies_to_process <- dplyr::bind_rows(opengwas_data_current_state)

  message(paste('Found', nrow(studies_to_process), 'new studies to process'))
  vroom::vroom_write(studies_to_process, paste0(pipeline_metadata_dir, "studies_to_process.tsv"))
}


calculate_state_opengwas_data <- function(entries) {
  expanded_directories <- apply(entries, 1, function(entry) {
    file_regex <- paste0(entry[['data_location']], "/", entry[['id_pattern']])
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
    study_name <- basename(directory)
    data_study_dir <- paste0(data_dir, 'study/', study_name, "/")
    study_metadata_json <- paste0(data_study_dir, "extraction_metadata.json")
    current_state <- vroom::vroom(paste0(results_dir, 'current_state.tsv'), delim='\t', show_col_types=F)

    if (file.exists(study_metadata_json)) {
      json_data <- readLines(study_metadata_json)
      if (jsonlite::validate(json_data)) {
        extraction_metadata <- jsonlite::fromJSON(study_metadata_json)
        p_value <- as.numeric(extraction_metadata$p_value_threshold)

        already_processed <- any(current_state$study == study_name)
        if (p_value <= extraction_p_value && already_processed) {
          return(data.frame())
        }
      }
    }

    study_metadata <- jsonlite::fromJSON(paste0(directory, "/", study_name, ".json"))
    ancestry <- study_metadata$population
    category <- study_metadata$category
    if (is.null(category)) category <- NA
    if (is.null(ancestry) || ancestry != ancestry_map[[entry[['ancestry']]]]) {
      return(data.frame())
    }

    return(data.frame(
      data_type = entry[['data_type']],
      study_name = study_name,
      trait = study_metadata$trait,
      ancestry = reverse_ancestry_map[[ancestry]],
      sample_size = study_metadata$sample_size,
      category = category,
      study_location = directory,
      extracted_location = data_study_dir,
      p_value_threshold = format(extraction_p_value, scientific=FALSE),
      associated_gene = NA,
      script = entry[['script']]
    ))
  }) |> dplyr::bind_rows()
}

main()
