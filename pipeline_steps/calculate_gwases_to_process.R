source("../common_cis_variants/constants.R")
gwas_list <- vroom::vroom("../common_cis/variants/data/gwas_list.csv")
if(!dir.exists(pipeline_metadata_dir)) dir.create(pipeline_metadata_dir)
extraction_p_value <- 5e-8


#TODO: calculate all the different files per study that you might need (imputed file, finampped file)

calculate_state_opengwas_data <- function(entries) {
  expanded_directories <- apply(entries, 1, function(entry) {
    file_regex <- paste0(entry[['data_location']], "/", entry[['id_pattern']])
    print(file_regex)
    print(entry[['data_type']])
    print(entry[['ancestry']])
    all_directories <- Sys.glob(file_regex)
    return(data.frame(
      data_type = entry[['data_type']],
      directory = all_directories,
      ancestry = entry[['ancestry']]
    ))
  }) |> dplyr::bind_rows()

  processing_information <- apply(expanded_directories, 1, function(entry) {
    directory <- entry[['directory']]
    study_name <- basename(directory)
    data_study_dir <- paste0(data_dir, 'study/', study_name, "/")
    study_metadata_json <- paste0(data_study_dir, "extraction_metadata.json")

    if (file.exists(study_metadata_json)) {
      extraction_metadata <- jsonlite::fromJSON(study_metadata_json)
      p_value <- as.numeric(extraction_metadata$p_value_threshold)
      if (p_value <= extraction_p_value) return(data.frame())
    } else {
      #TODO: maybe need to delete the extraction_metadata.json file, so it can be reprocessed?
    }

    study_metadata <- jsonlite::fromJSON(paste0(directory, "/", study_name, ".json"))
    ancestry <- study_metadata$population
    sample_size <- study_metadata$sample_size
    category <- study_metadata$category
    if (is.null(category)) category <- NA
    if (is.null(ancestry) || ancestry != ancestry_map[[entry[['ancestry']]]]) {
      return(data.frame())
    }

    return(data.frame(
      data_type = entry[['data_type']],
      study_name = study_name,
      ancestry = reverse_ancestry_map[[ancestry]],
      sample_size = sample_size,
      category = category,
      study_location = directory,
      extracted_location = data_study_dir,
      p_value_threshold = format(extraction_p_value, scientific=FALSE)
    ))
  }) |> dplyr::bind_rows()
}

opengwas_entries <- dplyr::filter(gwas_list, database == databases$opengwas)
opengwas_data_current_state <- calculate_state_opengwas_data(opengwas_entries)

#other state calculated here, then we can dplyr::bind_rows()
current_state <- dplyr::bind_rows(opengwas_data_current_state)

vroom::vroom_write(current_state, paste0(pipeline_metadata_dir, "gwases_to_process.tsv"))
