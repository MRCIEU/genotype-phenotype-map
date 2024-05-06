source("constants.R")
gwas_list <- vroom::vroom("data/test_list.csv")

calculate_state_opengwas_data <- function(entries) {
  print(entries)
  expanded_directories <- apply(entries, 1, function(entry) {
    print(entry)
    file_regex <- paste0(entry[['data_location']], entry[['id_pattern']])
    all_directories <- Sys.glob(file_regex)
    return(data.frame(
      data_type = entry[['data_type']],
      directory = all_directories
    ))
  }) |> dplyr::bind_rows()

  processing_information <- apply(expanded_directories, 1, function(entry) {
    directory <- entry['directory']
    study_name <- basename(directory)
    data_study_dir <- paste0(data_dir, study_name, "/")
    finemap_dir <- paste0(data_study_dir, "finemap/")

    clumped_rsids <- data.table::fread(paste0(directory, "/clump.txt"), header = F)$V1
    metadata <- jsonlite::fromJSON(paste0(directory, "/", study_name, ".json"))
    ancestry <- metadata$pop #TODO: fix this

    already_extracted <- file.exists(data_study_dir, clumped_rsids, ".tsv.gz")
    already_finemapped <- lapply(clumped_rsids, function(rsid) {
      finemap_file_exists <- Sys.glob(paste0(finemap_dir, rsid, "*"))
      return(finemap_file_exists)
    })

    return(data.frame(
      data_type = entry[['data_type']],
      study_name = study_name,
      ancestry = ancestry,
      study_location = directory,
      extracted_location = data_study_dir,
      rsid = clumped_rsids,
      extracted = already_extracted,
      finemapped = already_finemapped
    ))
  }) |> dplyr::bind_rows()
}

print(gwas_list)
opengwas_entries <- dplyr::filter(gwas_list, database == databases$opengwas)
opengwas_data_current_state <- calculate_state_opengwas_data(opengwas_entries)
#other state calculated here, then we can dplyr::bind_rows()
current_state <- dplyr::bind_rows(opengwas_data_current_state)
vroom::vroom_write(current_state, paste0(pipeline_metadata_dir, "current_state.tsv"))

gwases_to_extract <- current_state[current_state$extracted == T, c("study_location", "extracted_location", "rsid")]
vroom::vroom_write(gwases_to_extract, paste0(pipeline_metadata_dir, "gwases_to_extract.tsv", col_names = F))

gwases_to_finemap <- current_state[current_state$finmapped == T, c("study_location", "extracted_location", "rsid")]
vroom::vroom_write(gwases_to_finemap, paste0(pipeline_metadata_dir, "gwases_to_finemap.tsv", col_names = F))
