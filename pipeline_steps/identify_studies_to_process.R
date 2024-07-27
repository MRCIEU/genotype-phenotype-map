setwd('pipeline_steps')
source('constants.R')

gene_name_map <- vroom::vroom(paste0(thousand_genomes_dir, 'gene_name_map.tsv'), show_col_types=F)
study_list <- vroom::vroom('data/study_list.csv', show_col_types=F)
studies_processed_file <- paste0(paste0(results_dir, 'studies_processed.tsv'))

if (!is.na(TEST_RUN)) {
  study_list <- vroom::vroom(paste0('data/', TEST_RUN, '_list.csv'), show_col_types=F)
  studies_processed_file <- paste0(paste0(results_dir, TEST_RUN, '_test_studies_processed.tsv'))
}
if (file.exists(studies_processed_file)) {
  studies_processed <- vroom::vroom(studies_processed_file, delim='\t', show_col_types=F)
} else {
  studies_processed <- data.frame(study_name=NA)
}

main <- function() {
  if(!dir.exists(pipeline_metadata_dir)) dir.create(pipeline_metadata_dir)

  opengwas_entries <- dplyr::filter(study_list, data_format == data_formats$opengwas)
  besd_entries <- dplyr::filter(study_list, data_format == data_formats$besd)
  opengwas_studies_to_process <- calculate_opengwas_studies_to_process(opengwas_entries)
  besd_studies_to_process <- calculate_besd_studies_to_process(besd_entries)

  #TODO: check if there exists a study with that study_name already.  Can't be duplicates
  studies_to_process <- dplyr::bind_rows(opengwas_studies_to_process, besd_studies_to_process) |>
    dplyr::filter(!study_name %in% studies_processed$study_name)

  lapply(studies_to_process$extracted_location, function(extracted_location) {
    dir.create(paste0(extracted_location, '/original'), showWarnings = F, recursive = T)
    dir.create(paste0(extracted_location, '/imputed'), showWarnings = F, recursive = T)
    dir.create(paste0(extracted_location, '/finemapped'), showWarnings = F, recursive = T)
  })

  message(paste('Found', nrow(studies_to_process), 'new studies to process'))
  vroom::vroom_write(studies_to_process, paste0(pipeline_metadata_dir, 'studies_to_process.tsv'))
}

#' calculate_besd_studies_to_process
#' study ids should be in the format of <data_source>-<specifier>-<gene>
#' ex: gtex-liver-ENSG00000269981
#' NOTE: all underscores _ will be turned into dashes -
calculate_besd_studies_to_process <- function(entries) {
  if (nrow(entries) == 0) return(data.frame())
  expanded_studies <- apply(entries, 1, function(entry) {
    file_regex <- paste0(entry[['data_location']], '/', entry[['id_pattern']])
    data_source <- basename(entry['data_location'])
    all_studies <- Sys.glob(file_regex)
    studies_without_extensions <- unique(tools::file_path_sans_ext(all_studies))

    return(data.frame(
      data_type = entry[['data_type']],
      data_format = entry[['data_format']],
      data_source = data_source,
      study = studies_without_extensions,
      directory = entry[['data_location']],
      ancestry = entry[['ancestry']]
    ))
  })|> dplyr::bind_rows()

  processing_information <- apply(expanded_studies, 1, function(besd_study) {
    all_files_present <- Sys.glob(paste0(besd_study['study'], '.*'))
    if (length(all_files_present) != 4) {
      stop(paste('BESD study must include besd, epi, esi, and json files:', besd_study['study']))
    }

    metadata <- jsonlite::fromJSON(paste0(besd_study['study'], '.json'))

    probes <- vroom::vroom(paste0(besd_study['study'], '.epi'), col_select = 2, col_names = F, show_col_types = F)$X2
    specifier <- basename(besd_study['study'])
    studies <- paste(besd_study['data_source'], specifier, probes, sep = '-')
    studies <- gsub('_', '-', studies)

    data_study_dir <- paste0(data_dir, 'study/', studies, '/')

    gene_names <- gene_name_map$GENE_NAME[match(probes, gene_name_map$ENSEMBL_ID)]
    gene_names[is.na(gene_names)] <- probes[is.na(gene_names)]
    traits <- paste(besd_study['data_source'], gsub('[-_]', ' ', specifier), gene_names)
    category <- ifelse(is.na(metadata$category), study_categories$continuous, metadata$category)

    return(data.frame(
      data_type = besd_study[['data_type']],
      data_format = besd_study[['data_format']],
      study_name = studies,
      trait = traits,
      ancestry = besd_study[['ancestry']],
      sample_size = metadata$sample_size,
      category = category,
      study_location = besd_study[['study']],
      extracted_location = data_study_dir,
      p_value_threshold = format(DEFAULT_P_VALUE_THRESHOLD, scientific=FALSE),
      gene = probes
    ))
  }) |> dplyr::bind_rows()
}

calculate_opengwas_studies_to_process <- function(entries) {
  if (nrow(entries) == 0) return(data.frame())
  expanded_directories <- apply(entries, 1, function(entry) {
    file_regex <- paste0(entry[['data_location']], '/', entry[['id_pattern']])
    all_directories <- Sys.glob(file_regex)
    return(data.frame(
      data_type = entry[['data_type']],
      directory = all_directories,
      ancestry = entry[['ancestry']],
      data_format = entry[['data_format']]
    ))
  }) |> dplyr::bind_rows()

  processing_information <- apply(expanded_directories, 1, function(entry) {
    directory <- entry[['directory']]
    study <- basename(directory)
    data_study_dir <- paste0(data_dir, 'study/', study, '/')

    study_metadata <- jsonlite::fromJSON(paste0(directory, '/', study, '.json'))
    ancestry <- study_metadata$population
    category <- study_metadata$category
    if (is.null(category)) category <- NA
    if (is.null(ancestry) || ancestry != ancestry_map[[entry[['ancestry']]]]) {
      return(data.frame())
    }

    return(data.frame(
      data_type = entry[['data_type']],
      data_format = entry[['data_format']],
      study_name = study,
      trait = study_metadata$trait,
      ancestry = reverse_ancestry_map[[ancestry]],
      sample_size = study_metadata$sample_size,
      category = category,
      study_location = directory,
      extracted_location = data_study_dir,
      p_value_threshold = format(DEFAULT_P_VALUE_THRESHOLD, scientific=FALSE),
      gene = NA
    ))
  }) |> dplyr::bind_rows()
}

main()
