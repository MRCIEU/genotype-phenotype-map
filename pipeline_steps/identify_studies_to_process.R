setwd('pipeline_steps')
source('constants.R')

minimum_snps_in_opengwas_study <- 1000000

gene_name_map <- vroom::vroom(paste0(liftover_dir, 'gene_name_map.tsv'), show_col_types=F)
study_list <- vroom::vroom('data/study_list.csv', show_col_types=F)
studies_processed_file <- paste0(paste0(results_dir, 'studies_processed.tsv'))

if (!is.na(TEST_RUN)) {
  study_list <- vroom::vroom('data/test_list.csv', show_col_types=F)
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
    dir.create(extracted_location, showWarnings = F, recursive = T)
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
      bespoke_parsing = entry[['bespoke_parsing']],
      reference_build = entry[['reference_build']],
      data_source = data_source,
      study = studies_without_extensions,
      directory = entry[['data_location']],
      p_value_threshold = entry[['p_value_threshold']],
      ancestry = entry[['ancestry']]
    ))
  })|> dplyr::bind_rows()

  processing_information <- apply(expanded_studies, 1, function(besd_study) {
    all_files_present <- Sys.glob(paste0(besd_study['study'], '.*'))
    if (length(all_files_present) < 4) {
      stop(paste('BESD study must include besd, epi, esi, and json files:', besd_study['study']))
    }

    metadata <- jsonlite::fromJSON(paste0(besd_study['study'], '.json'))

    if (is.null(metadata$tissue) || is.null(metadata$sample_size)) {
      stop(paste('json file is missing tissue or sample_size', besd_study['study']))
    }

    epi <- vroom::vroom(paste0(besd_study['study'], '.epi'), col_names = F, show_col_types = F)
    probes <- epi$X2
    genes <- epi$X5
    specifier <- basename(besd_study['study'])
    studies <- paste(besd_study['data_source'], specifier, probes, sep = '-')
    studies <- gsub('_', '-', studies)
    studies <- gsub('\\.', '-', studies)

    data_study_dir <- paste0(data_dir, 'study/', studies, '/')

    if (besd_study['bespoke_parsing'] == bespoke_parsing_options$gtex_sqtl) {
      probe_strings <- sub('chr\\d+:(\\d+):(\\d+):clu_(\\d+):.*', 'cluster:\\3 bp:\\1-\\2', probes)
      traits <- paste(besd_study['data_source'], gsub('[-_]', ' ', specifier), genes, probe_strings)
    } else {
      traits <- paste(besd_study['data_source'], gsub('[-_]', ' ', specifier), genes)
    }
    category <- ifelse(is.null(metadata$category), study_categories$continuous, metadata$category)
    tissue <- ifelse(is.null(metadata$tissue), NA, metadata$tissue)

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
      reference_build = besd_study[['reference_build']],
      p_value_threshold = format(besd_study[['p_value_threshold']], scientific=FALSE),
      probe = probes,
      gene = genes,
      tissue = metadata$tissue
    ))
  }) |> dplyr::bind_rows()
}

calculate_opengwas_studies_to_process <- function(entries) {
  if (nrow(entries) == 0) return(data.frame())
  expanded_studies <- apply(entries, 1, function(entry) {
    file_regex <- paste0(entry[['data_location']], '/', entry[['id_pattern']])
    all_directories <- Sys.glob(file_regex)
    return(data.frame(
      data_type = entry[['data_type']],
      directory = all_directories,
      reference_build = entry[['reference_build']],
      ancestry = entry[['ancestry']],
      p_value_threshold = entry[['p_value_threshold']],
      data_format = entry[['data_format']]
    ))
  }) |> dplyr::bind_rows()

  processing_information <- apply(expanded_studies, 1, function(opengwas_study) {
    directory <- opengwas_study[['directory']]
    study <- basename(directory)
    data_study_dir <- paste0(data_dir, 'study/', study, '/')

    study_metadata <- jsonlite::fromJSON(paste0(directory, '/', study, '.json'))
    ancestry <- study_metadata$population
    category <- study_metadata$category
    if (is.null(category)) category <- NA
    if (is.null(ancestry) || ancestry != ancestry_map[[opengwas_study[['ancestry']]]]) {
      return(data.frame())
    } else if (as.numeric(study_metadata$nsnp) < minimum_snps_in_opengwas_study) {
      return(data.frame())
    }

    return(data.frame(
      data_type = opengwas_study[['data_type']],
      data_format = opengwas_study[['data_format']],
      study_name = study,
      trait = study_metadata$trait,
      ancestry = reverse_ancestry_map[[ancestry]],
      sample_size = study_metadata$sample_size,
      category = category,
      study_location = directory,
      extracted_location = data_study_dir,
      reference_build = opengwas_study[['reference_build']],
      p_value_threshold = format(opengwas_study[['p_value_threshold']], scientific=FALSE),
      gene = NA,
      probe = NA,
      tissue = NA
    ))
  }) |> dplyr::bind_rows()
}

main()
