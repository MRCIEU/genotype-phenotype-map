setwd('pipeline_steps')
source('constants.R')

minimum_snps_in_opengwas_study <- 1000000

study_list <- vroom::vroom('data/study_list.csv', show_col_types=F)
studies_to_ignore <- vroom::vroom('data/ignore_studies.tsv', delim='\t', show_col_types=F)
studies_processed_file <- glue::glue('{results_dir}studies_processed.tsv')

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
  validate_study_list(study_list)

  opengwas_entries <- dplyr::filter(study_list, data_format == data_formats$opengwas)
  besd_entries <- dplyr::filter(study_list, data_format == data_formats$besd)
  tsv_entries <- dplyr::filter(study_list, data_format == data_formats$tsv)

  opengwas_studies_to_process <- calculate_opengwas_studies_to_process(opengwas_entries)
  besd_studies_to_process <- calculate_besd_studies_to_process(besd_entries)
  tsv_studies_to_process <- calculate_tsv_studies_to_process(tsv_entries)


  # Filter out studies that have already been processed or are to be ignored
  studies_to_process <- dplyr::bind_rows(opengwas_studies_to_process, besd_studies_to_process, tsv_studies_to_process)

  duplicated_study_names <- duplicated(studies_to_process$study_name)
  if (any(duplicated_study_names)) {
    stop(glue::glue('Error: study names are not unique {studies_to_process$study_name[duplicated_study_names]}'))
  }

  studies_to_process <- studies_to_process |>
    dplyr::filter(!study_name %in% studies_processed$study_name) |>
    dplyr::filter(!study_name %in% studies_to_ignore$study)

  lapply(studies_to_process$extracted_location, function(extracted_location) {
    dir.create(extracted_location, showWarnings = F, recursive = T)
  })

  message(paste('Found', nrow(studies_to_process), 'new studies to process'))
  vroom::vroom_write(studies_to_process, glue::glue('{pipeline_metadata_dir}/studies_to_process.tsv'))
}

validate_study_list <- function(study_list) {
  study_sources <- vroom::vroom('data/study_sources.csv', show_col_types=F)

  if (!all(study_list$source %in% study_sources$source)) stop('Error: study_list and study_sources are not compatible')
  if (!all(study_list$data_type %in% ordered_data_types)) stop('Error: some data_type values in study_list are not valid')
  if (!all(study_list$data_format %in% data_formats)) stop('Error: some data_format values in study_list are not valid')
  if (!all(study_list$variant_type %in% variant_types)) stop('Error: some variant_type values in study_list are not valid')
  if (!all(study_list$reference_build %in% reference_builds)) stop('Error: some reference_build values in study_list are not valid')
}

#' calculate_besd_studies_to_process
#' study ids should be in the format of <data_source>-<specifier>-<gene>
#' ex: gtex-liver-ENSG00000269981
#' NOTE: all underscores _ will be turned into dashes -
calculate_besd_studies_to_process <- function(entries) {
  if (nrow(entries) == 0) return(data.frame())
  expanded_studies <- apply(entries, 1, function(entry) {
    file_regex <- glue::glue('{entry[["data_location"]]}/{entry[["id_pattern"]]}')
    data_source <- basename(entry['data_location'])
    all_studies <- Sys.glob(file_regex)
    all_studies <- all_studies[!file.info(all_studies)$isdir]
    studies_without_extensions <- unique(tools::file_path_sans_ext(all_studies))

    if (length(all_studies) == 0) return(data.frame())

    return(data.frame(
      data_type = entry[['data_type']],
      source = entry[['source']],
      data_format = entry[['data_format']],
      bespoke_parsing = entry[['bespoke_parsing']],
      reference_build = entry[['reference_build']],
      data_source = data_source,
      study = studies_without_extensions,
      directory = entry[['data_location']],
      p_value_threshold = entry[['p_value_threshold']],
      ancestry = entry[['ancestry']],
      variant_type = entry[['variant_type']]
    ))
  })|> dplyr::bind_rows()

  if (length(expanded_studies) == 0) return(data.frame())

  processing_information <- apply(expanded_studies, 1, function(besd_study) {
    all_files_present <- Sys.glob(glue::glue('{besd_study["study"]}.*'))
    if (length(all_files_present) < 4) {
      stop(paste('BESD study must include besd, epi, esi, and json files:', besd_study['study']))
    }

    metadata <- jsonlite::fromJSON(glue::glue('{besd_study["study"]}.json'))

    if (is.null(metadata$tissue) || is.null(metadata$sample_size)) {
      stop(paste('json file is missing tissue or sample_size', besd_study['study']))
    }

    epi <- vroom::vroom(glue::glue('{besd_study["study"]}.epi'), col_names = F, show_col_types = F)
    probes <- epi$X2
    genes <- epi$X5
    specifier <- basename(besd_study['study'])
    studies <- paste(besd_study['data_source'], specifier, probes, sep = '-')
    studies <- gsub('_', '-', studies)
    studies <- gsub('\\.', '-', studies)

    data_study_dir <- glue::glue('{data_dir}study/{studies}/')

    if (besd_study['bespoke_parsing'] == bespoke_parsing_options$gtex_sqtl) {
      probe_strings <- sub('chr\\d+:(\\d+):(\\d+):clu_(\\d+):.*', 'cluster:\\3 bp:\\1-\\2', probes)
      traits <- paste(besd_study['data_source'], gsub('[-_]', ' ', specifier), genes, probe_strings)
    } else {
      traits <- paste(besd_study['data_source'], gsub('[-_]', ' ', specifier), genes)
    }
    category <- ifelse(is.null(metadata$category), study_categories$continuous, tolower(metadata$category))
    tissue <- ifelse(is.null(metadata$tissue), NA, metadata$tissue)

    return(data.frame(
      data_type = besd_study[['data_type']],
      source = besd_study[['source']],
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
      variant_type = besd_study[['variant_type']],
      probe = probes,
      gene = genes,
      tissue = metadata$tissue
    ))
  }) |> dplyr::bind_rows()
}

calculate_opengwas_studies_to_process <- function(entries) {
  if (nrow(entries) == 0) return(data.frame())
  expanded_studies <- apply(entries, 1, function(entry) {
    file_regex <- glue::glue('{entry[["data_location"]]}/{entry[["id_pattern"]]}')
    all_directories <- Sys.glob(file_regex)
    return(data.frame(
      data_type = entry[['data_type']],
      source = entry[['source']],
      directory = all_directories,
      reference_build = entry[['reference_build']],
      ancestry = entry[['ancestry']],
      p_value_threshold = entry[['p_value_threshold']],
      data_format = entry[['data_format']],
      variant_type = entry[['variant_type']]
    ))
  }) |> dplyr::bind_rows()

  processing_information <- apply(expanded_studies, 1, function(opengwas_study) {
    directory <- opengwas_study[['directory']]
    study <- basename(directory)
    new_study_name <- gsub('_', '-', study)
    data_study_dir <- glue::glue('{data_dir}study/{new_study_name}/')

    study_metadata <- jsonlite::fromJSON(glue::glue('{directory}/{study}.json'))
    ancestry <- study_metadata$population
    category <- study_metadata$category
    if (is.null(category)) category <- NA
    if (is.null(ancestry) || ancestry != ancestry_map[[opengwas_study[['ancestry']]]]) {
      return(data.frame())
    } else if (as.numeric(study_metadata$nsnp) < minimum_snps_in_opengwas_study) {
      return(data.frame())
    }
    if (is.null(study_metadata$sample_size)) {
      return(data.frame())
    }

    return(data.frame(
      data_type = opengwas_study[['data_type']],
      data_format = opengwas_study[['data_format']],
      source = opengwas_study[['source']],
      study_name = new_study_name,
      trait = study_metadata$trait,
      ancestry = reverse_ancestry_map[[ancestry]],
      sample_size = study_metadata$sample_size,
      category = category,
      study_location = directory,
      extracted_location = data_study_dir,
      reference_build = opengwas_study[['reference_build']],
      p_value_threshold = format(opengwas_study[['p_value_threshold']], scientific=FALSE),
      variant_type = opengwas_study[['variant_type']],
      gene = NA,
      probe = NA,
      tissue = NA
    ))
  }) |> dplyr::bind_rows()
}

calculate_tsv_studies_to_process <- function(entries) {
  if (nrow(entries) == 0) return(data.frame())

  expanded_studies <- apply(entries, 1, function(entry) {
    metadata_file <- glue::glue('{entry[["data_location"]]}/metadata.tsv')
    if (!file.exists(metadata_file)) {
      stop(glue::glue('metadata file is missing: {metadata_file}'))
    }

    tsv_metadata <- vroom::vroom(metadata_file, show_col_types=F) |>
      dplyr::filter(grepl(entry[['id_pattern']], study_location))
    data_study_dir <- glue::glue('{data_dir}study/{tsv_metadata$study_name}/')

    return(data.frame(
      data_type = entry[['data_type']],
      data_format = entry[['data_format']],
      source = entry[['source']],
      study_name = tsv_metadata$study_name,
      trait = tsv_metadata$trait,
      ancestry = entry[['ancestry']],
      sample_size = tsv_metadata$sample_size,
      category = tsv_metadata$category,
      study_location = tsv_metadata$study_location,
      extracted_location = data_study_dir,
      reference_build = entry[['reference_build']],
      p_value_threshold = format(entry[['p_value_threshold']], scientific=FALSE),
      variant_type = entry[['variant_type']],
      gene = ifelse(entry[['data_type']] != ordered_data_types$phenotype, tsv_metadata$gene, NA),
      probe = NA,
      tissue = ifelse(entry[['data_type']] != ordered_data_types$phenotype, tsv_metadata$tissue, NA)
    ))
  }) |> dplyr::bind_rows()

  return(expanded_studies)
}

main()
