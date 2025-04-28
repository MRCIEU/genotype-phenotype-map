setwd('pipeline_steps')
source('constants.R')
source('common_extraction_functions.R')

library(futile.logger)

log_dir <- file.path(data_dir, "logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(log_dir, paste0("pipeline_worker_", format(Sys.time(), "%Y%m"), ".log"))

flog.appender(appender.tee(log_file))
flog.threshold(INFO)

parser <- argparser::arg_parser('Pipeline worker')
parser <- argparser::add_argument(parser, '--reprocess_dlq', help = 'Reprocess DLQ messages', type = 'logical', default = F)
parser <- argparser::add_argument(parser, '--custom_message_file', help = 'Custom message to process (if testing)', type = 'character', default = NULL)
args <- argparser::parse_args(parser)

redis_conn <- redux::hiredis(
  host = Sys.getenv("REDIS_HOST", "redis"),
  port = as.numeric(Sys.getenv("REDIS_PORT", 6379))
)

process_gwas <- 'process_gwas'
process_gwas_dlq <- glue::glue('{process_gwas}_dlq')

main <- function() {
  flog.info("Starting pipeline worker")
  
  if (args$reprocess_dlq) {
    flog.info("Reprocessing DLQ messages")
    retry_dlq_messages()
    return()
  }

  while(TRUE) {
    if (is.na(TEST_RUN)) {
      redis_message <- redis_conn$BRPOP(process_gwas, timeout = 0)
    } else {
      redis_message <- jsonlite::fromJSON(args$custom_message_file)
    }
    
    if (!is.null(redis_message)) {
      gwas_info <- jsonlite::fromJSON(redis_message[[2]])
      flog.info(paste("Received new message from queue with guid:", gwas_info$metadata$guid))

      processed <- process_message(gwas_info)

      if (!is.null(processed)) {
        flog.error(paste("Failed to process message with guid:", gwas_info$metadata$guid))

        message_and_error <- list(
          original_message = gwas_info,
          error = processed,
          timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
        )
        redis_conn$LPUSH(process_gwas_dlq, jsonlite::toJSON(message_and_error, auto_unbox = TRUE))
      } else {
        flog.info(paste("Successfully processed message with guid:", gwas_info$metadata$guid))
      }
    }
    Sys.sleep(5)
    if (!is.na(TEST_RUN)) {
      break
    }
  }
}

process_message <- function(gwas_info) {
  tryCatch({
    update_directories_for_worker(gwas_info$metadata$guid)

    dir.create(pipeline_metadata_dir, recursive = T, showWarnings = F)
    dir.create(extracted_study_dir, recursive = T, showWarnings = F)
    dir.create(ld_block_data_dir, recursive = T, showWarnings = F)

    gwas_info$metadata$file_location <- gwas_info$file_location
    if (grepl('\\.vcf', gwas_info$file_location)) {
      gwas_info$metadata$file_type <- 'vcf'
    } else {
      gwas_info$metadata$file_type <- 'csv'
    }

    flog.info('Verifying GWAS data')
    verification_result <- verify_gwas_data(gwas_info)
    if (!verification_result$valid) {
      flog.error(verification_result$error)
      send_update_gwas_upload(gwas_info, FALSE, paste("bad_data:",verification_result$error))
      return()
    }

    flog.info('Extracting regions')
    create_study_metadata_files(gwas_info)
    extract_regions <- glue::glue("Rscript extract_regions_from_summary_stats.R",
      " --worker_guid {gwas_info$metadata$guid}")
    system(extract_regions, wait = T)

    check_step_complete(glue::glue('{extracted_study_dir}/extracted_snps.tsv'), 'extracted_snps.tsv')
    extracted_regions <- vroom::vroom(glue::glue('{extracted_study_dir}/extracted_snps.tsv'), show_col_types = F)

    flog.info('Organising LD blocks')
    ld_blocks_to_colocalise_file <- glue::glue('{pipeline_metadata_dir}/updated_ld_blocks_to_colocalise.tsv')
    organise_ld_blocks <- glue::glue("Rscript organise_extracted_regions_into_ld_blocks.R",
      " --output_file {ld_blocks_to_colocalise_file}",
      " --worker_guid {gwas_info$metadata$guid}")
    system(organise_ld_blocks, wait = T, intern = T)
    check_step_complete(ld_blocks_to_colocalise_file, 'updated_ld_blocks_to_colocalise.tsv')
    ld_blocks_to_colocalise <- vroom::vroom(ld_blocks_to_colocalise_file, show_col_types = F)

    parallel_block_processing <- 4
    blocks <- head(ld_blocks_to_colocalise$ld_block, 10)
    lapply(blocks, function(block) {
      ld_info <- ld_block_dirs(block)

      output_files <- list(
        standardised = glue::glue('{ld_info$ld_block_data}/standardisation_complete'),
        imputed = glue::glue('{ld_info$ld_block_data}/imputation_complete'),
        finemapped = glue::glue('{ld_info$ld_block_data}/finemapping_complete'),
        coloc = glue::glue('{ld_info$ld_block_data}/colocalisation_complete')
      )

      flog.info(paste('Standardising regions for block:', block))
      standardise_regions <- glue::glue("Rscript standardise_studies_in_ld_block.R",
        " --ld_block {block} ", 
        " --completed_output_file {output_files$standardised}",
        " --worker_guid {gwas_info$metadata$guid}")
      system(standardise_regions, wait = T, intern = T)
      check_step_complete(output_files$standardised, block)

      flog.info(paste('Imputing regions for block:', block))
      impute_regions <- glue::glue("Rscript impute_studies_in_ld_block.R",
        " --ld_block {block} ",
        " --completed_output_file {output_files$imputed}",
        " --worker_guid {gwas_info$metadata$guid}")
      system(impute_regions, wait = T, intern = T)
      check_step_complete(output_files$imputed, block)

      flog.info(paste('Finemapping regions for block:', block))
      finemap_regions <- glue::glue("Rscript finemap_studies_in_ld_block.R",
        " --ld_block {block} ",
        " --completed_output_file {output_files$finemapped}",
        " --worker_guid {gwas_info$metadata$guid}")
      system(finemap_regions, wait = T, intern = T)
      check_step_complete(output_files$finemapped, block)

      flog.info(paste('Colocalising regions for block:', block))
      coloc_regions <- glue::glue("Rscript colocalise_studies_in_ld_block.R",
        " --ld_block {block} ",
        " --completed_output_file {output_files$coloc}",
        " --worker_guid {gwas_info$metadata$guid}")
      system(coloc_regions, wait = T, intern = T)
      check_step_complete(output_files$coloc, block)
    })

    flog.info('Compiling results')
    results <- compile_results()

    flog.info('Calling API to update GWAS upload')
    if (is.na(TEST_RUN)) {
      send_update_gwas_upload(gwas_info, TRUE, NULL, results$coloc_results, results$study_extractions)
    } 
    
  }, error = function(e) {
    flog.error(paste("Error processing message:", e$message))
    
    # Create error details
    error_details <- list(
      original_message = gwas_info,
      error = e$message,
      timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    )
    
    redis_conn$LPUSH(
      process_gwas_dlq, 
      jsonlite::toJSON(error_details, auto_unbox = TRUE)
    )

    send_update_gwas_upload(gwas_info, FALSE, e$message, NULL, NULL)

    return(e$message)
  })
}

verify_gwas_data <- function(gwas_info) {
  gwas <- vroom::vroom(gwas_info$file_location, show_col_types = F)
  gwas_info$metadata$column_names <- Filter(\(column) !is.null(column), gwas_info$metadata$column_names)
  gwas <- change_column_names(gwas, gwas_info$metadata$column_names)
  mandatory_columns <- c('CHR', 'BP', 'P', 'EA', 'OA', 'EAF')
  beta_columns <- c('BETA', 'SE')
  or_columns <- c('OR', 'OR_LB', 'OR_UB')
  
  has_beta <- all(beta_columns %in% colnames(gwas))
  has_or <- all(or_columns %in% colnames(gwas))
  
  if (!has_beta && !has_or) {
    return(list(valid = FALSE, error = "Must have either BETA and SE columns or OR, OR_LB and OR_UB columns"))
  }
  
  if (has_beta && has_or) {
    return(list(valid = FALSE, error = "Cannot have both BETA/SE and OR/OR_LB/OR_UB columns"))
  }

  if (has_or) {
    mandatory_columns <- c(mandatory_columns, or_columns)
  } else {
    mandatory_columns <- c(mandatory_columns, beta_columns) 
  }
  missing_columns <- mandatory_columns[!mandatory_columns %in% colnames(gwas)]

  if (length(missing_columns) > 0) {
    return(list(valid = FALSE, error = paste("Missing mandatory columns:", paste(missing_columns, collapse = ", "))))
  }

  if (any(gwas$P < 0, na.rm = T)) {
    return(list(valid = FALSE, error = "Negative P-values found"))
  }

  if (any(gwas$EAF < 0 | gwas$EAF > 1, na.rm = T)) {
    return(list(valid = FALSE, error = "EAF values out of range"))
  }

  if (has_beta && any(gwas$SE < 0, na.rm = T)) {
    return(list(valid = FALSE, error = "Negative SE values found"))
  }

  return(list(valid = TRUE, error = NULL))
}

create_study_metadata_files <- function(gwas_info) {
  study_to_process <- data.frame(
    data_type = ordered_data_types$phenotype,
    data_format = gwas_info$metadata$file_type,
    source = 'user',
    study_name = gwas_info$metadata$guid,
    trait = gwas_info$metadata$name,
    ancestry = gwas_info$metadata$ancestry,
    sample_size = gwas_info$metadata$sample_size,
    category = gwas_info$metadata$category,
    study_location = gwas_info$file_location,
    extracted_location = extracted_study_dir,
    reference_build = gwas_info$metadata$reference_build,
    p_value_threshold = lowest_p_value_threshold,
    variant_type = variant_types$common,
    gene = NA,
    probe = NA,
    tissue = NA
  )
  studies_to_process_file <- glue::glue('{pipeline_metadata_dir}/studies_to_process.tsv')
  vroom::vroom_write(study_to_process, studies_to_process_file)

  file.create(glue::glue('{extracted_study_dir}/study_metadata.json'))
  jsonlite::write_json(gwas_info$metadata, 
                      glue::glue('{extracted_study_dir}/study_metadata.json'),
                      auto_unbox = TRUE,
                      pretty = TRUE)
}

check_step_complete <- function(output_file, ld_block) {
  if (!file.exists(output_file)) {
    error_msg <- glue::glue('Step {output_file} failed')
    flog.error(error_msg)
    rlang::abort(
      message = error_msg,
      class = "pipeline_worker_error",
      data = list(
        output_file = output_file,
        ld_block = ld_block
      )
    )
  } else {
    flog.info(glue::glue('Step {output_file} completed successfully'))
  }
}

compile_results <- function() {
  ld_block_dirs <- list.dirs(ld_block_data_dir, recursive = TRUE, full.names = TRUE) |>
    {\(dirs) dirs[!dirs %in% dirname(dirs[-1])]}()
  
  compiled_coloc_results_file <- glue::glue('{extracted_study_dir}/compiled_coloc_results.tsv')
  compiled_study_extractions_file <- glue::glue('{extracted_study_dir}/compiled_extracted_studies.tsv')
  
  flog.info("Loading SNP annotations")
  snp_annotations <- vroom::vroom( file.path(variant_annotation_dir, "vep_annotations_hg38.tsv.gz"), show_col_types = FALSE) |>
    dplyr::rename_with(tolower) |>
    dplyr::mutate( snp = trimws(snp), chr = as.character(chr), bp = as.numeric(bp)) |>
    dplyr::select(snp, chr, bp)
  
  finemapped_studies_files <- Filter(
    function(file) file.exists(file), 
    glue::glue('{ld_block_dirs}/finemapped_studies.tsv')
  )
  
  if (length(finemapped_studies_files) > 0) {
    flog.info(paste("Processing", length(finemapped_studies_files), "finemapped study files"))
    study_extractions <- vroom::vroom( finemapped_studies_files, show_col_types = FALSE, col_types = finemapped_column_types) |>
      dplyr::filter(min_p <= p_value_threshold) |>
      dplyr::left_join(snp_annotations, by = c("chr" = "chr", "bp" = "bp")) |>
      dplyr::filter(!duplicated(unique_study_id)) |>
      dplyr::select(study, snp, unique_study_id, file, chr, bp, min_p, cis_trans, ld_block)
  } else {
    flog.warn("No finemapped study files found")
    study_extractions <- data.frame()
  }
  
  coloc_input_files <- Filter(
    function(file) file.exists(file), 
    glue::glue('{ld_block_dirs}/coloc_results.tsv')
  )
  
  if (length(coloc_input_files) > 0) {
    flog.info(paste("Processing", length(coloc_input_files), "coloc result files"))
    coloc_results <- vroom::vroom( coloc_input_files, delim = '\t', show_col_types = FALSE) |>
      dplyr::filter(!is.na(traits) & traits != 'None') |>
      dplyr::filter(
        stringr::str_detect( traits, paste(paste0("\\b", study_extractions$unique_study_id, "\\b"), collapse = "|"))
      )
    
    if (nrow(coloc_results) > 0) {
      coloc_results <- coloc_results |>
        dplyr::filter(posterior_prob > 0.5) |>
        dplyr::mutate( coloc_group_id = 1:dplyr::n(), candidate_snp = trimws(candidate_snp)) |>
        tidyr::separate_longer_delim(cols = traits, delim = ", ") |>
        dplyr::rename(unique_study_id = traits) |>
        dplyr::group_by(unique_study_id, candidate_snp) |>
        dplyr::slice_max(posterior_prob, n = 1, with_ties = FALSE) |>
        dplyr::ungroup()
    }
  } else {
    flog.warn("No coloc result files found")
    coloc_results <- data.frame()
  }
  
  flog.info("Writing compiled results to files")
  vroom::vroom_write(coloc_results, compiled_coloc_results_file)
  vroom::vroom_write(study_extractions, compiled_study_extractions_file)
  
  return(list(
    coloc_results = coloc_results,
    study_extractions = study_extractions
  ))
}


send_update_gwas_upload <- function(gwas_info, success, failure_reason, coloc_results, study_extractions) {
  api_url <- glue::glue("https://gpmap.opengwas.io/api/v1/gwas/{gwas_info$metadata$guid}")
  flog.info(paste("Sending update to API:", api_url))

  if (success) {
    put_body <- list(
      success = success,
      coloc_results = if (is.null(coloc_results) || is.na(coloc_results)) list() else coloc_results,
      study_extractions = if (is.null(study_extractions) || is.na(study_extractions)) list() else study_extractions
    )
  } else {  
    put_body <- list(
      success = success,
      failure_reason = failure_reason
    )
  }
  
  response <- httr::PUT(
    url = api_url,
    body = jsonlite::toJSON(put_body, auto_unbox = TRUE),
    httr::add_headers("Content-Type" = "application/json")
  )
  
  if (httr::status_code(response) != 200) {
    error_msg <- paste("Error updating GWAS:", httr::content(response, "text"))
    flog.error(error_msg)
    stop(error_msg)
  }
    
  flog.info("Successfully updated GWAS results")
}

main()
