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


if (is.na(TEST_RUN)) {
  redis_conn <- redux::hiredis(
    host = Sys.getenv("REDIS_HOST", "redis"),
    port = as.numeric(Sys.getenv("REDIS_PORT", 6379))
  )
} else {
  flog.appender(appender.console())
  redis_conn <- list(
    BRPOP = function(queue, timeout = 0) {
      return(jsonlite::fromJSON(args$custom_message_file))
    },
    LPUSH = function(queue, message) {
      return(1)
    }
  )
}

snp_annotations <- vroom::vroom(
  file.path(variant_annotation_dir, "vep_annotations_hg38.tsv.gz"),
  show_col_types = FALSE,
  col_select = c('snp', 'rsid', 'display_snp')
) |>
  dplyr::mutate(snp = trimws(snp))

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
    redis_message <- redis_conn$BRPOP(process_gwas, timeout = 0)
    
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
    output <-system(extract_regions, wait = T, intern = T)

    check_step_complete(glue::glue('{extracted_study_dir}/extracted_snps.tsv'), 'extracted_snps.tsv', output)
    extracted_regions <- vroom::vroom(glue::glue('{extracted_study_dir}/extracted_snps.tsv'), show_col_types = F)

    flog.info('Organising LD blocks')
    ld_blocks_to_colocalise_file <- glue::glue('{pipeline_metadata_dir}/updated_ld_blocks_to_colocalise.tsv')
    organise_ld_blocks <- glue::glue("Rscript organise_extracted_regions_into_ld_blocks.R",
      " --output_file {ld_blocks_to_colocalise_file}",
      " --worker_guid {gwas_info$metadata$guid}")
    output <- system(organise_ld_blocks, wait = T, intern = T)
    check_step_complete(ld_blocks_to_colocalise_file, 'updated_ld_blocks_to_colocalise.tsv', output)
    ld_blocks_to_colocalise <- vroom::vroom(ld_blocks_to_colocalise_file, show_col_types = F)

    parallel_block_processing <- 4
    blocks <- head(ld_blocks_to_colocalise$ld_block, 10)
    lapply(blocks, function(block) {
      ld_info <- ld_block_dirs(block)

      output_files <- list(
        standardised = glue::glue('{ld_info$ld_block_data}/standardisation_complete'),
        imputed = glue::glue('{ld_info$ld_block_data}/imputation_complete'),
        finemapped = glue::glue('{ld_info$ld_block_data}/finemapping_complete'),
        coloc = glue::glue('{ld_info$ld_block_data}/coloc_complete')
      )

      flog.info(paste('Standardising regions for block:', gwas_info$metadata$guid, block))
      standardise_regions <- glue::glue("Rscript standardise_studies_in_ld_block.R",
        " --ld_block {block} ", 
        " --completed_output_file {output_files$standardised}",
        " --worker_guid {gwas_info$metadata$guid}")
      output <- system(standardise_regions, wait = T, intern = T)
      check_step_complete(output_files$standardised, block, output)

      flog.info(paste('Imputing regions for block:', gwas_info$metadata$guid, block))
      impute_regions <- glue::glue("Rscript impute_studies_in_ld_block.R",
        " --ld_block {block} ",
        " --completed_output_file {output_files$imputed}",
        " --worker_guid {gwas_info$metadata$guid}")
      output <- system(impute_regions, wait = T, intern = T)
      check_step_complete(output_files$imputed, block, output)

      flog.info(paste('Finemapping regions for block:', gwas_info$metadata$guid, block))
      finemap_regions <- glue::glue("Rscript finemap_studies_in_ld_block.R",
        " --ld_block {block} ",
        " --completed_output_file {output_files$finemapped}",
        " --worker_guid {gwas_info$metadata$guid}")
      output <- system(finemap_regions, wait = T, intern = T)
      check_step_complete(output_files$finemapped, block, output)

      flog.info(paste('Colocalising regions for block:', gwas_info$metadata$guid, block))
      coloc_regions <- glue::glue("Rscript coloc_and_cluster_studies_in_ld_block.R",
        " --ld_block {block} ",
        " --completed_output_file {output_files$coloc}",
        " --worker_guid {gwas_info$metadata$guid}")
      output <- system(coloc_regions, wait = T, intern = T)
      check_step_complete(output_files$coloc, block, output)
    })

    flog.info(paste('Compiling results for:', gwas_info$metadata$guid))
    results <- compile_results(gwas_info)

    if (is.na(TEST_RUN)) {
      send_update_gwas_upload(gwas_info, TRUE, NULL, results$coloc_results, results$study_extractions)
    } 
    
  }, error = function(e) {
    error_msg <- if (!is.null(e$message) && nchar(e$message) > 0) {
      e$message
    } else {
      paste("Error:", toString(e))
    }
    
    flog.error(paste("Error processing message:", error_msg))
    flog.error(paste("Error class:", class(e)[1]))
    if (!is.null(e$call)) {
      flog.error(paste("Error call:", deparse(e$call)))
    }
    
    # Create error details
    error_details <- list(
      original_message = gwas_info,
      error = error_msg,
      timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    )
    
    redis_conn$LPUSH(
      process_gwas_dlq, 
      jsonlite::toJSON(error_details, auto_unbox = TRUE)
    )

    send_update_gwas_upload(gwas_info, FALSE, error_msg, NULL)

    return(error_msg)
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
    data_type = data_types$phenotype,
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
    p_value_threshold = gwas_info$metadata$p_value_threshold,
    variant_type = variant_types$common,
    gene = NA,
    probe = NA,
    tissue = NA,
    coverage = coverage_types$dense
  )
  studies_to_process_file <- glue::glue('{pipeline_metadata_dir}/studies_to_process.tsv')
  vroom::vroom_write(study_to_process, studies_to_process_file)

  file.create(glue::glue('{extracted_study_dir}/study_metadata.json'))
  jsonlite::write_json(gwas_info$metadata, 
                      glue::glue('{extracted_study_dir}/study_metadata.json'),
                      auto_unbox = TRUE,
                      pretty = TRUE)
}

check_step_complete <- function(output_file, ld_block, output) {
  if (!file.exists(output_file)) {
    error_msg <- glue::glue('Step {output_file} failed')
    flog.error(error_msg)
    flog.error(output)
    rlang::abort(
      message = error_msg,
      class = "pipeline_worker_error",
      data = list(
        output_file = output_file,
        ld_block = ld_block
      )
    )
  }
}

compile_results <- function(gwas_info) {
  ld_block_dirs <- list.dirs(ld_block_data_dir, recursive = TRUE, full.names = TRUE) |>
    {\(dirs) dirs[!dirs %in% dirname(dirs[-1])]}()
  
  compiled_coloc_results_file <- glue::glue('{extracted_study_dir}/compiled_coloc_results.tsv')
  compiled_study_extractions_file <- glue::glue('{extracted_study_dir}/compiled_extracted_studies.tsv')
  compiled_coloc_clustered_results_file <- glue::glue('{extracted_study_dir}/compiled_coloc_clustered_results.tsv')
  
  finemapped_studies_files <- Filter(
    function(file) file.exists(file), 
    glue::glue('{ld_block_dirs}/finemapped_studies.tsv')
  )
  
  if (length(finemapped_studies_files) > 0) {
    flog.info(paste("Processing", length(finemapped_studies_files), "finemapped study files for", gwas_info$metadata$guid))
    tryCatch({
      study_extractions <- vroom::vroom(finemapped_studies_files, show_col_types = FALSE, col_types = finemapped_column_types) |>
        dplyr::filter(min_p <= gwas_info$metadata$p_value_threshold) |>
        dplyr::left_join(snp_annotations, by = "snp") |>
        dplyr::select(study, chr, bp, snp, rsid, display_snp, unique_study_id, file, min_p, ld_block)
      flog.info(paste("Successfully processed", nrow(study_extractions), "study extractions"))
    }, error = function(e) {
      flog.error(paste("Error processing finemapped studies files:", e$message))
      flog.error(paste("Files:", paste(finemapped_studies_files, collapse = ", ")))
      stop(e)
    })
  } else {
    flog.warn("No finemapped study files found")
    study_extractions <- data.frame()
  }
  
  coloc_pairwise_results_files <- Filter(
    function(file) file.exists(file), 
    glue::glue('{ld_block_dirs}/coloc_pairwise_results.tsv.gz')
  )
  flog.info(glue::glue('{ld_block_dirs}'))
  flog.info(coloc_pairwise_results_files)
  if (length(coloc_pairwise_results_files) > 0) {
    flog.info(paste("Processing", length(coloc_pairwise_results_files), "coloc pairwise result files for", gwas_info$metadata$guid))
    coloc_results <- vroom::vroom(coloc_pairwise_results_files, delim = '\t', show_col_types = FALSE) |>
      dplyr::filter(study_a == gwas_info$metadata$guid | study_b == gwas_info$metadata$guid)
  } else {
    coloc_results <- data.frame()
  }

  coloc_clustered_results_files <- Filter(
    function(file) file.exists(file), 
    glue::glue('{ld_block_dirs}/coloc_clustered_results.tsv.gz')
  )
  if (length(coloc_clustered_results_files) > 0) {
    flog.info(paste("Processing", length(coloc_clustered_results_files), "coloc clustered result files for", gwas_info$metadata$guid))
    coloc_clustered_results <- vroom::vroom(coloc_clustered_results_files, delim = '\t', show_col_types = FALSE)
  } else {
    coloc_clustered_results <- data.frame()
  }
  
  flog.info(paste("Writing compiled results to files for", gwas_info$metadata$guid))
  vroom::vroom_write(coloc_results, compiled_coloc_results_file)
  vroom::vroom_write(study_extractions, compiled_study_extractions_file)
  vroom::vroom_write(coloc_clustered_results, compiled_coloc_clustered_results_file)
  
  return(list(
    coloc_results = coloc_results,
    study_extractions = study_extractions,
    coloc_clustered_results = coloc_clustered_results
  ))
}


send_update_gwas_upload <- function(gwas_info, success, failure_reason, results = NULL) {
  api_url <- glue::glue("https://gpmap.opengwas.io/api/v1/gwas/{gwas_info$metadata$guid}")
  flog.info(paste("Sending update to API:", api_url))

  if (success) {
    put_body <- list(
      success = success,
      results = results
    )
  } else {  
    put_body <- list(
      success = success,
      failure_reason = failure_reason,
      results = results
    )
  }
  if (!is.na(TEST_RUN)) return()

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
    
  flog.info(paste("Successfully updated GWAS results for", gwas_info$metadata$guid))
}

main()
