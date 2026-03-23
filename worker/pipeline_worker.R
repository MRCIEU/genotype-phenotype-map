setwd("pipeline_steps")
source("../worker/redis_client.R")
source("../worker/compile_pipeline_results.R")
source("constants.R")
source("gwas_calculations.R")
source("common_extraction_functions.R")

library(futile.logger)

log_dir <- file.path(data_dir, "logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(log_dir, paste0("pipeline_worker_", format(Sys.time(), "%Y%m"), ".log"))

flog.appender(appender.tee(log_file))
flog.threshold(INFO)

parser <- argparser::arg_parser("Pipeline worker")
parser <- argparser::add_argument(
  parser,
  "--custom_message_file",
  help = "Custom message to process (if testing)",
  type = "character",
  default = NULL
)
args <- argparser::parse_args(parser)
parallel_block_processing <- 6

if (!is_test_run) {
  tryCatch(
    {
      redis_conn <- connect_to_redis()
    },
    error = function(e) {
      flog.error(paste("Critical: Failed to establish Redis connection. Exiting. Error:", e$message))
      q(status = 1, save = "no")
    }
  )
} else {
  flog.appender(appender.console())
  redis_conn <- list(
    BRPOP = function(queue, timeout = 0.1) {
      if (queue == delete_gwas_queue) return(NULL)
      raw_payload <- paste(readLines(args$custom_message_file, warn = FALSE), collapse = "\n")
      parsed_payload <- tryCatch(jsonlite::fromJSON(raw_payload), error = function(e) NULL)
      if (!is.character(parsed_payload) || length(parsed_payload) < 2) {
        return(if (queue == process_gwas) list(process_gwas, raw_payload) else NULL)
      }
      file_queue_name <- parsed_payload[[1]]
      payload <- parsed_payload[[2]]
      if (queue == file_queue_name) return(list(queue, payload))
      return(NULL)
    },
    LPUSH = function(queue, message) {
      return(1)
    },
    LREM = function(queue, count, element) {
      return(0)
    }
  )
}

main <- function() {
  flog.info("Starting pipeline worker")

  while (TRUE) {
    # If the stop_processing file exists, hang forever
    stop_processing_file <- file.exists(glue::glue("{data_dir}/stop_processing"))
    if (stop_processing_file) {
      flog.info("Stop processing file found, hanging forever")
      system("tail -f /dev/null", wait = T)
    }

    if (!exists("redis_conn") || is.null(redis_conn)) {
      flog.error("Critical: redis_conn is not available. Exiting main loop.")
      break
    }

    tryCatch(
      {
        print("Getting delete message")
        delete_message <- get_from_delete_queue(redis_conn)
        print(delete_message)

        if (!is.null(delete_message)) {
          delete_info <- jsonlite::fromJSON(delete_message[[2]])
          flog.info(paste(delete_info$guid, "Received new message from delete queue"))
          delete_gwas(delete_info$guid)
          next
        }
        print("Getting process message")
        redis_message <- get_from_in_progress_queue(redis_conn)
        if (is.null(redis_message)) {
          redis_message <- get_from_process_queue(redis_conn)
        }
        print(redis_message)

        if (!is.null(redis_message)) {
          payload <- redis_message[[2]]
          add_to_in_progress_queue(redis_conn, payload)
          tryCatch(
            {
              flog.info(paste("Processing message with payload: ", payload))
              gwas_info <- jsonlite::fromJSON(payload)
              flog.info(paste(gwas_info$metadata$guid, "Received new message from queue"))

              processed <- process_message(gwas_info, payload)

              if (!is.null(processed)) {
                flog.error(paste(gwas_info$metadata$guid, "Failed to process message"))

                message_and_error <- list(
                  original_message = payload,
                  error = processed,
                  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
                )
                send_to_dlq(redis_conn, jsonlite::toJSON(message_and_error, auto_unbox = TRUE))
              } else {
                flog.info(paste(gwas_info$metadata$guid, "Successfully processed message"))
              }
            },
            finally = {
              remove_from_in_progress_queue(redis_conn, payload)
            }
          )
        }
        Sys.sleep(5)
      },
      error = function(e) {
        if (exists("gwas_info")) {
          flog.error(paste(gwas_info$metadata$guid, "Error processing message:", e$message))
          if (!is.null(e$call)) {
            flog.error(paste(gwas_info$metadata$guid, "Error call:", deparse(e$call)))
          }
        } else {
          flog.error(paste("Error processing message: ", e$message))
        }
        if (is_test_run) stop(e)
      }
    )
    if (is_test_run) break
  }
  return()
}

process_message <- function(original_gwas_info, original_payload = NULL) {
  return(tryCatch(
    {
      update_directories_for_worker(original_gwas_info$metadata$guid)

      if (dir.exists(extracted_study_dir)) {
        unlink(extracted_study_dir, recursive = TRUE)
        flog.info(paste(original_gwas_info$metadata$guid, "Cleaned up previous extracted_study_dir"))
      }
      if (dir.exists(ld_block_data_dir)) {
        unlink(ld_block_data_dir, recursive = TRUE)
        flog.info(paste(original_gwas_info$metadata$guid, "Cleaned up previous ld_block_data_dir"))
      }
      if (dir.exists(pipeline_metadata_dir)) {
        unlink(pipeline_metadata_dir, recursive = TRUE)
        flog.info(paste(original_gwas_info$metadata$guid, "Cleaned up previous pipeline_metadata_dir"))
      }

      dir.create(pipeline_metadata_dir, recursive = T, showWarnings = F)
      dir.create(extracted_study_dir, recursive = T, showWarnings = F)
      dir.create(ld_block_data_dir, recursive = T, showWarnings = F)

      gwas_data <- get_gwas_data_from_oracle(original_gwas_info)
      gwas_info <- gwas_data$gwas_info
      gwas <- gwas_data$gwas

      create_study_metadata_files(gwas_info)
      verification_result <- verify_gwas_data(gwas_info, gwas)
      if (!verification_result$valid) {
        flog.error(paste(gwas_info$metadata$guid, verification_result$error))
        send_update_gwas_upload(gwas_info, FALSE, paste("Caught error: ", verification_result$error))
        return()
      }

      updated_gwas <- change_column_names(gwas, gwas_info$metadata$column_names)
      if (!"EAF" %in% colnames(updated_gwas)) {
        gwas$EAF <- NA
        vroom::vroom_write(gwas, gwas_info$metadata$file_location)
        flog.info(paste(gwas_info$metadata$guid, "Added EAF column (NA) for LD panel fill-in during standardisation"))
      }

      flog.info(paste(gwas_info$metadata$guid, "Extracting regions"))
      extract_regions <- glue::glue(
        "Rscript extract_regions_from_summary_stats.R",
        " --worker_guid {gwas_info$metadata$guid} 2>&1"
      )
      output <- system(extract_regions, wait = T, intern = T)
      check_pipeline_step_complete(glue::glue("{extracted_study_dir}/extracted_snps.tsv"), "extracted_snps.tsv", output)

      if (file.info(glue::glue("{extracted_study_dir}/extracted_snps.tsv"))$size == 0) {
        flog.error(paste(gwas_info$metadata$guid, "No regions extracted"))
        stop(paste(
          "Caught error: No regions extracted from GWAS,",
          "please ensure the p-value threshold is set correctly"
        ))
      }

      flog.info(paste(gwas_info$metadata$guid, "Organising LD blocks"))
      ld_blocks_to_colocalise_file <- glue::glue("{pipeline_metadata_dir}/updated_ld_blocks_to_colocalise.tsv")
      organise_ld_blocks <- glue::glue(
        "Rscript organise_extracted_regions_into_ld_blocks.R",
        " --output_file {ld_blocks_to_colocalise_file}",
        " --worker_guid {gwas_info$metadata$guid} 2>&1"
      )
      output <- system(organise_ld_blocks, wait = T, intern = T)
      check_pipeline_step_complete(ld_blocks_to_colocalise_file, "updated_ld_blocks_to_colocalise.tsv", output)

      # cleaning up memory, so it utilises less in the parallel processes
      rm(gwas, updated_gwas, gwas_data)
      gc()

      ld_blocks_to_colocalise <- vroom::vroom(ld_blocks_to_colocalise_file, show_col_types = F)
      memory_intensive_blocks <- identify_memory_intensive_blocks(ld_blocks_to_colocalise$ld_block)

      blocks_parallel <- ld_blocks_to_colocalise$ld_block[
        !ld_blocks_to_colocalise$ld_block %in% memory_intensive_blocks
      ]
      blocks_sequential <- ld_blocks_to_colocalise$ld_block[
        ld_blocks_to_colocalise$ld_block %in% memory_intensive_blocks
      ]

      flog.info(paste(gwas_info$metadata$guid, "Processing", length(blocks_parallel), "blocks in parallel"))
      processed_blocks_parallel <- parallel::mclapply(
        blocks_parallel,
        mc.cores = parallel_block_processing,
        function(block) {
          return(process_single_block(block, gwas_info))
        }
      )

      flog.info(paste(
        gwas_info$metadata$guid,
        "Processing",
        length(blocks_sequential),
        "memory-intensive blocks sequentially"
      ))
      processed_blocks_sequential <- lapply(blocks_sequential, function(block) {
        return(process_single_block(block, gwas_info))
      })
      processed_blocks <- c(processed_blocks_parallel, processed_blocks_sequential)

      successful_block_names <- unlist(processed_blocks[!sapply(processed_blocks, is.null)])
      failed_blocks <- setdiff(ld_blocks_to_colocalise$ld_block, successful_block_names)

      if (length(failed_blocks) > 0) {
        flog.error(paste(gwas_info$metadata$guid, "Failed blocks:", paste(failed_blocks, collapse = ", ")))
        stop(paste(gwas_info$metadata$guid, "Failed to process all blocks"))
      }

      gc()
      flog.info(paste(gwas_info$metadata$guid, "Compiling results"))
      results <- compile_results(gwas_info)

      if (!is_test_run) {
        flog.info(paste(gwas_info$metadata$guid, "Uploading results"))
        upload_results(results, gwas_info)
        send_update_gwas_upload(gwas_info, TRUE, NULL, results)
      }
    },
    error = function(e) {
      error_msg <- e$message

      flog.error(paste(gwas_info$metadata$guid, "Error processing message:", error_msg))
      flog.error(paste(gwas_info$metadata$guid, "Error class:", class(e)[1]))
      if (!is.null(e$call)) {
        flog.error(paste(gwas_info$metadata$guid, "Error call:", deparse(e$call)))
      }

      # Create error details
      message_and_error <- list(
        original_message = original_payload,
        error = error_msg,
        timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      )

      send_to_dlq(redis_conn, jsonlite::toJSON(message_and_error, auto_unbox = TRUE))
      send_update_gwas_upload(original_gwas_info, FALSE, error_msg)

      return(error_msg)
    }
  ))
}


verify_gwas_data <- function(gwas_info, gwas) {
  gwas_info$metadata$column_names <- Filter(\(column) {
    return(!is.null(column) && length(column) > 0 &&
             (is.character(column) && nchar(column) > 0 || !is.character(column)))
  }, gwas_info$metadata$column_names)

  gwas <- change_column_names(gwas, gwas_info$metadata$column_names)
  mandatory_columns <- c("CHR", "BP", "P", "EA", "OA")

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

  error_checks <- ""

  if (any(is.na(gwas$CHR))) {
    error_checks <- paste(error_checks, "Some SNPs have missing CHR values, ")
  }

  if (any(is.na(gwas$BP))) {
    error_checks <- paste(error_checks, "Some SNPs have missing BP values, ")
  }

  if (any(gwas$P < 0, na.rm = T)) {
    error_checks <- paste(error_checks, "Some SNPs have negative P-values, ")
  }

  if ("EAF" %in% colnames(gwas) && any(gwas$EAF < 0 | gwas$EAF > 1, na.rm = T)) {
    error_checks <- paste(error_checks, "Some SNPs have EAF values out of range, ")
  }

  # Fail only when mix of present and missing EAF - all missing is OK (filled from LD panel)
  if ("EAF" %in% colnames(gwas) && !all(is.na(gwas$EAF)) && any(is.na(gwas$EAF))) {
    error_checks <- paste(error_checks, "Some SNPs have missing EAF values, ")
  }

  if (has_beta && any(gwas$SE < 0, na.rm = T)) {
    error_checks <- paste(error_checks, "Some SNPs have negative SE values, ")
  }

  if (nchar(error_checks) > 0) {
    return(list(valid = FALSE, error = error_checks))
  }

  return(list(valid = TRUE, error = NULL))
}

create_study_metadata_files <- function(gwas_info) {
  file.create(glue::glue("{extracted_study_dir}/study_metadata.json"))
  jsonlite::write_json(gwas_info$metadata,
    glue::glue("{extracted_study_dir}/study_metadata.json"),
    auto_unbox = TRUE,
    pretty = TRUE
  )

  study_to_process <- data.frame(
    data_type = data_types$phenotype,
    data_format = gwas_info$metadata$file_type,
    source = "user",
    study_name = gwas_info$metadata$guid,
    trait = gwas_info$metadata$name,
    ancestry = gwas_info$metadata$ancestry,
    sample_size = gwas_info$metadata$sample_size,
    category = gwas_info$metadata$category,
    study_location = gwas_info$metadata$file_location,
    extracted_location = extracted_study_dir,
    reference_build = gwas_info$metadata$reference_build,
    p_value_threshold = gwas_info$metadata$p_value_threshold,
    variant_type = variant_types$common,
    gene = NA,
    probe = NA,
    tissue = NA,
    coverage = coverage_types$dense
  )
  studies_to_process_file <- glue::glue("{pipeline_metadata_dir}/studies_to_process.tsv")
  vroom::vroom_write(study_to_process, studies_to_process_file)
  return()
}

identify_memory_intensive_blocks <- function(blocks, threshold = 10000, guid = NA) {
  memory_intensive_blocks <- c()

  for (block in blocks) {
    finemapped_file <- glue::glue("{data_dir}/ld_blocks/{block}/finemapped_studies.tsv")

    if (!file.exists(finemapped_file)) next

    dt <- data.table::fread(finemapped_file, select = "min_p", nThread = 1, verbose = FALSE, showProgress = FALSE)
    dt <- dt[min_p < min_p_allowed_for_worker]
    num_studies <- nrow(dt)

    if (num_studies > threshold) memory_intensive_blocks <- c(memory_intensive_blocks, block)
  }

  return(memory_intensive_blocks)
}

process_single_block <- function(block, gwas_info) {
  return(tryCatch(
    {
      start_time <- Sys.time()
      ld_info <- ld_block_dirs(block)

      output_files <- list(
        standardised = glue::glue("{ld_info$ld_block_data}/standardisation_complete"),
        imputed = glue::glue("{ld_info$ld_block_data}/imputation_complete"),
        finemapped = glue::glue("{ld_info$ld_block_data}/finemapping_complete"),
        coloc = glue::glue("{ld_info$ld_block_data}/coloc_complete")
      )

      flog.info(paste(gwas_info$metadata$guid, "Standardising regions for block:", block))
      standardise_regions <- glue::glue(
        "Rscript standardise_studies_in_ld_block.R",
        " --ld_block {block} ",
        " --completed_output_file {output_files$standardised}",
        " --worker_guid {gwas_info$metadata$guid} 2>&1"
      )
      output <- system(standardise_regions, wait = T, intern = T)
      check_pipeline_step_complete(output_files$standardised, block, output)
      gc()

      flog.info(paste(gwas_info$metadata$guid, "Imputing regions for block:", block))
      impute_regions <- glue::glue(
        "Rscript impute_studies_in_ld_block.R",
        " --ld_block {block} ",
        " --completed_output_file {output_files$imputed}",
        " --worker_guid {gwas_info$metadata$guid} 2>&1"
      )
      output <- system(impute_regions, wait = T, intern = T)
      check_pipeline_step_complete(output_files$imputed, block, output)
      gc()

      flog.info(paste(gwas_info$metadata$guid, "Finemapping regions for block:", block))
      finemap_regions <- glue::glue(
        "Rscript finemap_studies_in_ld_block.R",
        " --ld_block {block} ",
        " --completed_output_file {output_files$finemapped}",
        " --worker_guid {gwas_info$metadata$guid} 2>&1"
      )
      output <- system(finemap_regions, wait = T, intern = T)
      check_pipeline_step_complete(output_files$finemapped, block, output)
      gc()

      flog.info(paste(gwas_info$metadata$guid, "Colocalising regions for block:", block))
      gwas_upload_ids_to_compare <- gwas_info$metadata$compare_with_upload_guids
      if (is.null(gwas_upload_ids_to_compare)) gwas_upload_ids_to_compare <- character(0)
      gwas_upload_ids_to_compare <- as.character(unlist(gwas_upload_ids_to_compare))

      if (length(gwas_upload_ids_to_compare) > 0) {
        compare_val <- paste(gwas_upload_ids_to_compare, collapse = ",")
        compare_ids_arg <- glue::glue(' --gwas_upload_ids_to_compare {shQuote(compare_val, type = "sh")}')
      } else {
        compare_ids_arg <- ""
      }
      coloc_regions <- glue::glue(
        "Rscript coloc_and_cluster_studies_in_ld_block.R",
        " --ld_block {block} ",
        " --completed_output_file {output_files$coloc}",
        " --worker_guid {gwas_info$metadata$guid}",
        " --worker_p_value_threshold {gwas_info$metadata$p_value_threshold}",
        "{compare_ids_arg}",
        " 2>&1"
      )
      output <- system(coloc_regions, wait = T, intern = T)
      check_pipeline_step_complete(output_files$coloc, block, output)

      flog.info(paste(gwas_info$metadata$guid, "Time taken for block:", block, diff_time_taken(start_time)))
      return(block)
    },
    error = function(e) {
      flog.error(paste(gwas_info$metadata$guid, "Error processing block:", block, "-", e$message))
      return(NULL)
    }
  ))
}

check_pipeline_step_complete <- function(output_file, ld_block, output) {
  cmd_failed <- !is.null(attr(output, "status")) && attr(output, "status") != 0
  file_missing <- !file.exists(output_file)

  if (cmd_failed || file_missing) {
    fail_file <- glue::glue("{output_file}_failed.txt")
    writeLines(output, fail_file)
    error_msg <- glue::glue("{ld_block}: FAILED {output_file}")
    flog.error(error_msg)
    stop(error_msg)
  }
}

get_gwas_data_from_oracle <- function(gwas_info) {
  if (is_test_run) {
    gwas_info$metadata$file_location <- gwas_info$file_location
    gwas_info$metadata$file_type <- extraction_file_types$csv
    gwas <- vroom::vroom(gwas_info$file_location, show_col_types = F)
    return(list(gwas_info = gwas_info, gwas = gwas))
  }

  local_file_path <- file.path(extracted_study_dir, basename(gwas_info$file_location))
  dir.create(dirname(local_file_path), recursive = TRUE, showWarnings = FALSE)

  download_cmd <- paste(
    "oci os object get",
    "--auth", "instance_principal",
    "--bucket-name", shQuote(oracle_bucket_name),
    "--name", shQuote(gwas_info$file_location),
    "--file", shQuote(local_file_path)
  )

  download_status <- system(download_cmd, wait = TRUE)

  if (download_status != 0 || !file.exists(local_file_path)) {
    error_msg <- paste(gwas_info$metadata$guid, "Failed to download file from Oracle bucket:", gwas_info$file_location)
    flog.error(error_msg)
    send_update_gwas_upload(gwas_info, FALSE, error_msg)
    stop(error_msg)
  }

  gwas_info$metadata$file_location <- local_file_path
  gwas_info$file_location <- local_file_path

  if (grepl(".vcf", local_file_path, ignore.case = TRUE)) {
    # TODO: Implement VCF support
    stop("Caught error: VCF support not yet implemented")
  } else {
    gwas_info$metadata$file_type <- extraction_file_types$csv
    gwas <- vroom::vroom(local_file_path, show_col_types = F)
  }

  return(list(gwas_info = gwas_info, gwas = gwas))
}

upload_results <- function(results, gwas_info) {
  bucket_prefix <- glue::glue("gwas_upload/{gwas_info$metadata$guid}/")

  cmd <- paste(
    "oci os object sync",
    "--auth", "instance_principal",
    "--bucket-name", shQuote(oracle_bucket_name),
    "--src-dir", shQuote(extracted_study_dir),
    "--prefix", shQuote(bucket_prefix)
  )

  status <- system(cmd, wait = TRUE)

  if (status != 0) {
    error_msg <- paste(gwas_info$metadata$guid, "Failed to upload results to Oracle bucket")
    flog.error(error_msg)
    stop(error_msg)
  }

  flog.info(paste(gwas_info$metadata$guid, "Successfully uploaded results to Oracle bucket"))
  return()
}

send_update_gwas_upload <- function(gwas_info, success, failure_reason, results = NULL) {
  if (is_test_run) {
    return()
  }

  api_url <- glue::glue("https://gpmap.opengwas.io/api/v1/gwas/{gwas_info$metadata$guid}")
  flog.info(paste(gwas_info$metadata$guid, "Sending update to API:", api_url))

  if (success) {
    put_body <- list(
      success = success,
      coloc_pairs = results$coloc_pairwise_results,
      study_extractions = results$study_extractions,
      coloc_groups = results$coloc_clustered_results,
      associations = results$associations
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
    httr::add_headers("Content-Type" = "application/json"),
    httr::timeout(1200)
  )
  if (httr::status_code(response) != 200) {
    error_msg <- paste(gwas_info$metadata$guid, "Error updating GWAS:", httr::content(response, "text"))
    flog.error(error_msg)
    stop(error_msg)
  }
}

delete_gwas <- function(guid) {
  flog.info(paste(guid, "Deleting GWAS"))
  all_delete_status <- c()
  delete_status <- system(paste0("rm -rf ", gwas_upload_dir, "gwas_upload/", guid), wait = TRUE)
  all_delete_status <- c(all_delete_status, delete_status)
  delete_status <- system(paste0("rm -rf ", gwas_upload_dir, "ld_blocks/gwas_upload/", guid), wait = TRUE)
  all_delete_status <- c(all_delete_status, delete_status)
  delete_status <- system(paste0("rm -rf ", gwas_upload_dir, "pipeline_metadata/gwas_upload/", guid), wait = TRUE)
  all_delete_status <- c(all_delete_status, delete_status)

  if (any(all_delete_status != 0)) {
    error_msg <- paste(guid, "Failed to delete GWAS")
    flog.error(error_msg)
  }
  return()
}

main()
