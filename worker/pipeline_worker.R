setwd("pipeline_steps")
source("../worker/redis_client.R")
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

              processed <- process_message(gwas_info)

              if (!is.null(processed)) {
                flog.error(paste(gwas_info$metadata$guid, "Failed to process message"))

                message_and_error <- list(
                  original_message = gwas_info,
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

process_message <- function(original_gwas_info) {
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

      gwas_data <- get_gwas_data(original_gwas_info)
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
      check_step_complete(glue::glue("{extracted_study_dir}/extracted_snps.tsv"), "extracted_snps.tsv", output)

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
      check_step_complete(ld_blocks_to_colocalise_file, "updated_ld_blocks_to_colocalise.tsv", output)

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
      error_details <- list(
        original_message = original_gwas_info,
        error = error_msg,
        timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      )

      send_to_dlq(redis_conn, jsonlite::toJSON(error_details, auto_unbox = TRUE))

      send_update_gwas_upload(original_gwas_info, FALSE, error_msg)

      return(error_msg)
    }
  ))
}

get_gwas_data <- function(gwas_info) {
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

verify_gwas_data <- function(gwas_info, gwas) {
  gwas_info$metadata$column_names <- Filter(\(column) {
    return(!is.null(column) && length(column) > 0 &&
             (is.character(column) && nchar(column) > 0 || !is.character(column)))
  }, gwas_info$metadata$column_names)

  gwas <- change_column_names(gwas, gwas_info$metadata$column_names)
  mandatory_columns <- c("CHR", "BP", "P", "EA", "OA")
  beta_columns <- c("BETA", "SE")
  or_columns <- c("OR", "OR_LB", "OR_UB")

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

  if ("EAF" %in% colnames(gwas) && any(is.na(gwas$EAF))) {
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
      check_step_complete(output_files$standardised, block, output)
      gc()

      flog.info(paste(gwas_info$metadata$guid, "Imputing regions for block:", block))
      impute_regions <- glue::glue(
        "Rscript impute_studies_in_ld_block.R",
        " --ld_block {block} ",
        " --completed_output_file {output_files$imputed}",
        " --worker_guid {gwas_info$metadata$guid} 2>&1"
      )
      output <- system(impute_regions, wait = T, intern = T)
      check_step_complete(output_files$imputed, block, output)
      gc()

      flog.info(paste(gwas_info$metadata$guid, "Finemapping regions for block:", block))
      finemap_regions <- glue::glue(
        "Rscript finemap_studies_in_ld_block.R",
        " --ld_block {block} ",
        " --completed_output_file {output_files$finemapped}",
        " --worker_guid {gwas_info$metadata$guid} 2>&1"
      )
      output <- system(finemap_regions, wait = T, intern = T)
      check_step_complete(output_files$finemapped, block, output)
      gc()

      flog.info(paste(gwas_info$metadata$guid, "Colocalising regions for block:", block))
      gwas_upload_ids_to_compare <- gwas_info$metadata$gwas_upload_ids_to_compare
      if (is.null(gwas_upload_ids_to_compare)) gwas_upload_ids_to_compare <- character(0)
      gwas_upload_ids_to_compare <- as.character(unlist(gwas_upload_ids_to_compare))
      gwas_upload_ids_to_compare <- gwas_upload_ids_to_compare[nchar(trimws(gwas_upload_ids_to_compare)) > 0]
      compare_ids_arg <- if (length(gwas_upload_ids_to_compare) > 0) {
        compare_val <- paste(gwas_upload_ids_to_compare, collapse = ",")
        glue::glue(' --gwas_upload_ids_to_compare {shQuote(compare_val, type = "sh")}')
      } else {
        ""
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
      check_step_complete(output_files$coloc, block, output)
      flog.info(paste(gwas_info$metadata$guid, "Coloc output:", paste(output, collapse = "\n")))

      flog.info(paste(gwas_info$metadata$guid, "Time taken for block:", block, diff_time_taken(start_time)))
      return(block)
    },
    error = function(e) {
      flog.error(paste(gwas_info$metadata$guid, "Error processing block:", block, "-", e$message))
      return(NULL)
    }
  ))
}

check_step_complete <- function(output_file, ld_block, output) {
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

compile_results <- function(gwas_info) {
  ld_block_dirs <- list.dirs(ld_block_data_dir, recursive = TRUE, full.names = TRUE) |>
    (\(dirs) dirs[!dirs %in% dirname(dirs[-1])])()

  compiled_coloc_pairwise_results_file <- glue::glue("{extracted_study_dir}/compiled_coloc_pairwise_results.tsv")
  compiled_study_extractions_file <- glue::glue("{extracted_study_dir}/compiled_extracted_studies.tsv")
  compiled_coloc_clustered_results_file <- glue::glue("{extracted_study_dir}/compiled_coloc_clustered_results.tsv")
  compiled_associations_file <- glue::glue("{extracted_study_dir}/compiled_associations.tsv")
  lbfs_concatenated_file <- glue::glue("{extracted_study_dir}/gwas_with_lbfs.tsv.gz")

  finemapped_studies_files <- Filter(
    function(file) file.exists(file),
    glue::glue("{ld_block_dirs}/finemapped_studies.tsv")
  )
  compare_guids <- gwas_info$metadata$gwas_upload_ids_to_compare
  if (is.null(compare_guids)) compare_guids <- character(0)
  compare_guids <- as.character(unlist(compare_guids))
  compare_guids <- setdiff(compare_guids, gwas_info$metadata$guid)
  if (length(compare_guids) > 0) {
    ld_blocks <- sub(paste0("^", ld_block_data_dir), "", ld_block_dirs)
    compare_files <- as.character(outer(
      compare_guids,
      ld_blocks,
      function(guid, block) glue::glue("{gwas_upload_dir}ld_blocks/gwas_upload/{guid}/{block}/finemapped_studies.tsv")
    ))
    finemapped_studies_files <- c(finemapped_studies_files, Filter(file.exists, compare_files))
  }

  if (length(finemapped_studies_files) > 0) {
    all_finemapped_studies <- lapply(finemapped_studies_files, function(file) {
      se_result <- data.table::fread(
        file,
        showProgress = FALSE,
        colClasses = list(character = c("chr", "snp", "study", "unique_study_id"))
      ) |>
        dplyr::filter(min_p <= gwas_info$metadata$p_value_threshold)
      return(se_result)
    }) |>
      data.table::rbindlist(fill = TRUE)
    study_extractions <- all_finemapped_studies |> dplyr::filter(study == gwas_info$metadata$guid)
  } else {
    flog.warn(paste(gwas_info$metadata$guid, "No finemapped study files found"))
    study_extractions <- data.table::data.table()
    all_finemapped_studies <- data.table::data.table()
  }

  all_snps_in_ld_blocks <- study_extractions |>
    dplyr::distinct(snp) |>
    dplyr::pull(snp)

  snp_annotations <- data.table::fread(
    file.path(variant_annotation_dir, "vep_annotations_hg38.tsv.gz"),
    select = c("snp", "rsid", "display_snp"),
    showProgress = FALSE,
    colClasses = list(character = c("snp", "rsid", "display_snp"))
  ) |>
    dplyr::mutate(snp = trimws(snp)) |>
    dplyr::filter(snp %in% all_snps_in_ld_blocks)

  if (nrow(study_extractions) > 0) {
    study_extractions <- merge(study_extractions, snp_annotations, by = "snp", all.x = TRUE)
  }

  study_extractions <- study_extractions |>
    dplyr::select(study, unique_study_id, snp, chr, bp, min_p, ld_block, file, file_with_lbfs)
  vroom::vroom_write(study_extractions, compiled_study_extractions_file)

  coloc_clustered_results_files <- Filter(
    function(file) file.exists(file),
    glue::glue("{ld_block_dirs}/coloc_clustered_results.tsv.gz")
  )
  if (length(coloc_clustered_results_files) > 0) {
    coloc_clustered_results <- lapply(coloc_clustered_results_files, function(file) {
      return(data.table::fread(file, showProgress = FALSE))
    }) |>
      data.table::rbindlist(fill = TRUE)

    if (nrow(coloc_clustered_results) > 0) {
      worker_unique_study_ids <- study_extractions$unique_study_id
      components_with_worker <- coloc_clustered_results |>
        dplyr::filter(unique_study_id %in% worker_unique_study_ids) |>
        dplyr::distinct(ld_block, component)
      coloc_clustered_results <- coloc_clustered_results |>
        dplyr::inner_join(components_with_worker, by = c("ld_block", "component")) |>
        dplyr::group_by(ld_block, component) |>
        dplyr::mutate(coloc_group_id = dplyr::cur_group_id()) |>
        dplyr::ungroup() |>
        dplyr::arrange(coloc_group_id) |>
        dplyr::select(unique_study_id, coloc_group_id, snp, ld_block, h4_connectedness, h3_connectedness) |>
        data.table::as.data.table()
    }
  } else {
    coloc_clustered_results <- data.table::data.table()
  }
  vroom::vroom_write(coloc_clustered_results, compiled_coloc_clustered_results_file)

  coloc_pairwise_results_files <- Filter(
    function(file) file.exists(file),
    glue::glue("{ld_block_dirs}/coloc_pairwise_results.tsv.gz")
  )
  if (length(compare_guids) > 0) {
    ld_blocks <- sub(paste0("^", ld_block_data_dir), "", ld_block_dirs)
    compare_pairwise_files <- as.character(outer(
      compare_guids,
      ld_blocks,
      function(guid, block) {
        return(paste0(gwas_upload_dir, "ld_blocks/gwas_upload/", guid, "/", block, "/coloc_pairwise_results.tsv.gz"))
      }
    ))
    coloc_pairwise_results_files <- c(coloc_pairwise_results_files, Filter(file.exists, compare_pairwise_files))
  }

  if (length(coloc_pairwise_results_files) > 0) {
    coloc_pairwise_results <- lapply(coloc_pairwise_results_files, function(file) {
      cp_result <- data.table::fread(file, showProgress = FALSE) |>
        dplyr::filter(study_a == gwas_info$metadata$guid | study_b == gwas_info$metadata$guid) |>
        dplyr::select(
          unique_study_a,
          unique_study_b,
          PP.H3.abf,
          PP.H4.abf,
          ld_block,
          false_positive,
          false_negative,
          ignore
        )
      return(cp_result)
    }) |>
      data.table::rbindlist(fill = TRUE)

    if (nrow(coloc_pairwise_results) > 0) {
      data.table::setnames(coloc_pairwise_results,
        old = c("unique_study_a", "unique_study_b", "PP.H3.abf", "PP.H4.abf"),
        new = c("unique_study_id_a", "unique_study_id_b", "h3", "h4")
      )
    }
  } else {
    coloc_pairwise_results <- data.table::data.table()
  }
  vroom::vroom_write(coloc_pairwise_results, compiled_coloc_pairwise_results_file)

  concatenate_file_with_lbfs(gwas_info, study_extractions)

  associations <- find_associations_for_coloc_clustered_snps(
    gwas_info,
    coloc_clustered_results,
    all_finemapped_studies,
    snp_annotations
  )

  vroom::vroom_write(associations, compiled_associations_file)

  return(list(
    coloc_pairwise_results = coloc_pairwise_results,
    study_extractions = study_extractions,
    coloc_clustered_results = coloc_clustered_results,
    associations = associations
  ))
}

find_associations_for_coloc_clustered_snps <- function(
  gwas_info,
  coloc_clustered_results,
  all_finemapped_studies,
  snp_annotations
) {
  flog.info(paste(gwas_info$metadata$guid, "Finding associations for coloc clustered SNPs"))
  if (nrow(coloc_clustered_results) == 0 || nrow(all_finemapped_studies) == 0) {
    return(data.frame())
  }

  finemapped_studies <- all_finemapped_studies |>
    dplyr::select(unique_study_id, study, file, ld_block) |>
    dplyr::filter(unique_study_id %in% coloc_clustered_results$unique_study_id)

  coloc_snps <- coloc_clustered_results |>
    dplyr::filter(unique_study_id %in% finemapped_studies$unique_study_id) |>
    dplyr::select(unique_study_id, snp, ld_block)

  if (nrow(coloc_snps) == 0 || nrow(finemapped_studies) == 0) {
    return(data.frame())
  }

  coloc_snps_with_files <- coloc_snps |>
    dplyr::left_join(finemapped_studies, by = c("unique_study_id", "ld_block"))

  files_to_read <- coloc_snps_with_files |>
    dplyr::group_by(file) |>
    dplyr::group_split()

  snp_mapping <- data.table::as.data.table(
    dplyr::distinct(coloc_snps_with_files, snp, study)
  )
  data.table::setkey(snp_mapping, snp)

  process_one_assoc_file <- function(file_group) {
    file_path <- file_group$file[1]
    study_for_file <- file_group$study[1]
    target_snps <- unique(file_group$snp)

    if (!grepl("^/", file_path) && !grepl("^study", file_path)) {
      file_path <- file.path(data_dir, file_path)
    } else if (grepl("^study", file_path)) {
      file_path <- file.path(data_dir, file_path)
    }

    if (!file.exists(file_path)) return(NULL)

    avail_cols <- names(data.table::fread(file_path, nrows = 0, showProgress = FALSE))
    required <- c("SNP", "BETA", "SE", "EAF")
    if (!all(required %in% avail_cols)) return(NULL)
    cols_to_read <- c(required, if ("IMPUTED" %in% avail_cols) "IMPUTED")

    gwas <- data.table::fread(
      file_path,
      select = cols_to_read,
      showProgress = FALSE,
      nThread = 1
    )

    gwas <- gwas[SNP %in% target_snps]
    if (nrow(gwas) == 0) return(NULL)

    gwas[, p := beta_se_to_p(BETA, SE)]
    if (!"IMPUTED" %in% colnames(gwas)) gwas[, IMPUTED := FALSE]
    snp_mapping_this_study <- snp_mapping[study == study_for_file]
    gwas <- gwas[snp_mapping_this_study, on = c(SNP = "snp"), nomatch = 0]

    if (nrow(gwas) == 0) return(NULL)
    return(gwas[, .(SNP, study, BETA, SE, p, EAF, IMPUTED)])
  }

  assoc_chunks <- parallel::mclapply(
    files_to_read,
    process_one_assoc_file,
    mc.cores = parallel_block_processing
  )
  valid_chunks <- Filter(Negate(is.null), assoc_chunks)
  if (length(valid_chunks) == 0) {
    return(data.frame())
  }
  associations <- data.table::rbindlist(valid_chunks)

  if (nrow(associations) == 0) {
    return(data.frame())
  }

  associations <- associations[, .(
    snp = SNP,
    study_name = study,
    beta = BETA,
    se = SE,
    p,
    eaf = EAF,
    imputed = IMPUTED
  )]

  flog.info(paste(
    gwas_info$metadata$guid,
    "Extracted",
    nrow(associations),
    "associations for coloc clustered SNPs"
  ))
  return(as.data.frame(associations))
}

concatenate_file_with_lbfs <- function(gwas_info, study_extractions) {
  lbfs_concatenated_file <- glue::glue("{extracted_study_dir}/gwas_with_lbfs.tsv.gz")
  if (nrow(study_extractions) > 0) {
    flog.info(paste(gwas_info$metadata$guid, "Concatenating file_with_lbfs files"))

    lbf_files <- unique(study_extractions$file_with_lbfs)
    lbf_files <- lbf_files[!is.na(lbf_files)]

    if (length(lbf_files) > 0) {
      if (file.exists(lbfs_concatenated_file)) unlink(lbfs_concatenated_file)

      for (i in seq_along(lbf_files)) {
        file_path <- glue::glue("{data_dir}/{lbf_files[i]}")
        if (!file.exists(file_path)) next

        dt <- data.table::fread(file_path, showProgress = FALSE)
        if (nrow(dt) > 0) {
          data.table::fwrite(dt, lbfs_concatenated_file,
            append = TRUE,
            compress = "gzip", sep = "\t"
          )
        }
        rm(dt)
        gc(verbose = FALSE)
      }
    }
  }
  return(NULL)
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
    httr::add_headers("Content-Type" = "application/json")
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
