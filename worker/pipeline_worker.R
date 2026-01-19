setwd('pipeline_steps')
source('../worker/redis_client.R')
source('constants.R')
source('common_extraction_functions.R')

library(futile.logger)

log_dir <- file.path(data_dir, "logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(log_dir, paste0("pipeline_worker_", format(Sys.time(), "%Y%m"), ".log"))

flog.appender(appender.tee(log_file))
flog.threshold(INFO)

parser <- argparser::arg_parser('Pipeline worker')
parser <- argparser::add_argument(parser, '--custom_message_file', help = 'Custom message to process (if testing)', type = 'character', default = NULL)
args <- argparser::parse_args(parser)

if (is.na(TEST_RUN)) {
  tryCatch({
    redis_conn <- connect_to_redis()
  }, error = function(e) {
    flog.error(paste("Critical: Failed to establish Redis connection. Exiting. Error:", e$message))
    q(status = 1, save = "no")
  })
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



main <- function() {
  flog.info("Starting pipeline worker")
  
  while(TRUE) {
    if (!exists("redis_conn") || is.null(redis_conn)) {
      flog.error("Critical: redis_conn is not available. Exiting main loop.")
      break
    }

    tryCatch({
      delete_message <- get_from_delete_queue(redis_conn)

      if (!is.null(delete_message)) {
        delete_info <- jsonlite::fromJSON(delete_message[[2]])
        flog.info(paste(delete_info$guid, "Received new message from delete queue"))
        delete_gwas(delete_info$guid)
        next 
      }

      redis_message <- get_from_process_queue(redis_conn)
      
      if (!is.null(redis_message)) {
        gwas_info <- jsonlite::fromJSON(redis_message[[2]])
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
      }
      Sys.sleep(5)
      stop_processing_file <- file.exists(glue::glue('{data_dir}/stop_processing'))
      if (!is.na(TEST_RUN) || stop_processing_file) {
        break
      }
    }, error = function(e) {
      if (exists("gwas_info")) {
        flog.error(paste(gwas_info$metadata$guid, "Error processing message:", e$message))
        if (!is.null(e$call)) {
          flog.error(paste(gwas_info$metadata$guid, "Error call:", deparse(e$call)))
        }
      } else {
        flog.error(paste("Error processing message: ", e$message))
      }
    })
  }
}

process_message <- function(original_gwas_info) {
  tryCatch({
    update_directories_for_worker(original_gwas_info$metadata$guid)

    dir.create(pipeline_metadata_dir, recursive = T, showWarnings = F)
    dir.create(extracted_study_dir, recursive = T, showWarnings = F)
    dir.create(ld_block_data_dir, recursive = T, showWarnings = F)

    gwas_data <- get_gwas_data(original_gwas_info)
    gwas_info <- gwas_data$gwas_info
    gwas <- gwas_data$gwas

    verification_result <- verify_gwas_data(gwas_info, gwas)
    if (!verification_result$valid) {
      flog.error(paste(gwas_info$metadata$guid, verification_result$error))
      send_update_gwas_upload(gwas_info, FALSE, paste("Caught error: ",verification_result$error))
      return()
    }

    flog.info(paste(gwas_info$metadata$guid, 'Extracting regions'))
    create_study_metadata_files(gwas_info)
    extract_regions <- glue::glue("Rscript extract_regions_from_summary_stats.R",
      " --worker_guid {gwas_info$metadata$guid} 2>&1")
    output <- system(extract_regions, wait = T, intern = T)
    check_step_complete(glue::glue('{extracted_study_dir}/extracted_snps.tsv'), 'extracted_snps.tsv', output)

    if (file.info(glue::glue('{extracted_study_dir}/extracted_snps.tsv'))$size == 0) {
      flog.error(paste(gwas_info$metadata$guid, 'No regions extracted'))
      stop(paste('Caught error: No regions extracted from GWAS, please ensure the p-value threshold is set correctly'))
    }

    flog.info(paste(gwas_info$metadata$guid, 'Organising LD blocks'))
    ld_blocks_to_colocalise_file <- glue::glue('{pipeline_metadata_dir}/updated_ld_blocks_to_colocalise.tsv')
    organise_ld_blocks <- glue::glue("Rscript organise_extracted_regions_into_ld_blocks.R",
      " --output_file {ld_blocks_to_colocalise_file}",
      " --worker_guid {gwas_info$metadata$guid} 2>&1")
    output <- system(organise_ld_blocks, wait = T, intern = T)
    check_step_complete(ld_blocks_to_colocalise_file, 'updated_ld_blocks_to_colocalise.tsv', output)

    #cleaning up memory, so it utilises less in the parallel processes
    rm(gwas)
    rm(gwas_data)
    gc()

    ld_blocks_to_colocalise <- vroom::vroom(ld_blocks_to_colocalise_file, show_col_types = F)
    memory_intensive_blocks <- identify_memory_intensive_blocks(ld_blocks_to_colocalise$ld_block)

    blocks_parallel <- ld_blocks_to_colocalise$ld_block[!ld_blocks_to_colocalise$ld_block %in% memory_intensive_blocks]
    blocks_sequential <- ld_blocks_to_colocalise$ld_block[ld_blocks_to_colocalise$ld_block %in% memory_intensive_blocks]

    flog.info(paste(gwas_info$metadata$guid, 'Processing', length(blocks_parallel), 'blocks in parallel'))
    parallel_block_processing <- 6
    processed_blocks_parallel <- parallel::mclapply(blocks_parallel, mc.cores = parallel_block_processing, function(block) {
      return(process_single_block(block, gwas_info))
    })

    flog.info(paste(gwas_info$metadata$guid, 'Processing', length(blocks_sequential), 'memory-intensive blocks sequentially'))
    processed_blocks_sequential <- lapply(blocks_sequential, function(block) {
      return(process_single_block(block, gwas_info))
    })
    processed_blocks <- c(processed_blocks_parallel, processed_blocks_sequential)

    successful_block_names <- unlist(processed_blocks[!sapply(processed_blocks, is.null)])
    failed_blocks <- setdiff(ld_blocks_to_colocalise$ld_block, successful_block_names)

    if (length(failed_blocks) > 0) {
      flog.error(paste(gwas_info$metadata$guid, 'Failed blocks:', paste(failed_blocks, collapse = ", ")))
      stop(paste(gwas_info$metadata$guid, 'Failed to process all blocks'))
    }

    flog.info(paste(gwas_info$metadata$guid, 'Compiling results'))
    results <- compile_results(gwas_info)

    if (is.na(TEST_RUN)) {
      flog.info(paste(gwas_info$metadata$guid, 'Uploading results'))
      upload_results(results, gwas_info)
      send_update_gwas_upload(gwas_info, TRUE, NULL, results)
    } 

  }, error = function(e) {
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
  })
}

get_gwas_data <- function(gwas_info) {
  if (!is.na(TEST_RUN)) {
    gwas_info$metadata$file_location <- gwas_info$file_location
    gwas_info$metadata$file_type <- extraction_file_types$csv
    gwas <- vroom::vroom(gwas_info$file_location, show_col_types = F)
    return(list(gwas_info = gwas_info, gwas = gwas))
  }

  local_file_path <- file.path(extracted_study_dir, basename(gwas_info$file_location))
  dir.create(dirname(local_file_path), recursive = TRUE, showWarnings = FALSE)

  download_cmd <- paste('oci os object get',
                        '--auth', 'instance_principal',
                        '--bucket-name', shQuote(oracle_bucket_name),
                        '--name', shQuote(gwas_info$file_location),
                        '--file', shQuote(local_file_path))

  download_status <- system(download_cmd, wait = TRUE)

  if (download_status != 0 || !file.exists(local_file_path)) {
    error_msg <- paste(gwas_info$metadata$guid, 'Failed to download file from Oracle bucket:', gwas_info$file_location)
    flog.error(error_msg)
    send_update_gwas_upload(gwas_info, FALSE, error_msg)
    stop(error_msg)
  }

  gwas_info$metadata$file_location <- local_file_path
  gwas_info$file_location <- local_file_path

  if (grepl('.vcf', local_file_path, ignore.case = TRUE)) {
    #TODO: Implement VCF support
    stop("Caught error: VCF support not yet implemented")
  } else {
    gwas_info$metadata$file_type <- extraction_file_types$csv
    gwas <- vroom::vroom(local_file_path, show_col_types = F)
  }

  return(list(gwas_info = gwas_info, gwas = gwas))
}

verify_gwas_data <- function(gwas_info, gwas) {
  gwas_info$metadata$column_names <- Filter(\(column) {
    !is.null(column) && length(column) > 0 &&
    (is.character(column) && nchar(column) > 0 || !is.character(column))
  }, gwas_info$metadata$column_names)

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
  file.create(glue::glue('{extracted_study_dir}/study_metadata.json'))
  jsonlite::write_json(gwas_info$metadata, 
                      glue::glue('{extracted_study_dir}/study_metadata.json'),
                      auto_unbox = TRUE,
                      pretty = TRUE)

  study_to_process <- data.frame(
    data_type = data_types$phenotype,
    data_format = gwas_info$metadata$file_type,
    source = 'user',
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
  studies_to_process_file <- glue::glue('{pipeline_metadata_dir}/studies_to_process.tsv')
  vroom::vroom_write(study_to_process, studies_to_process_file)
}

identify_memory_intensive_blocks <- function(blocks, threshold = 10000, guid = NA) {
  memory_intensive_blocks <- c()
  
  for (block in blocks) {
    finemapped_file <- glue::glue('{data_dir}/ld_blocks/{block}/finemapped_studies.tsv')

    if (!file.exists(finemapped_file)) next

    dt <- data.table::fread(finemapped_file, select = 'min_p', nThread = 1, verbose = FALSE, showProgress = FALSE)
    dt <- dt[min_p < min_p_allowed_for_worker]
    num_studies <- nrow(dt)
      
    if (num_studies > threshold) memory_intensive_blocks <- c(memory_intensive_blocks, block)
  }
  
  return(memory_intensive_blocks)
}

process_single_block <- function(block, gwas_info) {
  tryCatch({
    start_time <- Sys.time()
    ld_info <- ld_block_dirs(block)

    output_files <- list(
      standardised = glue::glue('{ld_info$ld_block_data}/standardisation_complete'),
      imputed = glue::glue('{ld_info$ld_block_data}/imputation_complete'),
      finemapped = glue::glue('{ld_info$ld_block_data}/finemapping_complete'),
      coloc = glue::glue('{ld_info$ld_block_data}/coloc_complete')
    )

    flog.info(paste(gwas_info$metadata$guid, 'Standardising regions for block:', block))
    standardise_regions <- glue::glue("Rscript standardise_studies_in_ld_block.R",
      " --ld_block {block} ", 
      " --completed_output_file {output_files$standardised}",
      " --worker_guid {gwas_info$metadata$guid} 2>&1")
    output <- system(standardise_regions, wait = T, intern = T)
    check_step_complete(output_files$standardised, block, output)
    gc()

    flog.info(paste(gwas_info$metadata$guid, 'Imputing regions for block:', block))
    impute_regions <- glue::glue("Rscript impute_studies_in_ld_block.R",
      " --ld_block {block} ",
      " --completed_output_file {output_files$imputed}",
      " --worker_guid {gwas_info$metadata$guid} 2>&1")
    output <- system(impute_regions, wait = T, intern = T)
    check_step_complete(output_files$imputed, block, output)
    gc()

    flog.info(paste(gwas_info$metadata$guid, 'Finemapping regions for block:', block))
    finemap_regions <- glue::glue("Rscript finemap_studies_in_ld_block.R",
      " --ld_block {block} ",
      " --completed_output_file {output_files$finemapped}",
      " --worker_guid {gwas_info$metadata$guid} 2>&1")
    output <- system(finemap_regions, wait = T, intern = T)
    check_step_complete(output_files$finemapped, block, output)
    gc()

    flog.info(paste(gwas_info$metadata$guid, 'Colocalising regions for block:', block))
    coloc_regions <- glue::glue("Rscript coloc_and_cluster_studies_in_ld_block.R",
      " --ld_block {block} ",
      " --completed_output_file {output_files$coloc}",
      " --worker_guid {gwas_info$metadata$guid}",
      " --worker_p_value_threshold {gwas_info$metadata$p_value_threshold} 2>&1")
    output <- system(coloc_regions, wait = T, intern = T)
    check_step_complete(output_files$coloc, block, output)

    flog.info(paste(gwas_info$metadata$guid, 'Time taken for block:', block, diff_time_taken(start_time)))
    return(block)
  }, error = function(e) {
    flog.error(paste(gwas_info$metadata$guid, 'Error processing block:', block, '-', e$message))
    return(NULL)
  })
}

check_step_complete <- function(output_file, ld_block, output) {
  cmd_failed <- !is.null(attr(output, "status")) && attr(output, "status") != 0
  file_missing <- !file.exists(output_file)

  if (cmd_failed || file_missing) {
    fail_file <- glue::glue('{output_file}_failed.txt')
    writeLines(output, fail_file)
    error_msg <- glue::glue('{ld_block}: FAILED {output_file}')
    flog.error(error_msg)
    stop(error_msg)
  }
}

compile_results <- function(gwas_info) {
  ld_block_dirs <- list.dirs(ld_block_data_dir, recursive = TRUE, full.names = TRUE) |>
    {\(dirs) dirs[!dirs %in% dirname(dirs[-1])]}()

  snp_annotations <- vroom::vroom(
    file.path(variant_annotation_dir, "vep_annotations_hg38.tsv.gz"),
    altrep = FALSE,
    show_col_types = FALSE,
    col_select = c('snp', 'rsid', 'display_snp')
  ) |>
    dplyr::mutate(snp = trimws(snp))
  
  compiled_coloc_pairwise_results_file <- glue::glue('{extracted_study_dir}/compiled_coloc_pairwise_results.tsv')
  compiled_study_extractions_file <- glue::glue('{extracted_study_dir}/compiled_extracted_studies.tsv')
  compiled_coloc_clustered_results_file <- glue::glue('{extracted_study_dir}/compiled_coloc_clustered_results.tsv')
  compiled_associations_file <- glue::glue('{extracted_study_dir}/compiled_associations.tsv')
  lbfs_concatenated_file <- glue::glue('{extracted_study_dir}/gwas_with_lbfs.tsv.gz')

  finemapped_studies_files <- Filter(
    function(file) file.exists(file), 
    glue::glue('{ld_block_dirs}/finemapped_studies.tsv')
  )
  
  finemapped_studies_full <- data.frame()
  if (length(finemapped_studies_files) > 0) {
    study_extractions <- lapply(finemapped_studies_files, function(file) {
      return(vroom::vroom(file, show_col_types = FALSE, col_types = finemapped_column_types))
    }) |>
      dplyr::bind_rows() |>
      dplyr::filter(study == gwas_info$metadata$guid) |>
      dplyr::filter(min_p <= gwas_info$metadata$p_value_threshold) |>
      dplyr::left_join(snp_annotations, by = "snp")
  } else {
    flog.warn(paste(gwas_info$metadata$guid, "No finemapped study files found"))
    study_extractions <- data.frame()
  }

  coloc_clustered_results_files <- Filter(
    function(file) file.exists(file), 
    glue::glue('{ld_block_dirs}/coloc_clustered_results.tsv.gz')
  )
  if (length(coloc_clustered_results_files) > 0) {
    coloc_clustered_results <- lapply(coloc_clustered_results_files, function(file) {
      return(vroom::vroom(file, delim = '\t', show_col_types = FALSE, col_types = coloc_clustered_results_column_types))
    }) |>
      dplyr::bind_rows() |>
      dplyr::group_by(ld_block, component) |>
      dplyr::mutate(coloc_group_id = dplyr::cur_group_id()) |>
      dplyr::ungroup() |>
      dplyr::arrange(coloc_group_id) |>
      dplyr::select(unique_study_id, coloc_group_id, snp, ld_block, h4_connectedness, h3_connectedness)
  } else {
    coloc_clustered_results <- data.frame()
  }

  coloc_pairwise_results_files <- Filter(
    function(file) file.exists(file), 
    glue::glue('{ld_block_dirs}/coloc_pairwise_results.tsv.gz')
  )
  if (length(coloc_pairwise_results_files) > 0) {
    coloc_pairwise_results <- lapply(coloc_pairwise_results_files, function(file) {
      return(vroom::vroom(file, delim = '\t', show_col_types = FALSE, col_types = coloc_pairwise_results_column_types))
    }) |>
      dplyr::bind_rows() |>
      dplyr::filter(study_a == gwas_info$metadata$guid | study_b == gwas_info$metadata$guid) |>
      dplyr::select(unique_study_a, unique_study_b, PP.H3.abf, PP.H4.abf, ld_block, false_positive, false_negative, ignore) |>
      dplyr::rename(unique_study_id_a = unique_study_a, unique_study_id_b = unique_study_b, h3 = PP.H3.abf, h4 = PP.H4.abf)
  } else {
    coloc_pairwise_results <- data.frame()
  }
  
  # associations <- find_associations_for_coloc_clustered_snps(gwas_info, coloc_clustered_results, study_extractions, snp_annotations)
  lbfs_concatenated <- concatenate_file_with_lbfs(gwas_info, study_extractions)

  study_extractions <- study_extractions |>
    dplyr::select(study, unique_study_id, snp, file, chr, bp, min_p, ld_block)

  flog.info(paste(gwas_info$metadata$guid, "Writing compiled results to files"))
  vroom::vroom_write(coloc_pairwise_results, compiled_coloc_pairwise_results_file)
  vroom::vroom_write(study_extractions, compiled_study_extractions_file)
  vroom::vroom_write(coloc_clustered_results, compiled_coloc_clustered_results_file)
  vroom::vroom_write(lbfs_concatenated, lbfs_concatenated_file)
  # vroom::vroom_write(associations, compiled_associations_file)
  
  return(list(
    coloc_pairwise_results = coloc_pairwise_results,
    study_extractions = study_extractions,
    coloc_clustered_results = coloc_clustered_results,
    associations = NULL 
  ))
}

#TODO: go through this, make sure it's working... it's NOT
#TOOD: How do we get the summary stats?  I think we need to pull down the file_with_lbfs from the oracle bucket... ewwwww.
#TODO: Either that or also save BETA and P into the _1 files.
find_associations_for_coloc_clustered_snps <- function(gwas_info, coloc_clustered_results, study_extractions, snp_annotations) {
  if (nrow(coloc_clustered_results) > 0 && length(study_extractions) > 0) {
    finemapped_studies <- study_extractions |>
      dplyr::select(unique_study_id, file_with_lbfs, ld_block) |>
      dplyr::filter(unique_study_id %in% coloc_clustered_results$unique_study_id)
    
    coloc_snps <- coloc_clustered_results |>
      dplyr::filter(unique_study_id %in% finemapped_studies$unique_study_id) |>
      dplyr::select(unique_study_id, snp, coloc_group_id, ld_block)
    
    if (nrow(coloc_snps) > 0 && nrow(finemapped_studies) > 0) {
      # Join to get file paths
      coloc_snps_with_files <- coloc_snps |>
        dplyr::left_join(finemapped_studies, by = c("unique_study_id", "ld_block"))
      
      # Group by file to read each file once
      files_to_read <- coloc_snps_with_files |>
        dplyr::group_by(file_with_lbfs) |>
        dplyr::group_split()
      
      # Read each file and extract relevant SNPs
      associations_list <- lapply(files_to_read, function(file_group) {
        file_path <- file_group$file_with_lbfs[1]
        snp_mapping <- file_group |>
          dplyr::select(snp, unique_study_id, coloc_group_id) |>
          dplyr::distinct()
        target_snps <- unique(snp_mapping$snp)
        
        # Handle relative paths
        if (!grepl('^/', file_path) && !grepl('^study', file_path)) {
          file_path <- file.path(data_dir, file_path)
        } else if (grepl('^study', file_path)) {
          file_path <- file.path(data_dir, file_path)
        }
        
        if (!file.exists(file_path)) {
          flog.warn(paste(gwas_info$metadata$guid, "Association file not found:", file_path))
          return(data.frame())
        }
        
        file_associations <- vroom::vroom(file_path, show_col_types = FALSE) |>
          dplyr::filter(SNP %in% target_snps) |>
          dplyr::select(dplyr::any_of(c('SNP', 'CHR', 'BP', 'EA', 'OA', 'EAF', 'BETA', 'SE', 'P', 'IMPUTED', 'Z'))) |>
          dplyr::rename(snp = SNP)
        
        # Add unique_study_id and coloc_group_id by joining with snp_mapping
        file_associations <- file_associations |>
          dplyr::left_join(snp_mapping, by = "snp") |>
          dplyr::filter(!is.na(unique_study_id))
        
        return(file_associations)
      })
      
      associations <- dplyr::bind_rows(associations_list)
      
      if (nrow(associations) > 0) {
        # Join with snp_annotations to add display_snp and rsid
        associations <- associations |>
          dplyr::left_join(snp_annotations, by = "snp") |>
          dplyr::select(unique_study_id, coloc_group_id, snp, dplyr::any_of(c('rsid', 'display_snp', 'CHR', 'BP', 'EA', 'OA', 'EAF', 'BETA', 'SE', 'P', 'IMPUTED', 'Z')))
        
        flog.info(paste(gwas_info$metadata$guid, "Extracted", nrow(associations), "associations for coloc clustered SNPs"))
        return(associations)
      }
    }
  }
}

concatenate_file_with_lbfs <- function(gwas_info, study_extractions) {
  lbfs_concatenated_file <- glue::glue('{extracted_study_dir}/gwas_with_lbfs.tsv.gz')
  if (nrow(study_extractions) > 0) {
    flog.info(paste(gwas_info$metadata$guid, "Concatenating file_with_lbfs files"))

    lbf_files <- unique(study_extractions$file_with_lbfs)
    if (length(lbf_files) > 0) {
      lbfs_data_list <- lapply(lbf_files, function(file_path) {
        file_data <- vroom::vroom(glue::glue('{data_dir}/{file_path}'), show_col_types = FALSE)
        return(file_data)
      })

      lbfs_data_list <- lbfs_data_list[sapply(lbfs_data_list, nrow) > 0]
      lbfs_concatenated <- dplyr::bind_rows(lbfs_data_list)
    } else {
      lbfs_concatenated <- data.frame()
    }
  } else {
    lbfs_concatenated <- data.frame()
  }
  return(lbfs_concatenated)
}

upload_results <- function(results, gwas_info) {
  bucket_prefix <- glue::glue('gwas_upload/{gwas_info$metadata$guid}/')

  cmd <- paste('oci os object sync',
               '--auth', 'instance_principal',
               '--bucket-name', shQuote(oracle_bucket_name),
               '--src-dir', shQuote(extracted_study_dir),
               '--prefix', shQuote(bucket_prefix))

  status <- system(cmd, wait = TRUE)

  if (status != 0) {
    error_msg <- paste(gwas_info$metadata$guid, 'Failed to upload results to Oracle bucket')
    flog.error(error_msg)
    stop(error_msg)
  }

  flog.info(paste(gwas_info$metadata$guid, 'Successfully uploaded results to Oracle bucket'))
}

send_update_gwas_upload <- function(gwas_info, success, failure_reason, results = NULL) {
  if (!is.na(TEST_RUN)) return()

  api_url <- glue::glue("https://gpmap.opengwas.io/api/v1/gwas/{gwas_info$metadata$guid}")
  flog.info(paste(gwas_info$metadata$guid, "Sending update to API:", api_url))

  if (success) {
    put_body <- list(
      success = success,
      coloc_pairs = results$coloc_pairwise_results,
      study_extractions = results$study_extractions,
      coloc_groups = results$coloc_clustered_results
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
  delete_status <- system(glue::glue('rm -rf {gwas_upload_dir}gwas_upload/{guid}'), wait = TRUE)
  all_delete_status <- c(all_delete_status, delete_status)
  delete_status <- system(glue::glue('rm -rf {gwas_upload_dir}ld_blocks/gwas_upload/{guid}'), wait = TRUE)
  all_delete_status <- c(all_delete_status, delete_status)
  delete_status <- system(glue::glue('rm -rf {gwas_upload_dir}pipeline_metadata/gwas_upload/{guid}'), wait = TRUE)
  all_delete_status <- c(all_delete_status, delete_status)

  if (any(all_delete_status != 0)) {
    error_msg <- paste(guid, "Failed to delete GWAS")
    flog.error(error_msg)
  }
}

main()
