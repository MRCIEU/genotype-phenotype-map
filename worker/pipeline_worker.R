setwd('pipeline_steps')
source('constants.R')

parser <- argparser::arg_parser('Pipeline worker')
parser <- argparser::add_argument(parser, '--reprocess_dlq', help = 'Reprocess DLQ messages', type = 'logical', default = F)
args <- argparser::parse_args(parser)

redis_conn <- redux::hiredis(
  host = Sys.getenv("REDIS_HOST", "redis"),
  port = as.numeric(Sys.getenv("REDIS_PORT", 6379))
)

process_gwas <- 'process_gwas'
process_gwas_dlq <- glue::glue('{process_gwas}_dlq')

main <- function() {
  if (args$reprocess_dlq) {
    retry_dlq_messages()
    return()
  }

  while(TRUE) {
    if (is.na(TEST_RUN)) {
      redis_message <- redis_conn$BRPOP(process_gwas, timeout = 0)
    } else {
      redis_message <- '../tests/test_data/redis_message.json'
    }
    
    if (!is.null(redis_message)) {
      gwas_info <- jsonlite::fromJSON(redis_message[[2]])
      message("Received new message from queue with guid: ", gwas_info$metadata$guid)

      processed <- process_message(gwas_info)

      if (!is.null(processed)) {
        message("Failed to process message with guid: ", gwas_info$metadata$guid)

        message_and_error <- list(
          original_message = gwas_info,
          error = processed,
          timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
        )
        redis_conn$LPUSH(process_gwas_dlq, jsonlite::toJSON(message_and_error, auto_unbox = TRUE))
      } else {
        message("Successfully processed message with guid: ", gwas_info$metadata$guid)
      }
    }
    Sys.sleep(5)
    if (!is.na(TEST_RUN)) {
      break
    }
  }
}

# Function to process a single message
process_message <- function(gwas_info) {
  tryCatch({
    update_directories_for_worker(gwas_info$metadata$guid)

    dir.create(pipeline_metadata_dir, recursive = T, showWarnings = F)
    dir.create(extracted_study_dir, recursive = T, showWarnings = F)
    dir.create(ld_block_data_dir, recursive = T, showWarnings = F)

    message('extracting regions')
    
    gwas_info$metadata$file_location <- gwas_info$file_location
    if (grepl('\\.vcf', gwas_info$file_location)) {
      gwas_info$metadata$file_type <- 'vcf'
    } else {
      gwas_info$metadata$file_type <- 'csv'
    }

    create_study_metadata_files(gwas_info)
    extract_regions <- glue::glue("Rscript extract_regions_from_summary_stats.R",
      " --worker_guid {gwas_info$metadata$guid}")
    system(extract_regions, wait = T)

    check_step_complete(glue::glue('{extracted_study_dir}/extracted_snps.tsv'), 'extracted_snps.tsv')
    extracted_regions <- vroom::vroom(glue::glue('{extracted_study_dir}/extracted_snps.tsv'), show_col_types = F)

    message('organising ld blocks')
    ld_blocks_to_colocalise_file <- glue::glue('{pipeline_metadata_dir}/updated_ld_blocks_to_colocalise.tsv')
    organise_ld_blocks <- glue::glue("Rscript organise_extracted_regions_into_ld_blocks.R",
      " --output_file {ld_blocks_to_colocalise_file}",
      " --worker_guid {gwas_info$metadata$guid}")
    system(organise_ld_blocks, wait = T, intern = T)
    check_step_complete(ld_blocks_to_colocalise_file, 'updated_ld_blocks_to_colocalise.tsv')
    ld_blocks_to_colocalise <- vroom::vroom(ld_blocks_to_colocalise_file, show_col_types = F)


    parallel_block_processing <- 4
    blocks <- head(ld_blocks_to_colocalise$ld_block, 2)
    # parallel::mclapply(ld_blocks_to_colocalise$ld_block, mc.cores=parallel_block_processing, function(block) {
    lapply(blocks, function(block) {
      ld_info <- ld_block_dirs(block)

      output_files <- list(
        standardised = glue::glue('{ld_info$ld_block_data}/standardisation_complete'),
        imputed = glue::glue('{ld_info$ld_block_data}/imputation_complete'),
        finemapped = glue::glue('{ld_info$ld_block_data}/finemapping_complete'),
        coloc = glue::glue('{ld_info$ld_block_data}/colocalisation_complete')
      )

      message('standardising regions')
      standardise_regions <- glue::glue("Rscript standardise_studies_in_ld_block.R",
        " --ld_block {block} ", 
        " --completed_output_file {output_files$standardised}",
        " --worker_guid {gwas_info$metadata$guid}")
      system(standardise_regions, wait = T, intern = T)
      check_step_complete(output_files$standardised, block)

      message('imputing regions')
      impute_regions <- glue::glue("Rscript impute_studies_in_ld_block.R",
        " --ld_block {block} ",
        " --completed_output_file {output_files$imputed}",
        " --worker_guid {gwas_info$metadata$guid}")
      system(impute_regions, wait = T, intern = T)
      check_step_complete(output_files$imputed, block)

      message('finemapping regions')
      finemap_regions <- glue::glue("Rscript finemap_studies_in_ld_block.R",
        " --ld_block {block} ",
        " --completed_output_file {output_files$finemapped}",
        " --worker_guid {gwas_info$metadata$guid}")
      system(finemap_regions, wait = T, intern = T)
      check_step_complete(output_files$finemapped, block)

      message('colocalising regions')
      coloc_regions <- glue::glue("Rscript colocalise_studies_in_ld_block.R",
        " --ld_block {block} ",
        " --completed_output_file {output_files$coloc}",
        " --worker_guid {gwas_info$metadata$guid}")
      system(coloc_regions, wait = T, intern = T)
      check_step_complete(output_files$coloc, block)
    })

    message('sending email')
    send_email(gwas_info)
    
    return()
    
  }, error = function(e) {
    message("Error processing message: ", e$message)
    
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

    return(e$message)
  })
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

send_email <- function(message) {
  tryCatch({
    email_body <- list(
      sendmailR::mime_part_html(
        glue::glue(
          "<html>",
          "<body>",
          "<p>Thanks for using the Genotype-Phenotype Map!</p>",
          "<p>Your results are ready to view. <a href='{gpm_website_data$url}/results/{message$guid}'>Click here</a> to view your results.</p>",
          "<p>If you have any questions, please <a href='{gpm_website_data$contact}'>contact us here</a>.</p>",
          "<p>Best regards,<br>The Genotype-Phenotype Map Team</p>",
          "</body>",
          "</html>"
        ),
        type = "text/html; charset=UTF-8"
      )
    )

    sendmailR::sendmail(
      from = Sys.getenv("SMTP_FROM", "noreply@bristol.ac.uk"),
      to = message$email,
      subject = glue::glue("Genotype-Phenotype Map Results for {message$name}"),
      msg = email_body,
      control = list(smtpServer = "localhost")
    )
  
    return(TRUE)
  }, error = function(e) {
    message(sprintf("Error sending email: %s", e$message))
    return(FALSE)
  })
}

retry_dlq_messages <- function(max_retries = 10) {
  log_info("Starting DLQ retry process")
  retried_count <- 0
  
  dlq_message <- redis_conn$BRPOP(process_gwas_dlq, timeout = 1)
  
  if (is.null(dlq_message)) {
    log_info("No more messages in DLQ")
    break
  }
  
  tryCatch({
    error_details <- jsonlite::fromJSON(dlq_message[[2]])
    redis_conn$LPUSH(process_gwas, error_details$original_message)
    
    log_info(sprintf("Retried message from DLQ (original error: %s, timestamp: %s)", 
                    error_details$error, 
                    error_details$timestamp))
    
    retried_count <- retried_count + 1
    
  }, error = function(e) {
    message("Error processing DLQ message: ", e$message)
    redis_conn$LPUSH(process_gwas_dlq, dlq_message[[2]])
  })
}

check_step_complete <- function(output_file, ld_block) {
  if (!file.exists(output_file)) {
    rlang::abort(
      message = glue::glue('Step {output_file} failed'),
      class = "pipeline_worker_error",
      data = list(
        output_file = output_file,
        ld_block = ld_block
      )
    )
  }
}

main()
