#!/usr/bin/env Rscript

source('../pipeline_steps/constants.R')
library(jsonlite)

parser <- argparser::arg_parser('Sync data files to Oracle server')
parser <- argparser::add_argument(parser, '--all_files', help = 'Sync all files', type = 'logical', flag = TRUE)
parser <- argparser::add_argument(parser, '--db_files', help = 'Sync database files', type = 'logical', flag = TRUE)
parser <- argparser::add_argument(parser, '--svg_files', help = 'Sync SVG files', type = 'logical', flag = TRUE)
parser <- argparser::add_argument(parser, '--upload_server_data', help = 'Sync upload server data', type = 'logical', flag = TRUE)
parser <- argparser::add_argument(parser, '--summary_stats', help = 'Sync summary statistics', type = 'logical', flag = TRUE)

args <- argparser::parse_args(parser)

remote_server_data_dir <- '/oradiskvdb1/data'
remote_server_study_dir <- file.path(remote_server_data_dir, 'study')
remote_server_svg_dir <- '/oradiskvdb1/static/svgs_new'
remote_server_db_dir <- '/oradiskvdb1/db'

main <- function() {
  if (args$all_files || args$db_files) sync_db_files()
  if (args$all_files || args$svg_files) sync_svg_files()
  if (args$all_files || args$upload_server_data) sync_upload_server_data()
  if (args$all_files || args$summary_stats) sync_summary_stats()

  message('Sync completed successfully')
}

run_rsync <- function(src, server = NULL, dest, flags = '-av', extra_args = NULL) {
  cmd_parts <- c('rsync', flags)
  
  if (!is.null(extra_args)) {
    cmd_parts <- c(cmd_parts, extra_args)
  }
  
  if (!is.null(server)) {
    cmd_parts <- c(cmd_parts, shQuote(src), paste0('opc@',server, ':', shQuote(dest)))
  } else {
    cmd_parts <- c(cmd_parts, shQuote(src), shQuote(dest))
  }
  
  cmd <- paste(cmd_parts, collapse = ' ')
  message(paste('Running:', cmd))
  status <- system(cmd)
  if (status != 0) {
    stop(paste('rsync command failed with status', status))
  }
}

sync_db_files <- function() {
  message('copying db files to oracle server')
  
  db_files <- list(
    'studies.db' = 'studies_new.db',
    'coloc_pairs.db' = 'coloc_pairs_new.db',
    'associations.db' = 'associations_new.db',
    'ld.db' = 'ld_new.db'
  )
  
  for (src_file in names(db_files)) {
    src_path <- file.path(latest_results_dir, src_file)
    dest_file <- db_files[[src_file]]
    dest_path <- file.path(remote_server_db_dir, dest_file)
    
    if (!file.exists(src_path)) {
      warning(paste('Source file does not exist:', src_path))
      next
    }
    
    run_rsync(src_path, oracle_api_server, dest_path)
  }
}

sync_svg_files <- function() {
  message('copying svg files to oracle server')
  
  traits_src <- file.path(svg_dir, 'traits', '')
  groups_src <- file.path(svg_dir, 'groups', '')
  
  if (dir.exists(file.path(svg_dir, 'traits'))) {
    run_rsync(traits_src, oracle_api_server, file.path(remote_server_svg_dir, 'traits'))
  } else {
    warning(paste('Source directory does not exist:', traits_src))
  }
  
  if (dir.exists(file.path(svg_dir, 'groups'))) {
    run_rsync(groups_src, oracle_api_server, file.path(remote_server_svg_dir, 'groups'))
  } else {
    warning(paste('Source directory does not exist:', groups_src))
  }
}

sync_upload_server_data <- function() {
  message('VEP annotation files to upload server')
  vep_annotation_file <- file.path(variant_annotation_dir, "vep_annotations_hg38.tsv.gz")
  run_rsync(
    src = vep_annotation_file,
    server = oracle_upload_server,
    dest = file.path(remote_server_data_dir, 'variant_annotation', basename(vep_annotation_file)),
    flags = '-av --mkpath',
  )

  message('preparing ld_blocks directory')
  ld_blocks_rsync_dir <- file.path(server_sync_dir, 'ld_blocks')
  dir.create(ld_blocks_rsync_dir, recursive = TRUE, showWarnings = FALSE)

  old_wd <- getwd()
  setwd(ld_block_data_dir)
  run_rsync(
    src = './',
    dest = ld_blocks_rsync_dir,
    flags = '-aRv --prune-empty-dirs',
    extra_args = c('--include=*/', '--include=*coloc_pairwise_results.tsv.gz', '--exclude=*')
  )
  setwd(old_wd)
  
  finemapped_src_files <- list.files(
    ld_block_data_dir, 
    pattern = 'finemapped_studies.tsv$', 
    full.names = TRUE, 
    recursive = TRUE
  )

  # for (file in finemapped_src_files) {
  #   finemapped_studies <- vroom::vroom(file, show_col_types = F, delim = '\t')
  #   finemapped_studies$file <- gsub(paste0('^', data_dir), '', finemapped_studies$file)
  #   finemapped_studies$svg_file <- gsub(paste0('^', data_dir), '', finemapped_studies$svg_file)
  #   finemapped_studies$file_with_lbfs <- gsub(paste0('^', data_dir), '', finemapped_studies$file_with_lbfs)

  #   rsync_file <- gsub(ld_block_data_dir, ld_blocks_rsync_dir, file)
  #   print(rsync_file)
  #   vroom::vroom_write(finemapped_studies, rsync_file)
  # }
  
  message('rsyncing ld_blocks to oracle server')
  run_rsync(
    paste0(ld_blocks_rsync_dir, '/'),
    oracle_upload_server,
    file.path(remote_server_data_dir, 'ld_blocks_new'),
    flags = '-av'
  )

  message('rsyncing ld reference panel to oracle server')
  run_rsync(
    paste0(ld_reference_panel_dir, '/'),
    oracle_upload_server,
    file.path(remote_server_data_dir, 'ld_reference_panel_hg38'),
    flags = '-av'
  )

  message('rsyncing finemapped files to oracle server')
  old_wd <- getwd()
  setwd(study_dir)

  run_rsync(
    src = current_dir,
    server = oracle_upload_server,
    dest = remote_server_study_dir,
    flags = '-av --prune-empty-dirs',
    extra_args = c('--include=*/', '--include=*finemapped/*_{1,2,3,4,5,6,7,8,9,10}.tsv.gz', '--exclude=*')
  )
  setwd(old_wd)
}

# Helper function to check if OCI CLI is configured
check_oci_config <- function() {
  # Check if OCI CLI is available
  oci_check <- system('which oci', ignore.stdout = TRUE, ignore.stderr = TRUE)
  if (oci_check != 0) {
    stop('OCI CLI not found. Please install oci-cli.')
  }
  
  # Check if OCI config exists or if we need to use environment variables
  config_file <- Sys.getenv('OCI_CLI_CONFIG_FILE', '~/.oci/config')
  if (!file.exists(path.expand(config_file))) {
    # Try to use environment variables for authentication
    if (Sys.getenv('OCI_CLI_USER') == '' || Sys.getenv('OCI_CLI_TENANCY') == '') {
      stop('OCI CLI not configured. Please set OCI_CLI_CONFIG_FILE or OCI_CLI_USER/OCI_CLI_TENANCY environment variables.')
    }
  }
}

# Helper function to get object metadata from Oracle bucket
get_remote_object_info <- function(bucket_name, namespace, object_name) {
  cmd <- paste('oci os object head',
               '--bucket-name', shQuote(bucket_name),
               '--namespace-name', shQuote(namespace),
               '--name', shQuote(object_name),
               '--output json 2>/dev/null')
  
  result <- tryCatch({
    output <- system(cmd, intern = TRUE, ignore.stderr = TRUE)
    if (length(output) > 0) {
      jsonlite::fromJSON(paste(output, collapse = ''))
    } else {
      NULL
    }
  }, error = function(e) NULL)
  
  return(result)
}

# Helper function to upload file to Oracle bucket (only if needed)
upload_to_oracle_bucket <- function(local_file, bucket_name, namespace, object_name) {
  # Get local file info
  if (!file.exists(local_file)) {
    warning(paste('Local file does not exist:', local_file))
    return(FALSE)
  }
  
  local_size <- file.info(local_file)$size
  local_mtime <- file.info(local_file)$mtime
  
  # Check if object exists in bucket
  remote_info <- get_remote_object_info(bucket_name, namespace, object_name)
  
  # Determine if upload is needed
  needs_upload <- TRUE
  if (!is.null(remote_info)) {
    remote_size <- as.numeric(remote_info$contentLength)
    remote_mtime <- as.POSIXct(remote_info$lastModified, format = "%Y-%m-%dT%H:%M:%S")
    
    # Only upload if size or modification time differs
    if (remote_size == local_size && abs(as.numeric(difftime(remote_mtime, local_mtime, units = "secs"))) < 2) {
      needs_upload <- FALSE
      message(paste('Skipping', basename(local_file), '- already up to date'))
    }
  }
  
  if (needs_upload) {
    message(paste('Uploading', basename(local_file), 'to Oracle bucket'))
    cmd <- paste('oci os object put',
                 '--bucket-name', shQuote(bucket_name),
                 '--namespace-name', shQuote(namespace),
                 '--name', shQuote(object_name),
                 '--file', shQuote(local_file),
                 '--force')
    
    status <- system(cmd)
    if (status != 0) {
      stop(paste('Failed to upload', local_file, 'to Oracle bucket'))
    }
    return(TRUE)
  }
  
  return(FALSE)
}

sync_summary_stats <- function() {
  message('syncing finemapped GWAS results to Oracle bucket')
  
  # Get Oracle bucket configuration from environment variables
  oracle_bucket <- Sys.getenv('ORACLE_BUCKET_NAME')
  oracle_namespace <- Sys.getenv('ORACLE_NAMESPACE')
  
  if (oracle_bucket == '' || oracle_namespace == '') {
    stop('ORACLE_BUCKET_NAME and ORACLE_NAMESPACE environment variables must be set')
  }
  
  check_oci_config()
  
  # Find all files matching the pattern
  old_wd <- getwd()
  setwd(extracted_study_dir)
  on.exit(setwd(old_wd))
  
  # Find all finemapped files with _with_lbfs.tsv.gz pattern
  all_files <- list.files(
    path = '.',
    pattern = '_with_lbfs\\.tsv\\.gz$',
    full.names = TRUE,
    recursive = TRUE
  )
  
  # Filter to only include files in finemapped subdirectories
  all_files <- all_files[grepl('/finemapped/', all_files)]
  
  if (length(all_files) == 0) {
    message('No files found matching pattern')
    return()
  }
  
  message(paste('Found', length(all_files), 'files to check'))
  
  uploaded_count <- 0
  skipped_count <- 0
  
  for (local_file in all_files) {
    # Create object name in bucket (preserve directory structure)
    object_name <- gsub('^\\./', '', local_file)
    
    if (upload_to_oracle_bucket(local_file, oracle_bucket, oracle_namespace, object_name)) {
      uploaded_count <- uploaded_count + 1
    } else {
      skipped_count <- skipped_count + 1
    }
  }
  
  message(paste('Upload complete:', uploaded_count, 'uploaded,', skipped_count, 'skipped'))
}

main()