#!/usr/bin/env Rscript

source('../pipeline_steps/constants.R')
library(jsonlite)

parser <- argparser::arg_parser('Sync data files to Oracle server')
parser <- argparser::add_argument(parser, '--all_files', help = 'Sync all files', type = 'logical', flag = TRUE)
parser <- argparser::add_argument(parser, '--db_files', help = 'Sync database files', type = 'logical', flag = TRUE)
parser <- argparser::add_argument(parser, '--svg_files', help = 'Sync SVG files', type = 'logical', flag = TRUE)
parser <- argparser::add_argument(parser, '--upload_server_data', help = 'Sync upload server data', type = 'logical', flag = TRUE)
parser <- argparser::add_argument(parser, '--summary_stats', help = 'Sync summary statistics', type = 'logical', flag = TRUE)
parser <- argparser::add_argument(parser, '--rdsf', help = 'Sync to RDSF', type = 'logical', flag = TRUE)

args <- argparser::parse_args(parser)

rdsf_dir <- '/mnt/GPMap_rdsf/Data-Bris/working/'
upload_server_data_dir <- '/oradiskvdb1/data'
api_server_svg_dir <- '/home/opc/gpmap_data/static/svgs'
api_server_db_dir <- '/home/opc/gpmap_data/db'

main <- function() {
  if (args$all_files || args$db_files) sync_db_files()
  if (args$all_files || args$svg_files) sync_svg_files()
  if (args$all_files || args$upload_server_data) sync_upload_server_data()
  if (args$all_files || args$summary_stats) sync_summary_stats()
  if (args$all_files || args$rdsf) sync_to_rdsf()

  message('Sync completed successfully')
}

run_rsync <- function(src, server = NULL, dest, flags = '-avu', extra_args = NULL) {
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
    dest_path <- file.path(api_server_db_dir, dest_file)
    
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
    run_rsync(traits_src, oracle_api_server, file.path(api_server_svg_dir, 'traits'))
  } else {
    warning(paste('Source directory does not exist:', traits_src))
  }
  
  if (dir.exists(file.path(svg_dir, 'groups'))) {
    run_rsync(groups_src, oracle_api_server, file.path(api_server_svg_dir, 'groups'))
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
    dest = file.path(upload_server_data_dir, 'variant_annotation', basename(vep_annotation_file)),
    flags = '-avu --mkpath',
  )

  message('liftover data to oracle server')
  run_rsync(
    src = liftover_dir,
    server = oracle_upload_server,
    dest = file.path(upload_server_data_dir, 'liftover')
  )

  message('preparing ld_blocks directory')
  ld_blocks_rsync_dir <- file.path(server_sync_dir, 'ld_blocks')
  dir.create(ld_blocks_rsync_dir, recursive = TRUE, showWarnings = FALSE)

  run_rsync(
    src = ld_block_data_dir,
    dest = ld_blocks_rsync_dir,
    flags = '-avu --prune-empty-dirs',
    extra_args = c('--include=*/', '--include=*coloc_pairwise_results.tsv.gz', '--exclude=*')
  )
  
  finemapped_src_files <- list.files(
    ld_block_data_dir, 
    pattern = 'finemapped_studies.tsv$', 
    full.names = TRUE, 
    recursive = TRUE
  )

  for (file in finemapped_src_files) {
    finemapped_studies <- vroom::vroom(file, show_col_types = F, delim = '\t')
    finemapped_studies$file <- gsub(paste0('^', data_dir), '', finemapped_studies$file)
    finemapped_studies$svg_file <- gsub(paste0('^', data_dir), '', finemapped_studies$svg_file)
    finemapped_studies$file_with_lbfs <- gsub(paste0('^', data_dir), '', finemapped_studies$file_with_lbfs)

    rsync_file <- gsub(ld_block_data_dir, ld_blocks_rsync_dir, file)
    print(rsync_file)
    vroom::vroom_write(finemapped_studies, rsync_file)
  }
  
  message('rsyncing ld_blocks to oracle server')
  old_wd <- getwd()
  setwd(ld_blocks_rsync_dir)

  run_rsync(
    '.',
    oracle_upload_server,
    file.path(upload_server_data_dir, 'ld_blocks')
  )
  setwd(old_wd)

  message('rsyncing ld reference panel to oracle server')
  run_rsync(
    paste0(ld_reference_panel_dir, '/'),
    oracle_upload_server,
    file.path(upload_server_data_dir, 'ld_reference_panel_hg38')
  )

  message('rsyncing finemapped files to oracle server')
  # rsync doesn't support brace expansion, so we need separate include patterns for each number
  extra_study_dir_args <- c('--include=*/')
  for (i in 1:10) {
    extra_study_dir_args <- c(extra_study_dir_args, paste0('--include=*finemapped/*_', i, '.tsv.gz'))
  }
  extra_study_dir_args <- c(extra_study_dir_args, '--exclude=*')
  
  run_rsync(
    src = extracted_study_dir,
    server = oracle_upload_server,
    dest = file.path(upload_server_data_dir, 'study'),
    flags = '-avu --prune-empty-dirs',
    extra_args = extra_study_dir_args
  )
}

sync_summary_stats <- function() {
  message('syncing finemapped GWAS results to Oracle bucket')

  cmd <- paste('oci os object sync',
               '--bucket-name', shQuote(oracle_bucket_name),
               '--src-dir', shQuote(extracted_study_dir),
               '--prefix', shQuote('study/'),
               '--include', shQuote('*/finemapped/*_with_lbf.tsv.gz'))
  
  message(paste('Running:', cmd))
  status <- system(cmd, wait = TRUE)
}

sync_to_rdsf <- function() {
  original_wd <- getwd()
  on.exit(setwd(original_wd))
  
  rdsf_data_dir <- file.path(rdsf_dir, 'data')
  dir.create(rdsf_data_dir, recursive = TRUE, showWarnings = FALSE)
  
  create_and_store_tar(liftover_dir, 'liftover.tar.zst', rdsf_data_dir)
  create_and_store_tar(svg_dir, 'svgs.tar.zst', rdsf_data_dir)
  create_and_store_tar(ld_block_data_dir, 'ld_blocks.tar.zst', rdsf_data_dir)
  create_and_store_tar(ld_reference_panel_dir, 'ld_reference_panel_hg38.tar.zst', rdsf_data_dir)
  create_and_store_tar(variant_annotation_dir, 'variant_annotation.tar.zst', rdsf_data_dir)
  create_and_store_tar(extracted_study_dir, 'study.tar.zst', rdsf_data_dir)
  # create_and_store_tar(study_sequencing_dir, 'study_sequencing.tar.zst', rdsf_data_dir)
  
  message('RDSF sync completed successfully')
}

# TO inflate a tar.zst file, tar -xf your_file.tar.zst (or tar --zstd -xf your_file.tar.zst)
create_and_store_tar <- function(source_dir, file_name, dest_dir, zip = TRUE) {
  if (!dir.exists(source_dir)) {
    warning(paste('Source directory does not exist:', source_dir))
    return(invisible(NULL))
  }
  
  tar_path <- file.path(data_dir, file_name)
  message(paste('zipping', tar_path))
  if (zip) {
    tar_cmd <- paste('tar', '-c --use-compress-program=zstd -f ', shQuote(tar_path), shQuote(source_dir))
  }
  else {
    tar_cmd <- paste('tar', '-c -f ', shQuote(tar_path), shQuote(source_dir))
  }
  status <- system(tar_cmd)
  if (status != 0) {
    stop(paste('Failed to create tar file:', basename(tar_path)))
  }
  fs::file_move(tar_path, dest_dir)
  
  message(paste('completed', basename(tar_path)))
}

main()