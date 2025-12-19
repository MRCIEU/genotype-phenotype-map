#!/usr/bin/env Rscript

source('../pipeline_steps/constants.R')

parser <- argparser::arg_parser('Sync data files to Oracle server')
parser <- argparser::add_argument(parser, '--all_files', help = 'Sync all files', type = 'logical', flag = TRUE)
parser <- argparser::add_argument(parser, '--db_files', help = 'Sync database files', type = 'logical', flag = TRUE)
parser <- argparser::add_argument(parser, '--svg_files', help = 'Sync SVG files', type = 'logical', flag = TRUE)
parser <- argparser::add_argument(parser, '--ld_blocks', help = 'Sync LD blocks', type = 'logical', flag = TRUE)
parser <- argparser::add_argument(parser, '--summary_stats', help = 'Sync summary statistics', type = 'logical', flag = TRUE)

args <- argparser::parse_args(parser)

remote_server_data_dir <- '/oradiskvdb1/data'
remote_server_study_dir <- file.path(remote_server_data_dir, 'study')
remote_server_svg_dir <- '/oradiskvdb1/static/svgs_new'
remote_server_db_dir <- '/oradiskvdb1/db'

main <- function() {
  if (args$all_files || args$db_files) sync_db_files()
  if (args$all_files || args$svg_files) sync_svg_files()
  if (args$all_files || args$ld_blocks) sync_ld_blocks()
  if (args$all_files || args$summary_stats) sync_summary_stats()

  message('Sync completed successfully')
}

run_rsync <- function(src, server, dest, flags = '-av', extra_args = NULL) {
  cmd_parts <- c('rsync', flags)
  
  if (!is.null(extra_args)) {
    cmd_parts <- c(cmd_parts, extra_args)
  }
  
  cmd_parts <- c(cmd_parts, shQuote(src), paste0(server, ':', shQuote(dest)))
  
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

sync_ld_blocks <- function() {
  message('preparing ld_blocks directory')
  
  ld_blocks_rsync_dir <- file.path(server_sync_dir, 'ld_blocks')
  dir.create(ld_blocks_rsync_dir, recursive = TRUE, showWarnings = FALSE)
  
  ld_blocks_src <- file.path(data_dir, 'ld_blocks')
  
  # Find and copy finemapped_studies.tsv files
  finemapped_files <- list.files(
    ld_blocks_src, 
    pattern = 'finemapped_studies.tsv', 
    full.names = TRUE, 
    recursive = TRUE
  )
  
  for (file in finemapped_files) {
    rel_path <- gsub(paste0('^', data_dir), '', file)
    dest_path <- file.path(server_sync_dir, rel_path)
    dest_dir <- dirname(dest_path)
    dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
    file.copy(file, dest_path, overwrite = TRUE)
  }
  
  # Find and copy coloc_pairwise_results.tsv.gz files
  coloc_files <- list.files(
    ld_blocks_src, 
    pattern = 'coloc_pairwise_results.tsv.gz', 
    full.names = TRUE, 
    recursive = TRUE
  )
  
  for (file in coloc_files) {
    rel_path <- gsub(paste0('^', data_dir), '', file)
    dest_path <- file.path(server_sync_dir, rel_path)
    dest_dir <- dirname(dest_path)
    dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
    file.copy(file, dest_path, overwrite = TRUE)
  }
  
  # Replace paths in finemapped_studies.tsv files
  finemapped_dest_files <- list.files(
    ld_blocks_rsync_dir, 
    pattern = 'finemapped_studies.tsv', 
    full.names = TRUE, 
    recursive = TRUE
  )
  
  for (file in finemapped_dest_files) {
    content <- readLines(file)
    content <- gsub(study_dir, remote_server_study_dir, content, fixed = TRUE)
    writeLines(content, file)
  }
  
  # Replace paths in coloc_pairwise_results.tsv.gz files
  coloc_dest_files <- list.files(
    ld_blocks_rsync_dir, 
    pattern = 'coloc_pairwise_results.tsv.gz', 
    full.names = TRUE, 
    recursive = TRUE
  )
  
  for (file in coloc_dest_files) {
    # Read compressed file, replace, write back
    con_in <- gzfile(file, 'rt')
    content <- readLines(con_in)
    close(con_in)
    
    content <- gsub(study_dir, remote_server_study_dir, content, fixed = TRUE)
    
    con_out <- gzfile(paste0(file, '.tmp'), 'wt')
    writeLines(content, con_out)
    close(con_out)
    
    file.rename(paste0(file, '.tmp'), file)
  }
  
  message('rsyncing ld_blocks to oracle server')
  run_rsync(
    paste0(ld_blocks_rsync_dir, '/'),
    oracle_pipeline_server,
    file.path(remote_server_data_dir, 'ld_blocks_new'),
    flags = '-aRv'
  )
}

sync_summary_stats <- function() {
  message('syncing finemapped GWAS results to oracle server')
  
  old_wd <- getwd()
  setwd(study_dir)
  on.exit(setwd(old_wd))
  
  current_dir <- paste0(getwd(), '/')
  
  run_rsync(
    src = current_dir,
    server = oracle_pipeline_server,
    dest = remote_server_study_dir,
    flags = '-aRv --prune-empty-dirs',
    extra_args = c('--include=*/', '--include=*finemapped/*_with_lbfs.tsv.gz', '--exclude=*')
  )
}

main()