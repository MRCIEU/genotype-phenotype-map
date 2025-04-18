source('constants.R')

parser <- argparser::arg_parser('Sync finemapped data to server')
parser <- argparser::add_argument(parser, '--studies_to_process', help = 'Studies to process', type = 'character')
parser <- argparser::add_argument(parser, '--study_extractions', help = 'Study extractions', type = 'character')

args <- argparser::parse_args(parser)

main <- function() {
  dir.create(server_sync_dir, showWarnings = FALSE, recursive = TRUE)

  finemapped_files <- list.files(ld_block_data_dir, recursive = TRUE, full.names = TRUE) |>
    grep(pattern = 'finemapped_studies.tsv', value = TRUE)

  # renamed_files <- lapply(finemapped_files, function(file) {
  #   print(file)
  #   finemapped <- vroom::vroom(file, show_col_types = FALSE)
  #   original_files <- finemapped$file

  #   finemapped <- finemapped |>
  #     dplyr::mutate(file = gsub(data_dir, oracle_data_dir, file)) |>
  #     dplyr::mutate(file = gsub('/finemapped/', '_finemapped_', file))
  #   new_files <- finemapped$file
    
  #   new_location <- sub('ld_blocks', 'rsync_to_server/ld_blocks', file)
  #   dir.create(dirname(new_location), showWarnings = FALSE, recursive = TRUE)
  #   # vroom::vroom_write(finemapped, new_location)
  #   return (data.frame(original_file = original_files, new_file = new_files))
  # }) |> dplyr::bind_rows()

  finemap_file_map <- file.path(server_sync_dir, 'finemap_file_map.csv')

  study_to_process <- vroom::vroom(args$studies_to_process, show_col_types = FALSE)
  study_extractions <- vroom::vroom(args$study_extractions, show_col_types = FALSE) |>
    dplyr::filter(study %in% study_to_process$study_name)
  
  already_synced <- vroom::vroom(file.path(server_sync_dir, 'all_study_files.txt'), show_col_types = FALSE, col_names = F, delim = ',')$X1
  
  renamed_files <- vroom::vroom(file.path(server_sync_dir, 'all_finemap_file_map.csv'), show_col_types = FALSE, col_names = FALSE)
  print(nrow(renamed_files))
  print(head(renamed_files))
  renamed_files <- renamed_files |>
    dplyr::filter(X2 %in% already_synced)

  print(nrow(renamed_files))
  print(head(renamed_files))
  vroom::vroom_write(renamed_files, finemap_file_map, col_names = FALSE, delim = ',')

  send_data_command <- glue::glue('while IFS=, read -r src dst; do
    tar -C "$(dirname "$src")" --transform="s|.*|$dst|" -cf - "$(basename "$src")"
  done < {finemap_file_map} | ssh {oracle_server} "tar -xf - -C {oracle_data_dir}/study/"
  ')

  print(send_data_command)
  system(send_data_command, wait = TRUE)
}

main()