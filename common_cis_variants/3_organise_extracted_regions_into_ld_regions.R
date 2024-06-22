source('constants.R')

ld_regions <- vroom::vroom("data/ld_regions.tsv", show_col_types = F)
current_state <- vroom::vroom(paste0(data_dir, "/pipeline_metadata/current_state.tsv"), show_col_types = F) |>
  dplyr::filter(extracted == F)

all_updated_ld_blocks <- apply(current_state, 1, function(study) {
  study_name <- study[["study_name"]]
  p_value_threshold <- study[["p_value_threshold"]]
  study_dir <- study[["extracted_location"]]
  category <- study[["category"]]
  sample_size <- study[["sample_size"]]
  extracted_snps <- vroom::vroom(paste0(study_dir, "/extracted_snps.tsv"), show_col_types = F)

  if(nrow(extracted_snps) == 0) {
    return(data.frame())
  }
  updated_ld_blocks <- apply(extracted_snps, 1, function(extracted) {
    bp <- as.numeric(extracted[['BP']])
    extracted_chr <- as.numeric(extracted[['CHR']])
    ancestry <- extracted[['ANCESTRY']]

    ld_block <- dplyr::filter(ld_regions, chr == extracted_chr & start < bp & stop > bp & pop == ancestry)
    if (nrow(ld_block) > 1) stop(paste("Error: More than 1 LD Block associated with", extracted_chr, bp))

    ld_block_data <- paste0(ld_block_data_dir, ancestry, "/", extracted_chr, "/", ld_block$start, "_", ld_block$stop)
    ld_block_results <- paste0(ld_block_results_dir, ancestry, "/", extracted_chr, "/", ld_block$start, "_", ld_block$stop)

    study_file <- paste0(study_dir, 'original/', ancestry, "_", extracted_chr, "_", bp, ".tsv")
    ld_block$data_dir <- ld_block_data
    ld_block$results_dir <- ld_block_results
    if(!dir.exists(ld_block_data)) dir.create(ld_block_data, recursive = T)
    if(!dir.exists(ld_block_results)) dir.create(ld_block_results, recursive = T)

    extracted_studies_file <- paste0(ld_block_data, "/extracted_studies.tsv")
    extracted_studies <- tibble::tribble(~study, ~file, ~chr, ~bp, ~p_value_threshold, ~category, ~sample_size,
                                         study_name, study_file, extracted_chr, bp, p_value_threshold, category, sample_size
    )
    if (file.exists(extracted_studies_file)) {
      existing_extracted_studies <- vroom::vroom(extracted_studies_file, show_col_types = F)
      extracted_studies <- rbind(existing_extracted_studies, extracted_studies)
      extracted_studies <- extracted_studies[!duplicated(extracted_studies), ]
    }
    vroom::vroom_write(extracted_studies, extracted_studies_file)

    study_in_ld_block <- paste0(ld_block_data, "original/", study_name, "_", extracted_chr, "_", bp, ".tsv")
    file.symlink(study_file, study_in_ld_block)

    return(ld_block)
  }) |> dplyr::bind_rows()
  return(updated_ld_blocks)
})

all_updated_ld_blocks <- dplyr::bind_rows(all_updated_ld_blocks) |> dplyr::distinct() |> dplyr::arrange(chr)
vroom::vroom_write(all_updated_ld_blocks, paste0(pipeline_metadata_dir, "updated_ld_blocks_to_colocalise.tsv"))
