source('constants.R')
library(argparser, quietly = TRUE)

ld_regions <- vroom::vroom("data/ld_regions.tsv")
current_state <- vroom::vroom(paste0(data_dir, "/pipeline_metadata/current_state.tsv"), show_col_types = F) |>
  dplyr::filter(extracted == F)

updated_ld_blocks <- apply(current_state, 1, function(study) {
  study_name <- study[["study_name"]]
  p_value_threshold <- study[["p_value_threshold"]]
  study_dir <- study[['extracted_location']]
  extracted_snps <- vroom::vroom(paste0(study_dir, "/extracted_snps.tsv"), show_col_types = F)

  updated_ld_blocks <- apply(extracted_snps, 1, function(extracted) {
    bp <- as.numeric(extracted[['BP']])
    extracted_chr <- as.numeric(extracted[['CHR']])
    ancestry <- extracted[['ANCESTRY']]

    ld_block <- dplyr::filter(ld_regions, chr == extracted_chr & start < bp & stop > bp & pop == ancestry)
    if (nrow(ld_block) > 1) stop(paste("Error: More than 1 LD Block associated with", extracted_chr, bp))

    ld_block_dir <- paste0(ld_block_dir, ancestry, "/", extracted_chr, "/", ld_block$start, "_", ld_block$stop)
    study_in_ld_block <- paste0(ld_block_dir, "/", study_dir, "_", extracted_chr, "_", bp, ".z")
    if(!dir.exists(ld_block_dir)) dir.create(ld_block_dir, recursive = T)

    extracted_studies_file <- paste0(ld_block_dir, "/extracted_studies.tsv")
    extracted_studies <- tibble::tribble(~study, ~p_value, ~chr, ~bp,
                                         study_name, p_value_threshold, chr, bp
    )
    if (file.exists(extracted_studies_file)) {
      existing_extracted_studies <- vroom::vroom(extracted_studies_file, show_col_types = F)
      extracted_studies <- dplyr::bind_rows(existing_extracted_studies, extracted_studies)
    }
    vroom::vroom_write(extracted_studies, extracted_studies_file)

    study_file <- paste0(study_dir, ancestry, "_", extracted_chr, "_", bp, ".z")

    file.symlink(study_file, study_in_ld_block)
    return(ld_block)
  }) |> dplyr::bind_rows()
  return(updated_ld_blocks)
})

updated_ld_blocks <- dplyr::bind_rows(updated_ld_blocks) |> dplyr::distinct()

vroom::vroom_write(updated_ld_blocks, paste0(pipeline_metadata_dir, "updated_ld_blocks_to_colocalise.tsv"))
