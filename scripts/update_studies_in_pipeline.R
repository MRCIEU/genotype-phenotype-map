source("../pipeline_steps/constants.R")

ld_blocks <- vroom::vroom("../pipeline_steps/data/ld_blocks.tsv")
ld_info <- construct_ld_block(ld_blocks$ancestry, ld_blocks$chr, ld_blocks$start, ld_blocks$stop)
ld_info <- dplyr::filter(ld_info, dir.exists(ld_block_data))

main <- function() {
  # update_study_dirs()
  update_ld_blocks()
  # update_studies_processed()
  return()
}

# TODO: update this method accordingly with how you want to change the data in already ingested studies
update_method <- function(studies_file, type, unique_study_ids_to_update = NULL) {
  if (!file.exists(studies_file)) {
    print(paste("FILE MISSING:", studies_file))
    return()
  }
  # ...
}

update_study_dirs <- function() {
  print("starting update study dirs")
  ukb_d_dirs <- list.dirs(path = extracted_study_dir, recursive = FALSE, full.names = TRUE)

  # cores_to_use <- 64
  # update_data_in_study_dirs <- parallel::mclapply((X=ukb_d_dirs), mc.cores=cores_to_use, FUN=function(study_dir) {
  result <- lapply(tail(ukb_d_dirs, 1), function(study_dir) {
    print(paste("changing:", study_dir))
    study_id <- basename(study_dir)
    extracted_snps_file <- paste0(study_dir, "/extracted_snps.tsv")
    if (!file.exists(extracted_snps_file)) {
      print(paste("EXTRACTION MISSING: ", extracted_snps_file))
      return()
    }
    # ...
    return()
  })
  return()
}

update_ld_blocks <- function() {
  source("../pipeline_steps/gwas_calculations.R")
  # snp_annotations <- vroom::vroom(file.path(variant_annotation_dir, "vep_annotations_hg38.tsv.gz"), show_col_types =  F) |> # nolint: line_length_linter.
  # dplyr::select(chr, bp, snp)
  blocks <- ld_info$ld_block_data

  cores_to_use <- 15
  update_data_in_ld_blocks <- parallel::mclapply(X = blocks, mc.cores = cores_to_use, FUN = function(ld_block) {
    # update_data_in_ld_blocks <- lapply(blocks, function(ld_block) {
    print(glue::glue("block: {ld_block}\n"))
    return()
    # extracted_studies_file <- glue::glue('{ld_block}/extracted_studies.tsv')
    # update_method(extracted_studies_file, type='extracted')

    # standardised_studies_file <- glue::glue('{ld_block}/standardised_studies.tsv')
    # update_method(standardised_studies_file, type='standardised')

    # imputed_studies_file <- glue::glue('{ld_block}/imputed_studies.tsv')
    # update_method(imputed_studies_file, type='imputed')

    # finemapped_studies_file <- glue::glue('{ld_block}/finemapped_studies.tsv')
    # update_method(finemapped_studies_file, type='finemapped')

    # coloc_pairwise_results_file <- glue::glue('{ld_block}/coloc_pairwise_results.tsv.gz')
    # update_method(coloc_pairwise_results_file, type='coloc_pairwise', unique_study_ids_to_update = unique_study_ids_to_update) # nolint: line_length_linter.

    # coloc_clustered_results_file <- glue::glue('{ld_block}/coloc_clustered_results.tsv.gz')
    # update_method(coloc_clustered_results_file, type='coloc_clustered')

    # compare_rare_results_file <- glue::glue('{ld_block}/compare_rare_results.tsv')
    # update_method(compare_rare_results_file, type='compare_rare')
  })
  return()
}

update_studies_processed <- function() {
  print("starting stupdate_studies_processed")
  studies_to_process_file <- glue::glue("{pipeline_metadata_dir}/studies_to_process.tsv")
  if (file.exists(studies_to_process_file)) {
    studies_to_process <- vroom::vroom(studies_to_process_file, show_col_types = F)
    # ...
    # vroom::vroom_write(studies_to_process, studies_to_process_file)
  }

  studies_processed_file <- glue::glue("{latest_results_dir}/studies_processed.tsv.gz")
  processed_studies <- vroom::vroom(studies_processed_file, show_col_types = F)
  processed_studies$variant_type <- variant_types$common
  # ...
  # vroom::vroom_write(processed_studies, studies_processed_file)
  return()
}

main()
