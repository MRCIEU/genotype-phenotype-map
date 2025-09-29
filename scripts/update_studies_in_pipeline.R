source('../pipeline_steps/constants.R')

ld_blocks <- vroom::vroom('../pipeline_steps/data/ld_blocks.tsv')
ld_info <- construct_ld_block(ld_blocks$ancestry, ld_blocks$chr, ld_blocks$start, ld_blocks$stop)
ld_info <- dplyr::filter(ld_info, dir.exists(ld_block_data))

main <- function() {
  # update_study_dirs()
  update_ld_blocks()
  # update_studies_processed()
}

# TODO: update this method accordingly with how you want to change the data in already ingested studies
update_method <- function(studies_file, imputed_studies_file, ld_block, type) {

  if (!file.exists(studies_file)) {
    print(paste('FILE MISSING:', studies_file))
    return()
  }

  ld_block <- sub('.*EUR', 'EUR', ld_block)
  flattened_ld_block <- gsub('[-/]', '_', ld_block)
  print(paste('flattened_ld_block:', flattened_ld_block))
  imputed_studies <- vroom::vroom(imputed_studies_file, show_col_types = F)
  studies <- vroom::vroom(studies_file, show_col_types = F)
  all_finemapped_files <- lapply(imputed_studies$study, function(study) {
    files <- Sys.glob(glue::glue('{extracted_study_dir}{study}/finemapped/{flattened_ld_block}_*.tsv.gz'))
    return(data.frame(study=study, file=files))
  })
  all_finemapped_files <- dplyr::bind_rows(all_finemapped_files)
  all_finemapped_files <- all_finemapped_files[!grepl('_with_lbf.tsv.gz$', all_finemapped_files$file), ]

  only_missing_files <- dplyr::filter(all_finemapped_files, !file %in% studies$file)
  print(paste('populating', nrow(only_missing_files), 'finemapped files'))

  missing_finemapped_entries <- lapply(seq_len(nrow(only_missing_files)), function(i) {
    missing_entry <- only_missing_files[i, ]
    imputed_info <- dplyr::filter(imputed_studies, study == missing_entry$study) |>
      dplyr::select(study, ancestry, chr, bp, p_value_threshold, sample_size, cis_trans, ld_block, variant_type, category)

    missing_gwas <- vroom::vroom(missing_entry$file, show_col_types = F)
    missing_gwas$Z <- convert_lbf_to_abs_z(missing_gwas$LBF, missing_gwas$SE)
    missing_gwas$P <- 2 * pnorm(-abs(missing_gwas$Z))
    min_p <- min(missing_gwas$P, na.rm = T)
    snp <- missing_gwas$SNP[which.min(missing_gwas$P)]

    svg_file_name <- sub('finemapped', 'svgs/extractions', missing_entry$file)
    svg_file_name <- sub('tsv.gz', 'svg', svg_file_name)
    create_svg_for_ld_block(missing_gwas, missing_entry$study, svg_file_name, ld_block, is_sparse = FALSE)

    imputed_info$svg_file <- svg_file_name
    imputed_info$file <- missing_entry$file
    specific_id <- sub('.*_(\\d+).tsv.gz', '\\1', imputed_info$file)
    imputed_info$unique_study_id <- glue::glue('{imputed_info$study}_{ld_block}_{specific_id}')
    imputed_info$min_p <- min_p
    imputed_info$finemap_message <- NA
    imputed_info$first_finemap_num_results <- NA
    imputed_info$second_finemap_num_results <- NA
    imputed_info$qc_step_run <- NA
    imputed_info$snps_removed_by_qc <- NA
    imputed_info$time_taken <- NA
    imputed_info$ignore <- FALSE
    imputed_info$file_with_lbfs <- glue::glue('{extracted_study_dir}{missing_entry$study}/finemapped/{flattened_ld_block}_with_lbf.tsv.gz')
    return(imputed_info)
  }) |>
    dplyr::bind_rows()

  print(paste('writing', nrow(missing_finemapped_entries), 'finemapped files'))
  studies <- dplyr::bind_rows(studies, missing_finemapped_entries)
  print(paste('now there are', nrow(studies), 'finemapped files'))
  vroom::vroom_write(studies, glue::glue('{studies_file}'))
}

update_study_dirs <- function() {
  print("starting update study dirs")
  ukb_d_dirs <- list.dirs(path = extracted_study_dir, recursive = FALSE, full.names = TRUE)

  # cores_to_use <- 64
  # update_data_in_study_dirs <- parallel::mclapply((X=ukb_d_dirs), mc.cores=cores_to_use, FUN=function(study_dir) {
  result <- lapply(tail(ukb_d_dirs, 1), function(study_dir) {
    print(paste('changing:', study_dir))
    study_id <- basename(study_dir)
    extracted_snps_file <- paste0(study_dir, '/extracted_snps.tsv')
    if (!file.exists(extracted_snps_file )) {
      print(paste('EXTRACTION MISSING: ', extracted_snps_file))
      return()
    }
    # ...
  })
}

update_ld_blocks <- function() {
  blocks <- ld_info$ld_block_data

  cores_to_use <- 10
  update_data_in_ld_blocks <- parallel::mclapply(X=blocks, mc.cores=cores_to_use, FUN=function(ld_block) {
  # update_data_in_ld_blocks <- lapply(blocks, function(ld_block) {
    print(glue::glue('block: {ld_block}\n'))
    # extracted_studies_file <- glue::glue('{ld_block}/extracted_studies.tsv')
    # update_method(extracted_studies_file, type='extracted')


    # standardised_studies_file <- glue::glue('{ld_block}/standardised_studies.tsv')
    # update_method(standardised_studies_file, type='standardised')

    imputed_studies_file <- glue::glue('{ld_block}/imputed_studies.tsv')
    # update_method(imputed_studies_file, type='imputed')

    finemapped_studies_file <- glue::glue('{ld_block}/finemapped_studies.tsv')
    update_method(finemapped_studies_file, imputed_studies_file, ld_block, type='finemapped')

    # coloc_pairwise_results_file <- glue::glue('{ld_block}/coloc_pairwise_results.tsv.gz')
    # update_method(coloc_pairwise_results_file, type='coloc_pairwise')

    # coloc_clustered_results_file <- glue::glue('{ld_block}/coloc_clustered_results.tsv.gz')
    # update_method(coloc_clustered_results_file, type='coloc_clustered')

    # compare_rare_results_file <- glue::glue('{ld_block}/compare_rare_results.tsv')
    # update_method(compare_rare_results_file, type='compare_rare')
  })
}

update_studies_processed <- function() {
  print("starting stupdate_studies_processed")
  studies_to_process_file <- glue::glue('{pipeline_metadata_dir}/studies_to_process.tsv') 
  if (file.exists(studies_to_process_file)) {
    studies_to_process <- vroom::vroom(studies_to_process_file, show_col_types = F)
    #...
    # vroom::vroom_write(studies_to_process, studies_to_process_file)
  }

  studies_processed_file <- glue::glue('{latest_results_dir}/studies_processed.tsv.gz') 
  processed_studies <- vroom::vroom(studies_processed_file, show_col_types = F)
  processed_studies$variant_type <- variant_types$common
  # ...
  # vroom::vroom_write(processed_studies, studies_processed_file)

}

main()