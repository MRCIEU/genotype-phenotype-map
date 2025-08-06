source('../pipeline_steps/constants.R')
source('../pipeline_steps/common_extraction_functions.R')
library(parallel)

ld_blocks <- vroom::vroom('../pipeline_steps/data/ld_blocks.tsv')
ld_info <- construct_ld_block(ld_blocks$ancestry, ld_blocks$chr, ld_blocks$start, ld_blocks$stop)
ld_info <- dplyr::filter(ld_info, dir.exists(ld_block_data))

main <- function() {
  # update_study_dirs()
  update_ld_blocks()
  # update_studies_processed()
}

# snp_annotations_file <- glue::glue('{variant_annotation_dir}/vep_annotations_hg38.tsv.gz')
# snp_annotations <- vroom::vroom(snp_annotations_file, show_col_types = F) |>
  # dplyr::select(snp, ref_panel_snp)

missing_snps_in_study_extractions <- vroom::vroom('/local-scratch/projects/genotype-phenotype-map/results/current/missing_snps_in_study_extractions.tsv')

standardise_alleles <- function(gwas) {
  gwas$EA <- toupper(as.character(gwas$EA))
  gwas$OA <- toupper(as.character(gwas$OA))

  gwas$EA <- ifelse(gwas$EA == T | gwas$EA == 'TRUE', 'T', gwas$EA)
  gwas$OA <- ifelse(gwas$OA == T | gwas$OA == 'TRUE', 'T', gwas$OA)

  columns_to_coerce <- c('EAF', 'BETA', 'SE')
  gwas <- dplyr::mutate(gwas, dplyr::across(dplyr::all_of(columns_to_coerce), as.numeric))

  to_flip <- (gwas$EA > gwas$OA) & (!gwas$EA %in% c("D", "I"))
  if (any(to_flip)) {
    if ('EAF' %in% names(gwas)) {
      gwas$EAF[to_flip] <- 1 - gwas$EAF[to_flip]
    }
    if ('BETA' %in% names(gwas)) {
      gwas$BETA[to_flip] <- -1 * gwas$BETA[to_flip]
    }
    if ('Z' %in% names(gwas)) {
      gwas$Z[to_flip] <- -1 * gwas$Z[to_flip]
    }

    temp <- gwas$OA[to_flip]
    gwas$OA[to_flip] <- gwas$EA[to_flip]
    gwas$EA[to_flip] <- temp
  }

  gwas$SNP <- format_unique_snp_string(gwas$CHR, gwas$BP, gwas$EA, gwas$OA)
  return(gwas)
}

# TODO: update this method accordingly with how you want to change the data in already ingested studies
update_method <- function(studies_file, ld_block, type) {
  if (!file.exists(studies_file)) {
    print(paste('FILE MISSING:', studies_file))
    return()
  }
  imputed_studies <- vroom::vroom(glue::glue('{ld_block}/imputed_studies.tsv'), show_col_types = F, altrep = F)
  finemapped_studies <- vroom::vroom(studies_file, show_col_types = F, altrep = F)

  if (type == 'finemapped') {
    lapply(seq_len(nrow(finemapped_studies)), function(i) {
      finemapped_study <- finemapped_studies[i, , drop = FALSE]
      gwas_with_lbf <- vroom::vroom(finemapped_study$file_with_lbfs, show_col_types = F, altrep = F)

      if (!"P" %in% names(gwas_with_lbf)) {
        imputed_gwas_file <- dplyr::filter(imputed_studies, study == finemapped_study$study)$file

        if (length(imputed_gwas_file) == 0) {
          print(paste('Imputed gwas entry missing, deleting finemapped study:', finemapped_study$file_with_lbfs))
          # finemapped_studies <- finemapped_studies[-i, ]
          # vroom::vroom_write(finemapped_studies, studies_file)
          return()
        }

        message(paste('Adding P values to', finemapped_study$file_with_lbfs, 'from', imputed_gwas_file))

        imputed_gwas <- vroom::vroom(imputed_gwas_file, show_col_types = F) |>
          dplyr::select(SNP, P)
        gwas_with_lbf <- dplyr::left_join(gwas_with_lbf, imputed_gwas, by = 'SNP')
        vroom::vroom_write(gwas_with_lbf, finemapped_study$file_with_lbfs)
      }
    })
  }

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
  # blocks <- blocks[20:length(blocks)]

  cores_to_use <- 20
  # update_data_in_ld_blocks <- parallel::mclapply(X=blocks, mc.cores=cores_to_use, FUN=function(ld_block) {
  update_data_in_ld_blocks <- lapply(blocks, function(ld_block) {
    print(glue::glue('block: {ld_block}\n'))
    # extracted_studies_file <- glue::glue('{ld_block}/extracted_studies.tsv')
    # update_method(extracted_studies_file, type='extracted')

    # standardised_studies_file <- glue::glue('{ld_block}/standardised_studies.tsv')
    # update_method(standardised_studies_file, type='standardised')

    # imputed_studies_file <- glue::glue('{ld_block}/imputed_studies.tsv')
    # update_method(imputed_studies_file, type='imputed')

    finemapped_studies_file <- glue::glue('{ld_block}/finemapped_studies.tsv')
    update_method(finemapped_studies_file, ld_block, type='finemapped')

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
  vroom::vroom_write(processed_studies, studies_processed_file)

}

main()