source('../pipeline_steps/constants.R')
options(dplyr.width = Inf)

flippy_flippy <- function(ld_block_matrix, gwas) {
  rsid_match <- match(gwas$RSID, ld_block_matrix$RSID)

  gwas$REF_EA <- ld_block_matrix$EA[rsid_match]
  gwas$REF_OA <- ld_block_matrix$OA[rsid_match]

  to_flip <- gwas$EA == gwas$REF_OA & gwas$OA == gwas$REF_EA
  right_way_around <- gwas$EA == gwas$REF_EA & gwas$OA == gwas$REF_OA
  weird <- nrow(gwas) - sum(to_flip) - sum(right_way_around)
  print(paste('flipping', sum(to_flip), 'eaf values. ', sum(right_way_around), 'are good. ', weird, ' are weird'))
  gwas$EAF[to_flip] <- 1 - gwas$EAF[to_flip]

  gwas <- dplyr::select(gwas, -REF_EA, -REF_OA)

  return(gwas)
}

update_existing_eafs <- function() {
  ld_block_matrices <- list()
  study_to_remove <- '/local-scratch/projects/genotype-phenotype-map/data/study/ebi-a-GCST000758/'

  all_studies <- Sys.glob(paste0(extracted_study_dir, '[a-zC-Z]*/'))
  print(length(all_studies))
  remove_to <- match(study_to_remove, all_studies)
  all_studies <- all_studies[-seq(remove_to)]
  print(length(all_studies))

  for (study in all_studies) {
    extracted_snps_file <- paste0(study, '/extracted_snps.tsv')
      if (!file.exists(extracted_snps_file )) {
        print(paste('EXTRACTION MISSING: ', extracted_snps_file))
        next
      }

    extracted_snps <- vroom::vroom(extracted_snps_file, show_col_types = F)
    if (nrow(extracted_snps) == 0) next

    apply(extracted_snps, 1, function(extraction) {
      ld_region <- extraction[['ld_region']]
      if (ld_region == '//_') return()

      if (is.null(ld_block_matrices[[ld_region]])) {
        ld_block_matrices[[ld_region]] <- vroom::vroom(paste0(ld_block_matrices_dir, ld_region, '.tsv'), show_col_types = F)
      }

      ld_block_matrix <- ld_block_matrices[[ld_region]]

      original_gwas <- vroom::vroom(extraction[['file']], show_col_types = F)

      if (!all(is.na(as.numeric(original_gwas$EAF)))) {
        print(paste('EAF already included: ', extraction[['file']]))
        return()
      }

      imputed_study <- sub('original', 'imputed', extraction[['file']])
      if (!file.exists(imputed_study)) {
        print(paste('IMPUTATION MISSING: ', imputed_study))
        return()
      }
      imputed_gwas <- vroom::vroom(imputed_study, show_col_types = F)
      print(imputed_study)
      updated_imputed_gwas <- flippy_flippy(ld_block_matrix, imputed_gwas)
      vroom::vroom_write(updated_imputed_gwas, imputed_study)

      finemapped_studies <- sub('original', 'finemapped', extraction[['file']])
      finemap_file_prefix <- sub('\\..*', '', finemapped_studies)
      all_finemaps <- Sys.glob(paste0(finemap_file_prefix, '*'))

      if (length(all_finemaps) == 0) {
        print(paste('FINEMAPPING MISSING: ', finemap_file_prefix))
        return()
      }

      for (finemap_study in all_finemaps) {
        finemap_gwas <- vroom::vroom(finemap_study, show_col_types = F)
        print(finemap_study)
        updated_finemap_gwas <- flippy_flippy(ld_block_matrix, finemap_gwas)
        vroom::vroom_write(updated_finemap_gwas, finemap_study)
      }

    })
  }
}

standardise_everything <- function() {
  #ld_matrix_info_files <- Sys.glob(paste0(ld_block_matrices_dir, 'EUR/*/*.tsv'))

  #print('standardising...')
  #for (file in ld_matrix_info_files) {
  #  print(file)
  #  ld_matrix_info <- vroom::vroom(file, show_col_types = F)
  #  ld_matrix_info <- standardise_alleles(ld_matrix_info)
  #  vroom::vroom_write(ld_matrix_info, file)
  #}

  all_studies <- Sys.glob(paste0(extracted_study_dir, '*/'))

  for (study in all_studies) {
    print(paste('standardising', study))
    extracted_snps_file <- paste0(study, '/extracted_snps.tsv')
    extracted_snps <- vroom::vroom(extracted_snps_file, show_col_types = F)
    if (nrow(extracted_snps) == 0) next

    apply(extracted_snps, 1, function(extraction) {
      imputed_study <- sub('original', 'imputed', extraction[['file']])
      if (!file.exists(imputed_study)) {
        print(paste('IMPUTATION MISSING: ', imputed_study))
        return()
      }
      imputed_gwas <- vroom::vroom(imputed_study, show_col_types = F)
      imputed_gwas <- standardise_alleles(imputed_gwas)
      vroom::vroom_write(imputed_gwas, imputed_study)

      finemapped_studies <- sub('original', 'finemapped', extraction[['file']])
      finemap_file_prefix <- sub('\\..*', '', finemapped_studies)
      all_finemaps <- Sys.glob(paste0(finemap_file_prefix, '*'))

      if (length(all_finemaps) == 0) {
        print(paste('FINEMAPPING MISSING: ', finemap_file_prefix))
        return()
      }

      for (finemap_study in all_finemaps) {
        finemap_gwas <- vroom::vroom(finemap_study, show_col_types = F)
        finemap_gwas <- standardise_alleles(finemap_gwas)
        vroom::vroom_write(finemap_gwas, finemap_study)
      }

    })
  }
}

standardise_z_scores <- function() {
  library(parallel)
  all_studies <- Sys.glob(paste0(extracted_study_dir, '*/'))
  study_to_remove <- '/local-scratch/projects/genotype-phenotype-map/data/study/GTEx-sqtl-cis-Cells-Cultured-fibroblasts-chr1-chr1:100921841:100961828:clu-44102:ENSG00000162695-11/'
  remove_to <- match(study_to_remove, all_studies)
  all_studies <- all_studies[-seq(remove_to)]
  print(length(all_studies))

  results <- mclapply(X=all_studies, mc.cores=100, FUN=function(study) {
    print(paste('standardising', study))
    extracted_snps_file <- paste0(study, '/extracted_snps.tsv')
    if (!file.exists(extracted_snps_file )) {
      print(paste('EXTRACTION MISSING: ', extracted_snps_file))
      return() 
    }
    extracted_snps <- vroom::vroom(extracted_snps_file, show_col_types = F)
    if (nrow(extracted_snps) == 0) return()

    apply(extracted_snps, 1, function(extraction) {
      standardised_imputed_study <- sub('original', 'imputed', extraction[['file']])
      old_imputed_study <- sub('study/', 'study_fix/study/', standardised_imputed_study)
      if (!file.exists(standardised_imputed_study)) {
        print(paste('IMPUTATION MISSING: ', standardised_imputed_study))
        return()
      }
      standardised_imputed_gwas <- vroom::vroom(standardised_imputed_study, show_col_types = F)
      old_imputed_gwas <- vroom::vroom(old_imputed_study, show_col_types = F)
      updated_imputed_gwas <- fix_z_score_direction(old_imputed_gwas, standardised_imputed_gwas)

      if (any((updated_imputed_gwas$Z * updated_imputed_gwas$BETA) < 0, na.rm=T)) {
        quit(paste('ERROR: ', standardised_imputed_study, 'BETA and Z still dont match'))
      }

      vroom::vroom_write(updated_imputed_gwas, standardised_imputed_study)

      finemapped_studies <- sub('original', 'finemapped', extraction[['file']])
      finemap_file_prefix <- sub('\\..*', '', finemapped_studies)
      all_finemaps <- Sys.glob(paste0(finemap_file_prefix, '*'))

      if (length(all_finemaps) == 0) {
        print(paste('FINEMAPPING MISSING: ', finemap_file_prefix))
        return()
      }

      for (standardised_finemap_study in all_finemaps) {
        old_finemapped_study <- sub('study/', 'study_fix/study/', standardised_finemap_study)

        standardised_finemap_gwas <- vroom::vroom(standardised_finemap_study, show_col_types = F)
        old_finemap_gwas <- vroom::vroom(old_finemapped_study, show_col_types = F)

        updated_finemap_gwas <- fix_z_score_direction(old_finemap_gwas, standardised_finemap_gwas)
        if (any((updated_finemap_gwas$Z * updated_finemap_gwas$BETA) < 0, na.rm=T)) {
          quit(paste('ERROR: ', standardised_finemap_gwas, 'BETA and Z still dont match'))
        }
        vroom::vroom_write(updated_finemap_gwas, standardised_finemap_study)
      }
    })
  })
}

fix_z_score_direction <- function(old_gwas, standardised_gwas) {
  to_flip <- (old_gwas$EA > old_gwas$OA) & (!old_gwas$EA %in% c("D", "I"))
  print(paste('flipping', sum(to_flip), 'z scores'))
  if (any(to_flip)) {
    standardised_gwas$Z[to_flip] <- -1 * standardised_gwas$Z[to_flip]
  }

  return(standardised_gwas)
}


update_studies_processed <- function() {
  studies_processed_file <- paste0(paste0(results_dir, 'studies_processed.tsv'))
  studies_processed <- vroom::vroom(studies_processed_file, delim='\t', show_col_types=F)
  studies_processed$study_name <- sub('\\.', '-', studies_processed$study_name)
  
  studies_processed$study_name <- sub('BrainMeta-cis-eqtl-BrainMeta-cis-eQTL', 'BrainMeta-cis-eQTL', studies_processed$study_name)
  studies_processed$study_name <- sub('Brain-eMeta-Brain-eMeta', 'Brain-eMeta-full', studies_processed$study_name)

  vroom::vroom_write(studies_processed, studies_processed_file)
}

update_extracted_snps <- function() {
  brain_studies <- Sys.glob(paste0(extracted_study_dir, 'BrainMeta*/'))

  dont_print <- lapply(brain_studies, function(study) {
    extracted_snps_file <- paste0(study, '/extracted_snps.tsv')
    extracted_snps <- vroom::vroom(extracted_snps_file, show_col_types = F)
    if (nrow(extracted_snps) > 0) print(extracted_snps_file)

    extracted_snps$file <- sub('BrainMeta-cis-eqtl-BrainMeta-cis-eQTL', 'BrainMeta-cis-eQTL', extracted_snps$file)
    extracted_snps$file <- sub('Brain-eMeta-Brain-eMeta', 'Brain-eMeta-full', extracted_snps$file)
    extracted_snps$file[grepl('BrainMeta-cis-eQTL', extracted_snps$file)] <- sub('\\.', '-', extracted_snps$file[grepl('BrainMeta-cis-eQTL', extracted_snps$file)])

    vroom::vroom_write(extracted_snps, extracted_snps_file)
  })
}


brain_meta_updates <- function() {
  ld_regions <- vroom::vroom('../pipeline_steps/data/ld_regions.tsv')
  ld_info <- construct_ld_block(ld_regions$ancestry, ld_regions$chr, ld_regions$start, ld_regions$stop)

  lapply(ld_info$ld_block_data, function(ld_block) {
    print(ld_block)
    extracted_studies_file <- paste0(ld_block, '/extracted_studies.tsv')
    if (!file.exists(extracted_studies_file)) return()

    extracted_studies <- vroom::vroom(extracted_studies_file, show_col_types = F)

    extracted_studies$ancestry <- 'EUR'
    extracted_studies$study <- sub('BrainMeta-cis-eqtl-BrainMeta-cis-eQTL', 'BrainMeta-cis-eQTL', extracted_studies$study)
    extracted_studies$study <- sub('Brain-eMeta-Brain-eMeta', 'Brain-eMeta-full', extracted_studies$study)
    extracted_studies$study <- sub('\\.', '-', extracted_studies$study)
    extracted_studies$file <- sub('BrainMeta-cis-eqtl-BrainMeta-cis-eQTL', 'BrainMeta-cis-eQTL', extracted_studies$file)
    extracted_studies$file <- sub('Brain-eMeta-Brain-eMeta', 'Brain-eMeta-full', extracted_studies$file)
    extracted_studies$file[grepl('BrainMeta-cis-eQTL', extracted_studies$file)] <- sub('\\.', '-', extracted_studies$file[grepl('BrainMeta-cis-eQTL', extracted_studies$file)])

    vroom::vroom_write(extracted_studies, extracted_studies_file)


    imputed_studies_file <- paste0(ld_block, '/imputed_studies.tsv')
    if (!file.exists(imputed_studies_file)) return()

    imputed_studies <- vroom::vroom(imputed_studies_file, show_col_types = F)

    imputed_studies$ancestry <- 'EUR'
    imputed_studies$study <- sub('BrainMeta-cis-eqtl-BrainMeta-cis-eQTL', 'BrainMeta-cis-eQTL', imputed_studies$study)
    imputed_studies$study <- sub('Brain-eMeta-Brain-eMeta', 'Brain-eMeta-full', imputed_studies$study)
    imputed_studies$study <- sub('\\.', '-', imputed_studies$study)
    imputed_studies$file <- sub('BrainMeta-cis-eqtl-BrainMeta-cis-eQTL', 'BrainMeta-cis-eQTL', imputed_studies$file)
    imputed_studies$file <- sub('Brain-eMeta-Brain-eMeta', 'Brain-eMeta-full', imputed_studies$file)
    imputed_studies$file[grepl('BrainMeta-cis-eQTL', imputed_studies$file)] <- sub('\\.', '-', imputed_studies$file[grepl('BrainMeta-cis-eQTL', imputed_studies$file)])

    vroom::vroom_write(imputed_studies, imputed_studies_file)

    finemapped_studies_file <- paste0(ld_block, '/finemapped_studies.tsv')
    if (!file.exists(finemapped_studies_file)) return()

    finemapped_studies <- vroom::vroom(finemapped_studies_file, show_col_types = F)

    finemapped_studies$ancestry <- 'EUR'
    finemapped_studies$study <- sub('BrainMeta-cis-eqtl-BrainMeta-cis-eQTL', 'BrainMeta-cis-eQTL', finemapped_studies$study)
    finemapped_studies$study <- sub('Brain-eMeta-Brain-eMeta', 'Brain-eMeta-full', finemapped_studies$study)
    finemapped_studies$study <- sub('\\.', '-', finemapped_studies$study)
    finemapped_studies$unique_study_id <- sub('BrainMeta-cis-eqtl-BrainMeta-cis-eQTL', 'BrainMeta-cis-eQTL', finemapped_studies$unique_study_id)
    finemapped_studies$unique_study_id <- sub('Brain-eMeta-Brain-eMeta', 'Brain-eMeta-full', finemapped_studies$unique_study_id)
    finemapped_studies$unique_study_id <- sub('\\.', '-', finemapped_studies$unique_study_id)
    finemapped_studies$unique_study_id <- gsub(' ', '', finemapped_studies$unique_study_id)
    fill_study_name <- grepl('^EUR', finemapped_studies$unique_study_id)
    finemapped_studies$unique_study_id[fill_study_name] <- paste0(finemapped_studies$study[fill_study_name], '_', finemapped_studies$unique_study_id[fill_study_name])

    finemapped_studies$file <- sub('BrainMeta-cis-eqtl-BrainMeta-cis-eQTL', 'BrainMeta-cis-eQTL', finemapped_studies$file)
    finemapped_studies$file <- sub('Brain-eMeta-Brain-eMeta', 'Brain-eMeta-full', finemapped_studies$file)

    chr_bp_version <- sub('^[^_]*_', '', finemapped_studies$unique_study_id[grepl('BrainMeta-cis-eQTL', finemapped_studies$unique_study_id)])
    move_to <- paste0(extracted_study_dir, finemapped_studies$study[grepl('BrainMeta-cis-eQTL', finemapped_studies$unique_study_id)], '/finemapped/EUR_', chr_bp_version, '.tsv.gz')
    file.copy(finemapped_studies$file[grepl('BrainMeta-cis-eQTL', finemapped_studies$file)], move_to)

    finemapped_studies$file[grepl('BrainMeta-cis-eQTL', finemapped_studies$file)] <- move_to

    fill_eur <- !grepl('EUR', finemapped_studies$unique_study_id)
    finemapped_studies$unique_study_id[fill_eur] <- sub('_', '_EUR_', finemapped_studies$unique_study_id[fill_eur])

    no_min_p_finemapped <- finemapped_studies[is.na(finemapped_studies$min_p), ]
    actual_min_ps <- lapply(no_min_p_finemapped$file, function(file) {
      return(min(vroom::vroom(file, show_col_types = F)$P))
    })
    finemapped_studies$min_p[is.na(finemapped_studies$min_p)] <- unlist(actual_min_ps)

    vroom::vroom_write(finemapped_studies, finemapped_studies_file)
  })

}

update_imputed_studies <- function() {
  ld_regions <- vroom::vroom('../pipeline_steps/data/ld_regions.tsv')
  ld_info <- construct_ld_block(ld_regions$ancestry, ld_regions$chr, ld_regions$start, ld_regions$stop)

  lapply(ld_info$ld_block_data, function(ld_block) {
    extracted_studies_file <- paste0(ld_block, '/extracted_studies.tsv')
    imputed_studies_file <- paste0(ld_block, '/imputed_studies.tsv')
    if (!file.exists(extracted_studies_file)) return()

    extracted_studies <- vroom::vroom(extracted_studies_file, show_col_types = F)
    imputed_studies <- vroom::vroom(imputed_studies_file, col_types = vroom::cols(
      chr = vroom::col_character(),
      bp = vroom::col_number(),
      p_value_threshold = vroom::col_number(),
      eaf_from_reference_panel = vroom::col_logical(),
      sample_size = vroom::col_number()
      ), show_col_types = F)
    if (nrow(extracted_studies) == nrow(imputed_studies)) return()
    print(paste(ld_block, 'sizes are different...', nrow(extracted_studies) - nrow(imputed_studies)))

    add_to_imputed <- apply(extracted_studies, 1, function(extracted_study) {
      already_in_imputed <- dplyr::filter(imputed_studies, study == extracted_study['study'])
      if (nrow(already_in_imputed) > 0) return()

      imputed_study <- sub('original', 'imputed', extracted_study[['file']])
      if (!file.exists(imputed_study)) return()

      eaf_from_reference_panel <- ifelse(grepl('GTEx', extracted_study['study']), TRUE, NA)

      return(data.frame(study=extracted_study['study'], file=extracted_study['file'], chr=extracted_study['chr'], bp=as.numeric(extracted_study['bp']),
                        p_value_threshold=as.numeric(extracted_study['p_value_threshold']), sample_size=as.numeric(extracted_study['sample_size']), cis_trans=extracted_study['cis_trans'],
                        rows_imputed=NA, eaf_from_reference_panel=eaf_from_reference_panel
                        )
      )

    }) |> dplyr::bind_rows()

    imputed_studies <- dplyr::bind_rows(imputed_studies, add_to_imputed) |> dplyr::distinct()
    vroom::vroom_write(imputed_studies, imputed_studies_file)
  })

}

update_finemapped_studies <- function() {
  ld_regions <- vroom::vroom('../pipeline_steps/data/ld_regions.tsv')
  ld_regions <- ld_regions[ld_regions$chr >= 12, ]
  ld_info <- construct_ld_block(ld_regions$ancestry, ld_regions$chr, ld_regions$start, ld_regions$stop)

  lapply(ld_info$ld_block_data, function(ld_block) {
    finemapped_studies_file <- paste0(ld_block, '/finemapped_studies.tsv')
    imputed_studies_file <- paste0(ld_block, '/imputed_studies.tsv')
    if (!file.exists(finemapped_studies_file)) return()
    print(ld_block)

    imputed_studies <- vroom::vroom(imputed_studies_file, show_col_types=F) 
    finemapped_studies <- vroom::vroom(finemapped_studies_file, col_types = vroom::cols(
      chr = vroom::col_character(),
      bp = vroom::col_number(),
      min_p = vroom::col_number(),
      p_value_threshold = vroom::col_number(),
      sample_size = vroom::col_number()
      ), show_col_types = F)

    if (nrow(imputed_studies) == 0) return()
    add_to_finemapped <- apply(imputed_studies, 1, function(imputed_study) {
      finemap_file_prefix <- sub('imputed', 'finemapped', imputed_study[['file']])
      finemap_file_prefix <- sub('\\..*', '', finemap_file_prefix)
      finemap_files <- Sys.glob(paste0(finemap_file_prefix, '*'))

      if (length(finemap_files) == 0) return()
      if (nrow(dplyr::filter(finemapped_studies, grepl(finemap_file_prefix, file))) > 0) return()

      unique_study_ids <- file_prefix(finemap_files)
      message <- if(length(finemap_files) > 1) 'success' else NA

      bps <- lapply(finemap_files, function(file) {
        finemap <- vroom::vroom(file, show_col_types=F)
        min <- finemap[which.min(finemap$P), ]
        return(data.frame(bp=min$BP, min_p=min$P))
      }) |> dplyr::bind_rows()

      print(unique_study_ids)
      print(finemap_files)
      print(bps$bp)
      print(bps$min_p)

      return(data.frame(study=imputed_study['study'], unique_study_id=unique_study_ids, file=finemap_files, chr=imputed_study['chr'], bp=as.numeric(bps$bp),
                        p_value_threshold=as.numeric(imputed_study['p_value_threshold']), sample_size=as.numeric(imputed_study['sample_size']), cis_trans=imputed_study['cis_trans'],
                        message=message, min_p=bps$min_p
                        )
      )


    }) |> dplyr::bind_rows()

    if (nrow(add_to_finemapped) == 0) return()
    print(paste('missing imputed from finemapped',nrow(add_to_finemapped)))

    finemapped_studies <- dplyr::bind_rows(finemapped_studies, add_to_finemapped) |> dplyr::distinct()
    vroom::vroom_write(finemapped_studies, finemapped_studies_file)
  })
}

populate_json_for_besd <- function() {
  gtex <- vroom::vroom('/local-scratch/projects/genotype-phenotype-map/GTEx_portal.csv')

  gtex$Tissue <- gsub('[ -]', '_', gtex$Tissue)
  gtex$Tissue <- gsub('_+', '_', gtex$Tissue)
  gtex$Tissue <- gsub('[()]', '', gtex$Tissue)

  apply(gtex, 1, function(tissue) {
    metadata <- list(sample_size=as.numeric(tissue['# RNASeq and Genotyped samples']), ancestry='EUR', cis_trans='cis', category='continuous')
    metadata_json <- jsonlite::toJSON(metadata, auto_unbox = T, pretty = T)
    file_name <- paste0('/local-scratch/data/GTEx-cis/', tissue['Tissue'], '.json')
    write(metadata_json, file_name)
  })

}

remove_imputed_and_finemapped_results_from_pipeline <- function(study_pattern) {
  NUM_PARALLEL_JOBS <- 100

  all_studies <- Sys.glob(paste0(extracted_study_dir, study_pattern, '*/'))
  # study_to_remove <- '/local-scratch/projects/genotype-phenotype-map/data/study/study-name'
  # remove_to <- match(study_to_remove, all_studies)
  # all_studies <- all_studies[-seq(remove_to)]

  results <- parallel::mclapply(X=all_studies, mc.cores=NUM_PARALLEL_JOBS, FUN=function(study_dir) {
    print(paste('changing:', study_dir))
    study_id <- basename(study_dir)
    extracted_snps_file <- paste0(study_dir, '/extracted_snps.tsv')
    if (!file.exists(extracted_snps_file )) {
      print(paste('EXTRACTION MISSING: ', extracted_snps_file))
      return()
    }
    extracted_snps <- vroom::vroom(extracted_snps_file, show_col_types = F)
    if (nrow(extracted_snps) == 0) return()

    apply(extracted_snps, 1, function(extraction) {
      ld_block <- paste0(ld_block_data_dir, extraction[['ld_region']])

      imputed_studies_file <- paste0(ld_block, '/imputed_studies.tsv')
      if (!file.exists(imputed_studies_file)) {
        print(paste('IMPUTED STUDIES FILE MISSING:', imputed_studies_file))
        return()
      }

      imputed_studies <- vroom::vroom(imputed_studies_file, show_col_types = F)
      entries <- nrow(imputed_studies)
      imputed_studies <- dplyr::filter(imputed_studies, !study == study_id)
      print(paste('removed', entries - nrow(imputed_studies), 'rows from imputed_studies'))
      vroom::vroom_write(imputed_studies, imputed_studies_file)
      imputed_files <- Sys.glob(paste0(study_dir, '/imputed/*'))
      file.remove(imputed_files)


      finemapped_studies_file <- paste0(ld_block, '/finemapped_studies.tsv')
      if (!file.exists(finemapped_studies_file)) {
        print(paste('FINEMAPPED STUDIES FILE MISSING:', finemapped_studies_file))
        return()
      }

      finemapped_studies <- vroom::vroom(finemapped_studies_file, show_col_types = F)
      entries <- nrow(finemapped_studies)
      finemapped_studies <- dplyr::filter(finemapped_studies, !study == study_id)
      print(paste('removed', entries - nrow(finemapped_studies), 'rows from finemapped_studies'))
      vroom::vroom_write(finemapped_studies, finemapped_studies_file)
    })
    finemapped_files <- Sys.glob(paste0(study_dir, '/finemapped/*'))
    file.remove(finemapped_files)
  })
}

print_alleles_to_flip <- function() {
  flipped_dir <- paste0(thousand_genomes_dir, '/flipped/')
  bims <- Sys.glob(paste0(flipped_dir, '*.bim'))

  lapply(bims, function(bim_file) {
    bim <- vroom::vroom(bim_file, col_names=F)
    print(names(bim))
    bim <- bim[(bim$X5 > bim$X6), ]
    vroom::vroom_write(dplyr::select(bim, X2), paste0(bim_file, '_toflip'), col_names=F)
  })

}

print_alleles_to_flip()
