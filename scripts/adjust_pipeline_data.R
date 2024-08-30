source('../pipeline_steps/constants.R')

remove_studies_from_pipeline <- function(study_pattern) {
  NUM_PARALLEL_JOBS <- 100
  all_studies <- Sys.glob(paste0(extracted_study_dir, study_pattern, '*/'))
  results <- parallel::mclapply(X=all_studies, mc.cores=NUM_PARALLEL_JOBS, FUN=function(study) {
    print(paste('changing:', study))
    extracted_snps_file <- paste0(study, '/extracted_snps.tsv')
    if (!file.exists(extracted_snps_file )) {
      print(paste('Extraction missing: ', extracted_snps_file))
      return()
    }
    extracted_snps <- vroom::vroom(extracted_snps_file, show_col_types = F)
    if (nrow(extracted_snps) == 0) return()

    apply(extracted_snps, 1, function(extraction) {
      ld_block <- extraction[['ld_region']]
      print('actually doing')

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
  })
}

cleanup_empty_dirs <- function(study_pattern) {
  NUM_PARALLEL_JOBS <- 100
  all_studies <- Sys.glob(paste0(extracted_study_dir, study_pattern, '*/'))
  results <- parallel::mclapply(X=all_studies, mc.cores=NUM_PARALLEL_JOBS, FUN=function(study) {
    print(paste('changing:', study))
    extracted_snps_file <- paste0(study, '/extracted_snps.tsv')
    if (!file.exists(extracted_snps_file )) {
      print(paste('Extraction missing: ', extracted_snps_file))
      return()
    }
    extracted_snps <- vroom::vroom(extracted_snps_file, show_col_types = F)
    if (nrow(extracted_snps) == 0) {
      file.remove(paste0(study, 'original'), recursive=T)
      file.remove(paste0(study, 'imputed'), recursive=T)
      file.remove(paste0(study, 'finemapped'), recursive=T)
    }
  })
}

update_extracted_studies <- function() {
  NUM_PARALLEL_JOBS <- 200
  ld_regions <- vroom::vroom('../pipeline_steps/data/ld_regions.tsv')
  ld_info <- construct_ld_block(ld_regions$ancestry, ld_regions$chr, ld_regions$start, ld_regions$stop)

  results <- parallel::mclapply(X=ld_info$ld_block_data, mc.cores=NUM_PARALLEL_JOBS, FUN=function(ld_block) {
    extracted_studies_file <- paste0(ld_block, '/extracted_studies.tsv')
    if (!file.exists(extracted_studies_file)) return()

    extracted_studies <- vroom::vroom(extracted_studies_file, show_col_types = F) |>
      dplyr::select(-data_type)
    extracted_studies$reference_build <- reference_builds$GRCh37

    vroom::vroom_write(extracted_studies, extracted_studies_file)
  })
}

update_ld_matrix_data <- function() {
  block_matrix_dirs <- paste0(ld_reference_panel_dir, 'EUR/', seq(1,22))
  for (dir in block_matrix_dirs) {
    setwd(dir)
    files <- list.files()
    new_files <- lapply(files, function(old_file) {
      new_file <- sub('-', '_', old_file)
      new_bp <- as.numeric(sub('.*_(\\d+)\\..*', '\\1', new_file))
      new_bp <- new_bp - 1
      new_file <- sub('(.*_)\\d+(\\..*)', paste0('\\1', new_bp, '\\2'), new_file)
      file.rename(old_file, new_file)
      return(new_file)
    })

    freq_files <- Filter(\(file) grepl('afreq$', file), files)
    bim_files <- Filter(\(file) grepl('bim$', file), files)
    if (length(freq_files) != length(bim_files)) {
      stop('bad')
    }
    mapply(function(freq_file, bim_file) {
      bim <- vroom::vroom(bim_file, delim='\t', col_names = F, show_col_types = F) |>
        dplyr::select(X2, X4)
      freq <- vroom::vroom(freq_file, show_col_types = F) |>
        dplyr::rename(SNP='ID', CHR='#CHROM', EA='ALT', OA='REF', EAF='ALT_FREQS') |>
        dplyr::select(SNP, CHR, EA, OA, EAF)
      freq$BP <- bim$X4[match(freq$SNP, bim$X2)]

      new_file_name <- sub('bim', 'tsv', bim_file)
      vroom::vroom_write(freq, new_file_name)
      print(new_file_name)
    },freq_files, bim_files)
  }
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
      ld_block <- extraction[['ld_region']]
      print('actually doing')

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

skip_steps <- function() {
  ld_blocks <- vroom::vroom('data/ld_regions.tsv', show_col_types=F)
  vroom::vroom_write(data.frame(), paste0(pipeline_metadata_dir, 'updated_ld_blocks_to_colocalise.tsv'))

  for (i in 1:nrow(ld_blocks)) {
    ld_block <- ld_blocks[i, ]
    ld_block_dir <- paste0(ld_block_data_dir, ld_block$ancestry, '/', ld_block$chr, '/', ld_block$start, '_', ld_block$stop)
    vroom::vroom_write(data.frame(), paste0(ld_block_dir, '/imputation_complete'))
    vroom::vroom_write(data.frame(), paste0(ld_block_dir, '/finemapping_complete'))
  }
}


