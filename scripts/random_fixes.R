source("../pipeline_steps/constants.R")
options(dplyr.width = Inf)

update_files_with_lbfs <- function() {
  ld_blocks <- vroom::vroom("../pipeline_steps/data/ld_blocks.tsv")
  ld_info <- construct_ld_block(ld_blocks$ancestry, ld_blocks$chr, ld_blocks$start, ld_blocks$stop)
  ld_info <- ld_info[dir.exists(ld_info$ld_block_data), ]
  blocks <- ld_info$block

  parallel::mclapply(blocks, mc.cores = 30, function(block) {
    tryCatch({
      finemapped_studies <- vroom::vroom(
        glue::glue("{ld_block_data_dir}/{block}/finemapped_studies.tsv"), show_col_types = F
      ) |>
        dplyr::filter(file.info(file_with_lbfs)$mtime < "2025-11-03")

      imputed_studies <- vroom::vroom(glue::glue("{ld_block_data_dir}/{block}/imputed_studies.tsv"), show_col_types = F)
      studies_to_update <- unique(finemapped_studies$study)

      safe_lapply(studies_to_update, function(imputed_study) {
        this_imputed_study <- dplyr::filter(imputed_studies, study == imputed_study)
        if (nrow(this_imputed_study) != 1) stop(glue::glue("Imputed study not found: {imputed_study} in {block}"))
        print(glue::glue("Updating {this_imputed_study$file}"))
        if (file.info(this_imputed_study$file)$size == 0) {
          return()
        }
        imputed_gwas_with_lbfs <- vroom::vroom(this_imputed_study$file,
          show_col_types = F,
          col_select = c("SNP", "CHR", "BP", "EA", "OA", "EAF", "Z", "BETA", "SE", "P")
        )
        if (!"IMPUTED" %in% colnames(imputed_gwas_with_lbfs)) {
          imputed_gwas_with_lbfs$IMPUTED <- FALSE
        }
        original_nrow <- nrow(imputed_gwas_with_lbfs)

        finemapped_for_study <- dplyr::filter(finemapped_studies, study == imputed_study)
        file_with_lbf_name <- finemapped_for_study$file_with_lbfs[1]
        for (i in seq_len(nrow(finemapped_for_study))) {
          this_finemap_study <- finemapped_for_study[i, ]
          finemap_number <- as.numeric(stringr::str_extract(this_finemap_study$unique_study_id, "(?<=_)\\d+$"))
          lbf_column_name <- glue::glue("LBF_{finemap_number}")
          finemapped_gwas <- vroom::vroom(this_finemap_study$file, show_col_types = F) |>
            dplyr::select(SNP, LBF)
          imputed_gwas_with_lbfs <- dplyr::left_join(imputed_gwas_with_lbfs, finemapped_gwas, by = "SNP") |>
            dplyr::rename(!!lbf_column_name := LBF)
        }
        if (any(is.na(imputed_gwas_with_lbfs$SNP) |
                is.na(imputed_gwas_with_lbfs$BETA) |
                is.na(imputed_gwas_with_lbfs$SE)
            ) | nrow(imputed_gwas_with_lbfs) != original_nrow) {
          stop(glue::glue("Missing SNP, BETA, or SE in {file_with_lbf_name} for {imputed_study} in {block}"))
        }
        vroom::vroom_write(imputed_gwas_with_lbfs, file_with_lbf_name)
        return()
      })
      return(invisible(NULL))
    },
    error = function(e) {
      print(glue::glue("Error in {block}: {e}"))
      stop(glue::glue("Error in {block}: {e}"))
    })
    return(invisible(NULL))
  })
  return()
}

#' Populate BETA column in finemapped summary stats files from file_with_lbfs.
#' The 'file' in finemapped_studies.tsv is a tsv of summary stats that may lack BETA;
#' this function adds BETA from the corresponding file_with_lbfs for each SNP and saves.
populate_beta_in_finemapped_files <- function() {
  ld_blocks <- vroom::vroom("../pipeline_steps/data/ld_blocks.tsv")
  ld_info <- construct_ld_block(ld_blocks$ancestry, ld_blocks$chr, ld_blocks$start, ld_blocks$stop)
  ld_info <- ld_info[dir.exists(ld_info$ld_block_data), ]
  blocks <- ld_info$block

  resolve_path <- function(path) {
    if (file.exists(path)) return(path)
    full_path <- file.path(data_dir, sub("^/", "", path))
    if (file.exists(full_path)) return(full_path)
    return(path)
  }

  parallel::mclapply(blocks, mc.cores = 30, function(block) {
    tryCatch({
      finemapped_studies <- vroom::vroom(
        glue::glue("{ld_block_data_dir}/{block}/finemapped_studies.tsv"),
        show_col_types = F
      )
      message(glue::glue("Processing {nrow(finemapped_studies)} finemapped studies in {block}"))

      for (i in seq_len(nrow(finemapped_studies))) {
        row <- finemapped_studies[i, ]
        if (is.na(row$file_with_lbfs) || is.na(row$file)) next
        file_path <- resolve_path(row$file)
        file_with_lbfs_path <- resolve_path(row$file_with_lbfs)

        if (!file.exists(file_path) || !file.exists(file_with_lbfs_path)) {
          message(glue::glue("Skipping {row$unique_study_id}: file or file_with_lbfs not found"))
          next
        }

        finemapped_gwas <- vroom::vroom(file_path, show_col_types = F)
        if ("BETA" %in% colnames(finemapped_gwas)) {
          next
        }

        lbfs_gwas <- vroom::vroom(
          file_with_lbfs_path,
          show_col_types = F,
          col_select = c("SNP", "BETA")
        )

        finemapped_gwas <- finemapped_gwas |>
          dplyr::left_join(lbfs_gwas, by = "SNP") |>
          dplyr::select(dplyr::any_of(c("SNP", "CHR", "BP", "BETA", "SE", "EAF", "IMPUTED", "LBF")), dplyr::everything())

        vroom::vroom_write(finemapped_gwas, file_path)
      }
      return(invisible(NULL))
    }, error = function(e) {
      message(glue::glue("Error in {block}: {e$message}"))
    })
  })
  return()
}

fix_missing_bp_values <- function() {
  ld_blocks <- vroom::vroom("../pipeline_steps/data/ld_blocks.tsv")
  ld_info <- construct_ld_block(ld_blocks$ancestry, ld_blocks$chr, ld_blocks$start, ld_blocks$stop)
  ld_info <- ld_info[dir.exists(ld_info$ld_block_data), ]
  blocks <- ld_info$block

  block <- "EUR/6/31282730-32518958"
  imputed_snps_to_remove <- vroom::vroom(
    glue::glue("{ld_block_data_dir}/{block}/imputation_snps_to_remove.tsv"),
    show_col_types = F
  ) |>
    dplyr::filter(!is.na(bps_to_remove))
  finemapped_snps_to_remove <- vroom::vroom(
    glue::glue("{ld_block_data_dir}/{block}/problematic_finemapped_snps.tsv"),
    show_col_types = F
  )
  orig_finemap <- vroom::vroom("../orig_finemapped_studies.tsv", show_col_types = F) |>
    dplyr::mutate(coverage = dplyr::case_when(
      grepl("godmc", study) ~ coverage_types$sparse,
      TRUE ~ coverage_types$dense
    ))
  current_finemap <- vroom::vroom(glue::glue("{ld_block_data_dir}/{block}/finemapped_studies.tsv"), show_col_types = F)
  updated_finemap <- current_finemap
  updated_finemap$p_value_threshold <- 1.5e-4
  updated_finemap$bp <- stringr::str_extract(updated_finemap$snp, "(?<=:)\\d+(?=_)")
  # vroom::vroom_write(updated_finemap, glue::glue('{ld_block_data_dir}/{block}/better_finemapped_studies.tsv'))

  updated_bps <- which(updated_finemap$bp != current_finemap$bp)
  unique_studies_to_redo <- updated_finemap$unique_study_id[updated_bps]

  pairwise_results <- vroom::vroom(
    glue::glue("{ld_block_data_dir}/{block}/coloc_pairwise_results.tsv.gz"),
    show_col_types = F
  )
  filtered_pairwise_results <- dplyr::filter(
    pairwise_results,
    !unique_study_a %in% unique_studies_to_redo & !unique_study_b %in% unique_studies_to_redo
  )
  vroom::vroom_write(
    filtered_pairwise_results,
    glue::glue("{ld_block_data_dir}/{block}/coloc_pairwise_results_to_redo.tsv.gz")
  )
  return()
}

recreate_svg_files <- function() {
  source("../pipeline_steps/svg_helpers.R")
  ld_blocks <- vroom::vroom("../pipeline_steps/data/ld_blocks.tsv")

  ld_info <- construct_ld_block(ld_blocks$ancestry, ld_blocks$chr, ld_blocks$start, ld_blocks$stop)
  ld_info <- ld_info[dir.exists(ld_info$ld_block_data), ]

  blocks <- ld_info$block
  # lapply(blocks, function(block) {
  parallel::mclapply(blocks, mc.cores = 30, function(block) {
    print(block)
    finemapped_studies <- vroom::vroom(
      glue::glue("{ld_block_data_dir}/{block}/finemapped_studies.tsv"),
      show_col_types = F
    ) |>
      dplyr::mutate(coverage = dplyr::case_when(
        grepl("godmc", study) ~ coverage_types$sparse,
        TRUE ~ coverage_types$dense
      ))
    vroom::vroom_write(finemapped_studies, glue::glue("{ld_block_data_dir}/{block}/finemapped_studies.tsv"))

    imputed_snps_to_remove <- vroom::vroom(
      glue::glue("{ld_block_data_dir}/{block}/imputation_snps_to_remove.tsv"),
      show_col_types = F
    ) |>
      dplyr::filter(!is.na(bps_to_remove))
    finemapped_snps_to_remove <- vroom::vroom(
      glue::glue("{ld_block_data_dir}/{block}/problematic_finemapped_snps.tsv"),
      show_col_types = F
    )

    finemapped_studies <- dplyr::filter(
      finemapped_studies,
      study %in% imputed_snps_to_remove$study | unique_study_id %in% finemapped_snps_to_remove$unique_study_id
    )
    print(glue::glue("{nrow(finemapped_studies)} svgs to update"))

    lapply(seq_len(nrow(finemapped_studies)), function(i) {
      tryCatch(
        {
          study <- finemapped_studies[i, ]
          is_sparse <- study$coverage == coverage_types$sparse
          gwas <- vroom::vroom(study$file, show_col_types = F)
          create_svg_for_ld_block(gwas, study$study, study$svg_file, block, is_sparse)
        },
        error = function(e) {
          print(study)
          print(glue::glue("Error in {block}: {e}"))
          stop(glue::glue("Error in {block}: {e}"))
        }
      )
      return()
    })
    return()
  })
  return()
}

find_dodgy_ld_blocks <- function() {
  ld_blocks <- vroom::vroom("../pipeline_steps/data/ld_blocks.tsv")
  ld_info <- construct_ld_block(ld_blocks$ancestry, ld_blocks$chr, ld_blocks$start, ld_blocks$stop)
  ld_info <- ld_info[dir.exists(ld_info$ld_block_data), ]
  blocks <- ld_info$block
  lapply(blocks, function(block) {
    imputed_snps_to_remove <- vroom::vroom(
      glue::glue("{ld_block_data_dir}/{block}/imputation_snps_to_remove.tsv"),
      show_col_types = F
    ) |>
      dplyr::filter(!is.na(bps_to_remove))
    imputed_studies <- vroom::vroom(
      glue::glue("{ld_block_data_dir}/{block}/imputed_studies.tsv"),
      show_col_types = F
    ) |>
      dplyr::filter(study %in% imputed_snps_to_remove$study)

    for (i in seq_len(nrow(imputed_studies))) {
      this_study <- imputed_studies[i, ]
      snps_to_remove <- dplyr::filter(imputed_snps_to_remove, study == this_study$study)
      backup_study <- sub("/data/", "/backup/data/", this_study$file)
      backup_study_size <- as.numeric(system(glue::glue("zcat {backup_study} | wc -l"), intern = TRUE))
      study_size <- as.numeric(system(glue::glue("zcat {this_study$file} | wc -l"), intern = TRUE))
      if (backup_study_size - study_size > (length(snps_to_remove$bps_to_remove) * 2)) {
        print(glue::glue("{block}"))
        return()
      }
    }
    return()
  })
  return()
}

copy_over_studies_from_backup <- function(blocks) {
  lapply(blocks, function(block) {
    print(block)
    if (!file.exists(glue::glue("{ld_block_data_dir}/{block}/imputation_snps_to_remove.tsv"))) {
      return()
    }
    if (!file.exists(glue::glue("{ld_block_data_dir}/{block}/problematic_finemapped_snps.tsv"))) {
      return()
    }
    imputed_snps_to_remove <- vroom::vroom(
      glue::glue("{ld_block_data_dir}/{block}/imputation_snps_to_remove.tsv"),
      show_col_types = F
    ) |>
      dplyr::filter(!is.na(bps_to_remove))
    finemapped_snps_to_remove <- vroom::vroom(
      glue::glue("{ld_block_data_dir}/{block}/problematic_finemapped_snps.tsv"),
      show_col_types = F
    )

    imputed_studies <- vroom::vroom(
      glue::glue("{ld_block_data_dir}/{block}/imputed_studies.tsv"),
      show_col_types = F
    ) |>
      dplyr::filter(study %in% unique(imputed_snps_to_remove$study))

    finemapped_studies <- vroom::vroom(
      glue::glue("{ld_block_data_dir}/{block}/finemapped_studies.tsv"),
      show_col_types = F
    ) |>
      dplyr::filter(
        unique_study_id %in% finemapped_snps_to_remove$unique_study_id |
          study %in% imputed_snps_to_remove$study
      )

    message(
      "Copying over ",
      nrow(imputed_studies),
      " imputed studies and ",
      nrow(finemapped_studies),
      " finemapped studies from backup"
    )

    for (i in seq_len(nrow(imputed_studies))) {
      study <- imputed_studies[i, ]
      imputed_study_file <- study$file
      backup_study_imputed_file <- sub("/data/", "/backup/data/", imputed_study_file)
      if (!file.exists(backup_study_imputed_file)) {
        print(glue::glue("Backup file does not exist: {backup_study_imputed_file}"))
        next
      }
      file.copy(backup_study_imputed_file, imputed_study_file, overwrite = TRUE)
    }

    for (i in seq_len(nrow(finemapped_studies))) {
      study <- finemapped_studies[i, ]
      finemapped_study_file <- study$file
      backup_study_finemapped_file <- sub("/data/", "/backup/data/", finemapped_study_file)
      if (!file.exists(backup_study_finemapped_file)) {
        print(glue::glue("Backup file does not exist: {backup_study_finemapped_file}"))
        next
      }
      file.copy(backup_study_finemapped_file, finemapped_study_file, overwrite = TRUE)
    }
    return()
  })
  return()
}

remove_additional_unneeded_colocs_pairs <- function() {
  ld_blocks <- vroom::vroom("../pipeline_steps/data/ld_blocks.tsv")
  ld_info <- construct_ld_block(ld_blocks$ancestry, ld_blocks$chr, ld_blocks$start, ld_blocks$stop)
  ld_info <- ld_info[dir.exists(ld_info$ld_block_data), ]
  blocks <- tail(ld_info$block)
  blocks <- c("EUR/1/153358890-156424168")

  # parallel::mclapply(blocks, mc.cores = 30, function(block) {
  dont_print <- lapply(blocks, function(block) {
    print(block)
    finemapped_studies <- vroom::vroom(
      glue::glue("{ld_block_data_dir}/{block}/finemapped_studies.tsv"),
      show_col_types = F
    )
    finemapped_studies_filtered <- dplyr::filter(finemapped_studies, min_p <= lowest_p_value_threshold)
    finemapped_studies_low_p <- dplyr::filter(finemapped_studies, min_p > lowest_p_value_threshold)
    print(head(dplyr::pull(finemapped_studies_low_p, unique_study_id), n = 30))
    print(
      glue::glue("{nrow(finemapped_studies) - nrow(finemapped_studies_filtered)}",
        " / {nrow(finemapped_studies)} finemapped studies that are weak"
      )
    ) # nolint: line_length_linter.
    coloc_pairwise_results <- vroom::vroom(
      glue::glue("{ld_block_data_dir}/{block}/coloc_pairwise_results.tsv.gz"),
      show_col_types = F
    )
    coloc_pairwise_results_filtered <- dplyr::filter(
      coloc_pairwise_results,
      unique_study_a %in% finemapped_studies_filtered$unique_study_id &
        unique_study_b %in% finemapped_studies_filtered$unique_study_id
    )
    if (nrow(coloc_pairwise_results_filtered) == nrow(coloc_pairwise_results)) {
      return()
    }
    print(glue::glue("Removing {nrow(coloc_pairwise_results) - nrow(coloc_pairwise_results_filtered)} coloc pairwise results out of {nrow(coloc_pairwise_results)}")) # nolint: line_length_linter.
    # vroom::vroom_write(coloc_pairwise_results, glue::glue('{ld_block_data_dir}/{block}/coloc_pairwise_results.tsv.gz')) # nolint: line_length_linter.
    return()
  })
  return()
}

remove_bad_studies_fix <- function() {
  source("../pipeline_steps/gwas_calculations.R")
  ld_blocks <- vroom::vroom("../pipeline_steps/data/ld_blocks.tsv")
  ld_info <- construct_ld_block(ld_blocks$ancestry, ld_blocks$chr, ld_blocks$start, ld_blocks$stop)
  ld_info <- ld_info[dir.exists(ld_info$ld_block_data), ]
  blocks <- ld_info$block
  blocks <- c("EUR/3/172148025-174557918", "EUR/5/7023201-9205437", "EUR/2/227191954-229419238")

  dont_print <- parallel::mclapply(blocks, mc.cores = 30, function(block) {
    # dont_print <- lapply(blocks, function(block) {
    print(block)
    if (!file.exists(glue::glue("{ld_block_data_dir}/{block}/imputation_snps_to_remove.tsv"))) {
      return()
    }
    imputed_snps_to_remove <- vroom::vroom(
      glue::glue("{ld_block_data_dir}/{block}/imputation_snps_to_remove.tsv"),
      show_col_types = F
    ) |>
      dplyr::filter(!is.na(bps_to_remove))
    finemapped_snps_to_remove <- vroom::vroom(
      glue::glue("{ld_block_data_dir}/{block}/problematic_finemapped_snps.tsv"),
      show_col_types = F
    )

    really_bad_studies <- table(imputed_snps_to_remove$study)
    really_bad_studies <- names(really_bad_studies[really_bad_studies > 150])

    imputed_studies <- vroom::vroom(
      glue::glue("{ld_block_data_dir}/{block}/imputed_studies.tsv"),
      show_col_types = F
    ) |>
      dplyr::filter(study %in% really_bad_studies)

    print(glue::glue("Updating {nrow(imputed_studies)} imputed studies"))

    dont_print <- lapply(seq_len(nrow(imputed_studies)), function(i) {
      study <- imputed_studies[i, ]
      imputed_gwas <- vroom::vroom(study$file, show_col_types = F)
      if (!"IMPUTED" %in% colnames(imputed_gwas)) {
        print(glue::glue("Updating {study$file}"))
        imputed_gwas$IMPUTED <- FALSE
        vroom::vroom_write(imputed_gwas, study$file)
      }
      return()
    })
    return()
  })
  return()
}

cleanup_bad_snps <- function(block) {
  if (!file.exists(glue::glue("{ld_block_data_dir}/{block}/imputation_snps_to_remove_updated.tsv"))) {
    more_imputed_snps_to_remove <- data.frame()
  }
  if (!file.exists(glue::glue("{ld_block_data_dir}/{block}/imputation_snps_to_remove.tsv"))) {
    imputed_snps_to_remove <- data.frame()
  }

  imputed_snps_to_remove <- vroom::vroom(
    glue::glue("{ld_block_data_dir}/{block}/imputation_snps_to_remove.tsv"),
    show_col_types = F
  ) |>
    dplyr::filter(!is.na(bps_to_remove))
  more_imputed_snps_to_remove <- vroom::vroom(
    glue::glue("{ld_block_data_dir}/{block}/imputation_snps_to_remove_updated.tsv"),
    show_col_types = F
  ) |>
    dplyr::filter(!is.na(bps_to_remove))

  imputed_snps_to_remove <- dplyr::bind_rows(imputed_snps_to_remove, more_imputed_snps_to_remove)

  unique_altered_studies <- unique(imputed_snps_to_remove$study)

  imputed_studies <- vroom::vroom(glue::glue("{ld_block_data_dir}/{block}/imputed_studies.tsv"), show_col_types = F) |>
    dplyr::filter(study %in% imputed_snps_to_remove$study)

  print(glue::glue("{block}: Updating {nrow(imputed_studies)} imputed studies"))

  altered_studies <- lapply(seq_len(nrow(imputed_studies)), function(i) {
    study <- imputed_studies[i, ]
    study_name <- study$study
    specific_to_remove <- dplyr::filter(imputed_snps_to_remove, study == study_name & !is.na(bps_to_remove))
    if (nrow(specific_to_remove) == 0) {
      return(NULL)
    }

    print(paste(block, study$file))
    gwas <- vroom::vroom(study$file, show_col_types = F)
    updated_gwas <- dplyr::filter(gwas, !BP %in% specific_to_remove$bps_to_remove)
    rows_to_remove <- nrow(gwas) - nrow(updated_gwas)
    if (rows_to_remove == 0) {
      return(NULL)
    }
    # print(glue::glue('{block}: Removed {rows_to_remove} SNPs in {study$study}'))

    vroom::vroom_write(updated_gwas, study$file)
    return(study_name)
  })
  altered_studies <- unlist(altered_studies)

  print(glue::glue("{block}: Altered {length(altered_studies)} studies out of {length(unique(imputed_snps_to_remove$study))}")) # nolint: line_length_linter.
  if (length(altered_studies) == 0) {
    return()
  }

  finemapped_studies <- vroom::vroom(
    glue::glue("{ld_block_data_dir}/{block}/finemapped_studies.tsv"),
    show_col_types = F
  )

  finemapped_studies_to_alter <- finemapped_studies |>
    dplyr::filter(study %in% imputed_snps_to_remove$study)

  print(glue::glue("{block}: Updating {nrow(finemapped_studies_to_alter)} finemapped studies"))

  unchanged_finemapped_studies <- finemapped_studies |>
    dplyr::filter(!study %in% imputed_snps_to_remove$study)

  updated_finemapped_studies <- lapply(seq_len(nrow(finemapped_studies_to_alter)), function(i) {
    study <- finemapped_studies_to_alter[i, ]
    study_name <- study$study
    gwas <- vroom::vroom(study$file, show_col_types = F)
    specific_to_remove <- dplyr::filter(imputed_snps_to_remove, study == study_name & !is.na(bps_to_remove))

    if (nrow(specific_to_remove) == 0) {
      return(study)
    }

    if (nrow(gwas) == 0) {
      print(glue::glue("No rows in file: {study$file}, remove from finemapped studies"))
      return(NULL)
    }

    updated_gwas <- dplyr::filter(gwas, !BP %in% specific_to_remove$bps_to_remove)
    rows_to_remove <- nrow(gwas) - nrow(updated_gwas)

    if (rows_to_remove != 0) {
      # print(glue::glue('{block}: Removed {rows_to_remove} SNPs in {study$unique_study_id}'))
      vroom::vroom_write(updated_gwas, study$file)
    }
    if (nrow(updated_gwas) == 0) {
      print(glue::glue("No rows in file: {study$file}, remove from finemapped studies"))
      return(NULL)
    }

    # if the snp is not in the updated gwas, find the new snp most significant SNP
    if (!study$snp %in% updated_gwas$SNP) {
      updated_gwas$NEW_Z <- convert_lbf_to_abs_z(updated_gwas$LBF, updated_gwas$SE)
      interesting_snp <- updated_gwas[which.max(updated_gwas$NEW_Z), ]
      study$bp <- interesting_snp$BP
      study$snp <- interesting_snp$SNP
      study$min_p <- 2 * pnorm(-abs(interesting_snp$NEW_Z))
      # print(glue::glue('{block}: Found new SNP {study$snp} in updated gwas for {study$unique_study_id}'))
    }

    return(study)
  }) |> dplyr::bind_rows()

  updated_finemapped_studies <- dplyr::bind_rows(updated_finemapped_studies, unchanged_finemapped_studies)
  print(glue::glue("{block}: Updated finemapped studies: {nrow(updated_finemapped_studies)}, original: {nrow(finemapped_studies)}")) # nolint: line_length_linter.
  if (nrow(updated_finemapped_studies) != nrow(finemapped_studies)) {
    stop(glue::glue("ERROR: Updated finemapped studies: {nrow(updated_finemapped_studies)}, original: {nrow(finemapped_studies)}")) # nolint: line_length_linter.
  }
  vroom::vroom_write(updated_finemapped_studies, glue::glue("{ld_block_data_dir}/{block}/finemapped_studies.tsv"))

  non_significant_finemapped_studies <- updated_finemapped_studies |>
    dplyr::filter(min_p > lowest_p_value_threshold)
  print(glue::glue("{block}: Too weak finemapped studies: {nrow(non_significant_finemapped_studies)}"))

  coloc_pairwise_results <- vroom::vroom(
    glue::glue("{ld_block_data_dir}/{block}/coloc_pairwise_results.tsv.gz"),
    show_col_types = F
  )
  new_coloc_pairwise_results <- coloc_pairwise_results |>
    dplyr::filter(!(
      (PP.H4.abf > 0.5 | PP.H3.abf > 0.5) &
        # (study_a %in% imputed_snps_to_remove$study | study_b %in% imputed_snps_to_remove$study |
        (study_a %in% altered_studies | study_b %in% altered_studies |
           unique_study_a %in% non_significant_finemapped_studies$unique_study_id |
           unique_study_b %in% non_significant_finemapped_studies$unique_study_id)
    ))
  print(glue::glue("{block}: Removing {nrow(coloc_pairwise_results) - nrow(new_coloc_pairwise_results)} coloc pairwise results out of {nrow(coloc_pairwise_results)}")) # nolint: line_length_linter.
  vroom::vroom_write(
    new_coloc_pairwise_results,
    glue::glue("{ld_block_data_dir}/{block}/coloc_pairwise_results.tsv.gz")
  )
  return()
}

cleanup_all_problematic_snps <- function(blocks) {
  source("../pipeline_steps/gwas_calculations.R")
  ld_blocks <- vroom::vroom("../pipeline_steps/data/ld_blocks.tsv")
  ld_info <- construct_ld_block(ld_blocks$ancestry, ld_blocks$chr, ld_blocks$start, ld_blocks$stop)
  ld_info <- ld_info[dir.exists(ld_info$ld_block_data), ]
  blocks <- ld_info$block

  blocks_to_skip <- vroom::vroom("finished_blocks.csv", show_col_types = F, delim = ",")$block
  blocks <- setdiff(blocks, blocks_to_skip)
  print(length(blocks))
  print(blocks)

  dont_print <- parallel::mclapply(blocks, mc.cores = 30, function(block) {
    # dont_print <- lapply(blocks, function(block) {
    print(block)
    tryCatch(
      {
        cleanup_bad_snps(block)
      },
      error = function(e) {
        print(glue::glue("Error in {block}: {e}"))
        stop(glue::glue("Error in {block}: {e}"))
      }
    )
    return()
  })
  return()
}

find_problematic_finemapped_results_for_ld_block <- function(ld_block) {
  finemapped_studies <- vroom::vroom(
    glue::glue("{ld_block_data_dir}/{ld_block}/finemapped_studies.tsv"),
    show_col_types = F
  )
  imputed_studies <- vroom::vroom(glue::glue("{ld_block_data_dir}/{ld_block}/imputed_studies.tsv"), show_col_types = F)

  all_problematic_finemapped_snps <- lapply(seq_len(nrow(finemapped_studies)), function(i) {
    finemapped_study <- finemapped_studies[i, ]
    imputed_study <- dplyr::filter(imputed_studies, study == finemapped_study$study)
    finemapped_gwas <- vroom::vroom(finemapped_study$file, show_col_types = F)

    if (nrow(finemapped_gwas) == 0) {
      print(glue::glue("ERROR: No rows in file: {finemapped_study$file}"))
      return(data.frame(unique_study_id = finemapped_study$unique_study_id, error = "No rows in file"))
    }

    finemapped_gwas <- dplyr::mutate(finemapped_gwas, P = 2 * pnorm(-abs(convert_lbf_to_abs_z(LBF, SE))))
    significant_finemapped_rows <- dplyr::filter(finemapped_gwas, P < lowest_p_value_threshold)
    if (nrow(significant_finemapped_rows) == 0) {
      return()
    }

    imputed_gwas <- vroom::vroom(imputed_study$file, show_col_types = F)

    if (nrow(imputed_gwas) == 0 || is.null(imputed_gwas$SNP) || is.na(imputed_gwas$SNP)) {
      print(glue::glue("{imputed_study$study} {finemapped_study$unique_study_id} has no SNPs?!?"))
      return()
    }

    problematic_finemapped_snps <- imputed_gwas |>
      dplyr::filter(SNP %in% significant_finemapped_rows$SNP & P > 1e-3) |>
      dplyr::select(SNP, BETA, SE, P) |>
      dplyr::mutate(unique_study_id = finemapped_study$unique_study_id) |>
      dplyr::left_join(dplyr::select(significant_finemapped_rows, SNP, LBF), by = "SNP")
    if (nrow(problematic_finemapped_snps) == 0) {
      return()
    }

    return(problematic_finemapped_snps)
  }) |>
    dplyr::bind_rows()

  if (nrow(all_problematic_finemapped_snps) == 0) {
    return()
  }
  print(glue::glue("{nrow(all_problematic_finemapped_snps)} SNPs over {length(unique(all_problematic_finemapped_snps$unique_study_id))} studies for {ld_block}")) # nolint: line_length_linter.
  vroom::vroom_write(
    all_problematic_finemapped_snps,
    glue::glue("{ld_block_data_dir}/{ld_block}/problematic_finemapped_snps.tsv")
  )
  return()
}

find_bad_imputations <- function() {
  source("../pipeline_steps/imputation_method.R")
  ld_blocks <- vroom::vroom("../pipeline_steps/data/ld_blocks.tsv")
  ld_info <- construct_ld_block(ld_blocks$ancestry, ld_blocks$chr, ld_blocks$start, ld_blocks$stop)
  ld_info <- ld_info[dir.exists(ld_info$ld_block_data), ]
  blocks <- ld_info$block
  # blocks <- blocks[60:length(blocks)]
  parallel::mclapply(blocks, mc.cores = 35, function(block) {
    print(block)
    investigate_bad_imputation(block, glue::glue("{ld_block_data_dir}/{block}/imputed_studies.tsv"))
    return()
  })
  return()
}

investigate_bad_imputation <- function(ld_block, imputed_file) {
  imputed_studies <- vroom::vroom(imputed_file, show_col_types = F)
  ld_matrix <- vroom::vroom(
    glue::glue("{ld_reference_panel_dir}/{ld_block}.unphased.vcor1.gz"),
    col_names = F,
    show_col_types = F
  )
  to_remove_file <- glue::glue("{ld_block_data_dir}/{ld_block}/imputation_snps_to_remove.tsv")

  if (!file.exists(to_remove_file)) {
    to_remove <- data.frame(study = character(), bps_to_remove = numeric())
  } else {
    to_remove <- vroom::vroom(to_remove_file, show_col_types = F)
  }

  updated_to_remove <- lapply(seq_len(nrow(imputed_studies)), function(i) {
    study <- imputed_studies[i, ]
    gwas <- vroom::vroom(study$file, show_col_types = F)

    result <- filter_imputation_results(gwas, ld_matrix, min(gwas$BP), max(gwas$BP))

    if (!is.null(result$bps_to_remove) && length(result$bps_to_remove) > 0) {
      candidate_rows <- data.frame(study = study$study, bps_to_remove = result$bps_to_remove, stringsAsFactors = FALSE)
      if (nrow(to_remove) > 0) {
        return(dplyr::anti_join(candidate_rows, to_remove, by = c("study", "bps_to_remove")))
      } else {
        return(candidate_rows)
      }
    }
  })

  updated_to_remove_file <- glue::glue("{ld_block_data_dir}/{ld_block}/imputation_snps_to_remove_updated.tsv")
  print(glue::glue("{ld_block}: {length(updated_to_remove)} SNPs to remove"))
  vroom::vroom_write(dplyr::bind_rows(updated_to_remove), updated_to_remove_file)
  return()
}

repopulate_missing_finemapped_results <- function() {
  ld_blocks <- vroom::vroom("../pipeline_steps/data/ld_blocks.tsv")
  ld_info <- construct_ld_block(ld_blocks$ancestry, ld_blocks$chr, ld_blocks$start, ld_blocks$stop)
  ld_info <- ld_info[dir.exists(ld_info$ld_block_data), ]
  ld_info <- ld_info[15:nrow(ld_info), , drop = F]
  return()
}

update_trait_names <- function() {
  traits_processed <- vroom::vroom("/local-scratch/projects/genotype-phenotype-map/results/latest/traits_processed.tsv.gz") # nolint: line_length_linter.
  studies_processed <- vroom::vroom("/local-scratch/projects/genotype-phenotype-map/results/latest/studies_processed.tsv.gz") # nolint: line_length_linter.
  traits_processed <- dplyr::left_join(
    traits_processed,
    studies_processed[, c("study_name", "gene")],
    by = "study_name"
  )

  new_methqtl_trait_names <- paste(
    traits_processed$gene[grepl("methQTL", traits_processed$trait) & !is.na(traits_processed$gene)],
    traits_processed$trait[grepl("methQTL", traits_processed$trait) & !is.na(traits_processed$gene)]
  )
  traits_processed$trait[grepl(
    "methQTL",
    traits_processed$trait
  ) & !is.na(traits_processed$gene)] <- new_methqtl_trait_names

  sqtl_study_names <- traits_processed[grepl("sQTL", traits_processed$trait), ]$study_name
  sqtl_study_names <- sub("^.*[- ](chr.*):clu.*$", "\\1", sqtl_study_names)
  new_names <- paste(traits_processed$trait[grepl("sQTL", traits_processed$trait)], sqtl_study_names)
  traits_processed$trait[grepl("sQTL", traits_processed$trait)] <- new_names

  vroom::vroom_write(
    traits_processed,
    "/local-scratch/projects/genotype-phenotype-map/results/latest/traits_processed_new.tsv.gz"
  )
  return()
}

populate_godmc_gene_info <- function() {
  besd_epi <- vroom::vroom("/local-scratch/data/hg38/godmc/methylation.epi", col_names = F)
  sp <- vroom::vroom(glue::glue("{latest_results_dir}/studies_processed.tsv.gz"))
  epic <- vroom::vroom(glue::glue("{variant_annotation_dir}/EPICv2.hg38.manifest.gencode.v41.tsv.gz"))
  msa <- vroom::vroom(glue::glue("{variant_annotation_dir}/MSA.hg38.manifest.gencode.v41.tsv.gz")) |>
    dplyr::filter(!probeID %in% epic$probeID)
  gene_info <- vroom::vroom(glue::glue("{variant_annotation_dir}/gene_info.tsv"), show_col_types = F)

  all_methylation <- dplyr::bind_rows(epic, msa) |>
    dplyr::filter(!is.na(genesUniq))

  all_methylation$probe <- sub("_.*", "", all_methylation$probeID)
  all_methylation$gene <- sub(";.*", "", all_methylation$genesUniq)

  all_methylation$ensg <- ifelse(
    !grepl("ENSG", all_methylation$gene),
    gene_info$ensembl_id[match(all_methylation$gene, gene_info$gene)],
    all_methylation$gene
  )

  all_methylation$ensg <- ifelse(
    !grepl("ENSG", all_methylation$ensg),
    NA,
    all_methylation$ensg
  )

  all_methylation <- all_methylation |> dplyr::filter(!is.na(ensg))


  epic_for_epi <- all_methylation[, c("probe", "ensg")] |>
    dplyr::rename(X2 = probe, X5 = ensg) |>
    dplyr::distinct(X2, .keep_all = TRUE)

  epic_for_sp <- all_methylation[, c("probe", "gene", "ensg")] |>
    dplyr::distinct(probe, .keep_all = TRUE)

  sp <- sp |>
    dplyr::rows_update(epic_for_sp, by = "probe", unmatched = "ignore")

  besd_epi <- besd_epi |>
    dplyr::rows_update(epic_for_epi, by = "X2", unmatched = "ignore")
  besd_epi$X5[!grepl("ENSG", besd_epi$X5)] <- NA

  vroom::vroom_write(sp, glue::glue("{latest_results_dir}/studies_processed.tsv.gz"))
  vroom::vroom_write(besd_epi, "/local-scratch/data/hg38/godmc/methylation.epi", col_names = F)
  return()
}

fix_all_svg_extractions <- function() {
  source("../pipeline_steps/svg_helpers.R")
  studies <- vroom::vroom(glue::glue("{latest_results_dir}/studies_processed.tsv.gz"))
  study_extractions <- vroom::vroom(glue::glue("{latest_results_dir}/study_extractions.tsv.gz")) |>
    dplyr::left_join(studies, by = c("study" = "study_name"))

  study_extractions <- study_extractions[630000:nrow(study_extractions), ]
  # lapply(seq_len(nrow(study_extractions)), function(i) {
  parallel::mclapply(seq_len(nrow(study_extractions)), mc.cores = 50, function(i) {
    study <- study_extractions[i, ]
    print(i)
    gwas <- vroom::vroom(study$file, show_col_types = F)
    is_sparse <- study$data_type == data_types$methylation
    create_svg_for_ld_block(gwas, study$study, study$svg_file, study$ld_block, is_sparse = is_sparse)
    return()
  })
  return()
}

fix_bad_full_svgs <- function() {
  source("../pipeline_steps/svg_helpers.R")
  bad_svg_studies <- vroom::vroom(glue::glue("{latest_results_dir}/studies_processed.tsv.gz")) |>
    dplyr::filter(data_type == "phenotype" & variant_type == "common")

  # Get all bad_svg_studies rows after the row where study == 'ebi-a-GCST90025995'
  idx <- which(bad_svg_studies$study == "ebi-a-GCST90025995")
  if (length(idx) > 0 && idx < nrow(bad_svg_studies)) {
    bad_svg_studies <- bad_svg_studies[(idx + 1):nrow(bad_svg_studies), ]
  } else if (length(idx) > 0 && idx == nrow(bad_svg_studies)) {
    bad_svg_studies <- bad_svg_studies[0, ]
  }

  lapply(seq_len(nrow(bad_svg_studies)), function(i) {
    study <- bad_svg_studies[i, ]
    vcf_file <- glue::glue("{study$extracted_location}/vcf/hg38.vcf.gz")
    bcf_query <- glue::glue(
      "/home/bcftools/bcftools query ",
      '--format "[%ID]\t[%CHROM]\t[%POS]\t[%LP]" ',
      "{vcf_file}"
    )
    gwas <- system(bcf_query, wait = T, intern = T)
    gwas <- data.table::fread(text = gwas)

    colnames(gwas) <- c("RSID", "CHR", "BP", "LP")
    gwas <- gwas |> dplyr::mutate(LP = as.numeric(LP))
    max_lp <- max(gwas$LP, na.rm = TRUE)
    if (!is.na(max_lp) && (is.infinite(max_lp) || max_lp >= 300)) {
      print(study$study_name)
      print(max_lp)
      create_svgs_from_gwas(study, gwas)
    }
    return()
  })
  return()
}

fix_missing_snps_in_rare_results <- function() {
  source("../pipeline_steps/common_extraction_functions.R")
  ld_blocks <- vroom::vroom("../pipeline_steps/data/ld_blocks.tsv")
  ld_info <- construct_ld_block(ld_blocks$ancestry, ld_blocks$chr, ld_blocks$start, ld_blocks$stop)
  ld_info <- ld_info[dir.exists(ld_info$ld_block_data), ]

  missing_snps <- vroom::vroom(glue::glue("{current_results_dir}/missing_snps_in_study_extractions.tsv")) |>
    dplyr::filter(grepl("-wes-", unique_study_id) & grepl("TRUE", snp))

  message(glue::glue("Fixing {nrow(missing_snps)} missing SNPs in study extractions"))
  dont_print <- lapply(seq_len(nrow(missing_snps)), function(i) {
    missing_snp <- missing_snps[i, ]
    message(glue::glue("Fixing {missing_snp$unique_study_id} {missing_snp$snp} {missing_snp$ld_block}"))
    gwas <- vroom::vroom(missing_snp$file, show_col_types = F) |>
      dplyr::mutate(
        EA = dplyr::case_when(
          grepl("TRUE", EA) ~ "T",
          T ~ "A"
        ),
        OA = dplyr::case_when(
          grepl("TRUE", OA) ~ "T",
          T ~ "A"
        ),
        SNP = dplyr::case_when(
          grepl("TRUE", SNP) ~ glue::glue("{CHR}:{BP}_{EA}_{OA}"),
          T ~ SNP
        )
      )

    # vroom::vroom_write(gwas, missing_snp$file)

    new_snp <- sub("TRUE", "T", missing_snp$snp)

    message(glue::glue("Updated {missing_snp$unique_study_id} to {new_snp}"))

    rare_results_file <- glue::glue("{ld_block_data_dir}/{missing_snp$ld_block}/compare_rare_results.tsv")
    rare_results <- vroom::vroom(rare_results_file, show_col_types = F)
    rare_results$candidate_snp[rare_results$candidate_snp == missing_snp$snp] <- new_snp
    message(glue::glue("Rows updated: {sum(rare_results$candidate_snp == missing_snp$snp)}"))

    vroom::vroom_write(rare_results, rare_results_file)
    return()
  })
  return()
}

fix_missing_snps_in_study_extractions <- function() {
  source("../pipeline_steps/common_extraction_functions.R")
  ld_blocks <- vroom::vroom("../pipeline_steps/data/ld_blocks.tsv")
  ld_info <- construct_ld_block(ld_blocks$ancestry, ld_blocks$chr, ld_blocks$start, ld_blocks$stop)
  ld_info <- ld_info[dir.exists(ld_info$ld_block_data), ]

  missing_snps <- vroom::vroom(glue::glue("{current_results_dir}/missing_snps_in_study_extractions.tsv")) |>
    dplyr::filter(!grepl("-wes-", unique_study_id))

  message(glue::glue("Fixing {nrow(missing_snps)} missing SNPs in study extractions"))
  lapply(seq_len(nrow(missing_snps)), function(i) {
    missing_snp <- missing_snps[i, ]
    message(glue::glue("Fixing {missing_snp$unique_study_id} {missing_snp$snp}"))
    message(missing_snp$file_with_lbfs)
    gwas <- vroom::vroom(file.path(data_dir, missing_snp$file_with_lbfs), show_col_types = F) |>
      dplyr::filter(!is.na(CHR) & !is.na(BP) & !is.na(EA) & !is.na(OA) & !is.na(EAF)) |>
      dplyr::mutate(COMPRESSED_SNP = format_compressed_allele_snp_string(CHR, BP, EA, OA))

    new_snp <- gwas |>
      dplyr::filter(COMPRESSED_SNP == missing_snp$snp) |>
      dplyr::pull(SNP)

    if (length(new_snp) == 0) {
      message(glue::glue("No SNP found for {missing_snp$unique_study_id} {missing_snp$snp}, using max abs_beta_over_se")) # nolint: line_length_linter.
      new_snp <- gwas |>
        dplyr::filter(Z == max(Z)) |>
        dplyr::pull(SNP)
    }

    message(glue::glue("Updated {missing_snp$unique_study_id} to {new_snp}"))

    finemapped_studies_file <- glue::glue("{ld_block_data_dir}/{missing_snp$ld_block}/finemapped_studies.tsv")
    finemapped_studies <- vroom::vroom(finemapped_studies_file, show_col_types = F)
    finemapped_studies$snp[finemapped_studies$unique_study_id == missing_snp$unique_study_id] <- new_snp
    message(glue::glue("Rows updated: {sum(finemapped_studies$unique_study_id == missing_snp$unique_study_id)}"))

    vroom::vroom_write(finemapped_studies, finemapped_studies_file)
    return()
  })
  return()
}

delete_some_coloc_complete <- function() {
  ld_blocks <- vroom::vroom("../pipeline_steps/data/ld_blocks.tsv")
  ld_info <- construct_ld_block(ld_blocks$ancestry, ld_blocks$chr, ld_blocks$start, ld_blocks$stop)
  ld_info <- ld_info[dir.exists(ld_info$ld_block_data), ]

  target_date_str <- "2025-07-19" # The cutoff date (YYYY-MM-DD)
  file_to_delete <- "coloc_complete" # The specific file name to look for
  file_to_check <- "coloc_clustered_results.tsv.gz" # The specific file name to look for
  target_date <- as.Date(target_date_str)
  all_dirs <- ld_info$ld_block_data

  # Loop through each directory
  for (dir_path in all_dirs) {
    # print(dir_path)
    file_path <- file.path(dir_path, file_to_delete)
    file_path_check <- file.path(dir_path, file_to_check)

    # Check if the file exists in the current directory
    if (file.exists(file_path_check)) {
      # Get file information, specifically the modification time
      file_info <- file.info(file_path_check)
      file_mtime <- as.Date(file_info$mtime) # Convert modification time to Date
      message(file_mtime, " ", target_date)

      if (file_mtime < target_date) {
        cat(paste0("  -> '", file_to_check, "' is OLDER than ", target_date_str, ". Deleting...\n", file_path))
        file.remove(file_path)
      }
    } else {
      cat(paste0("  -> '", file_to_check, "' is missing, delete."))
      file.remove(file_path)
    }
  }
  return()
}

remove_duplicate_unique_study_ids <- function() {
  # duplicates <- vroom::vroom(glue::glue('{current_results_dir}/study_extractions.tsv.gz'))
  # duplicates <- duplicates[duplicated(duplicates$unique_study_id), ]

  ld_blocks <- vroom::vroom("../pipeline_steps/data/ld_blocks.tsv")
  ld_info <- construct_ld_block(ld_blocks$ancestry, ld_blocks$chr, ld_blocks$start, ld_blocks$stop)
  ld_info <- ld_info[dir.exists(ld_info$ld_block_data), ]
  # ld_info <- ld_info[ld_info$block == 'EUR/1/24201011-28100360', ]

  dont_print <- lapply(ld_info$ld_block_data, function(ld_block) {
    finemapped_studies_file <- glue::glue("{ld_block}/finemapped_studies.tsv")
    finemapped_studies <- vroom::vroom(finemapped_studies_file, show_col_types = F)
    dedup_finemapped_studies <- finemapped_studies[!duplicated(finemapped_studies$unique_study_id), ]
    print(paste("deduped finemapped studies:", nrow(finemapped_studies), "to", nrow(dedup_finemapped_studies)))
    vroom::vroom_write(dedup_finemapped_studies, finemapped_studies_file)
    return()
  })
  return()
}

fix_rare_unique_study_ids <- function() {
  ld_blocks <- vroom::vroom("../pipeline_steps/data/ld_blocks.tsv")
  ld_info <- construct_ld_block(ld_blocks$ancestry, ld_blocks$chr, ld_blocks$start, ld_blocks$stop)
  ld_info <- ld_info[dir.exists(ld_info$ld_block_data), ]

  dont_print <- lapply(ld_info$ld_block_data, function(ld_block) {
    print(ld_block)
    rare_results_file <- glue::glue("{ld_block}/compare_rare_results.tsv")
    if (!file.exists(rare_results_file)) {
      print(paste("file missing:", rare_results_file))
      return()
    }
    rare_results <- vroom::vroom(rare_results_file, show_col_types = F)

    rare_results <- rare_results |>
      tidyr::separate_rows(traits, min_ps, genes, files, sep = ", ") |>
      dplyr::mutate(traits = mapply(function(trait, snp) {
        new_bp <- strsplit(snp, ":")[[1]][2]
        new_bp <- gsub("_", "-", new_bp)
        return(sub("(.*)_[0-9A-Za-z-]+", paste0("\\1_", new_bp), trait))
      }, traits, candidate_snp)) |>
      dplyr::group_by(candidate_snp, ld_block) |>
      dplyr::summarise(
        traits = paste(traits, collapse = ", "),
        min_ps = paste(min_ps, collapse = ", "),
        genes = paste(genes, collapse = ", "),
        files = paste(files, collapse = ", "),
        .groups = "drop"
      )
    vroom::vroom_write(
      rare_results |> dplyr::select(traits, candidate_snp, min_ps, genes, ld_block, files),
      rare_results_file
    )
    return()
  })
  return()
}

remove_na_coloc_results <- function() {
  ld_blocks <- vroom::vroom("../pipeline_steps/data/ld_blocks.tsv")
  ld_info <- construct_ld_block(ld_blocks$ancestry, ld_blocks$chr, ld_blocks$start, ld_blocks$stop)
  lapply(ld_info$ld_block_data, function(ld_block) {
    coloc_results_file <- glue::glue("{ld_block}/coloc_pairwise_results.tsv.gz")
    if (!file.exists(coloc_results_file)) {
      return()
    }
    coloc_results <- vroom::vroom(coloc_results_file, show_col_types = F)
    orig_size <- nrow(coloc_results)
    coloc_results <- dplyr::filter(coloc_results, !is.na(h4))
    new_size <- nrow(coloc_results)
    message(glue::glue("{ld_block}: removing {orig_size - new_size} rows: {(orig_size - new_size)/orig_size*100}%"))

    if (orig_size != new_size) {
      vroom::vroom_write(coloc_results, coloc_results_file)
    }
    return()
  })
  return()
}

fix_existing_min_p_values <- function() {
  ld_blocks <- vroom::vroom("../pipeline_steps/data/ld_blocks.tsv")
  ld_info <- construct_ld_block(ld_blocks$ancestry, ld_blocks$chr, ld_blocks$start, ld_blocks$stop)
  ld_info <- ld_info[dir.exists(ld_info$ld_block_data), ]
  ld_info <- ld_info[15:nrow(ld_info), , drop = F]

  silent <- lapply(seq_len(nrow(ld_info)), function(i) {
    gc()
    block <- ld_info[i, , drop = F]
    print(block$block)
    finemapped <- vroom::vroom(paste0(block$ld_block_data, "/finemapped_studies.tsv"), show_col_types = F)
    imputed <- vroom::vroom(paste0(block$ld_block_data, "/imputed_studies.tsv"), show_col_types = F)

    updated_p_vals <- lapply(seq_len(nrow(finemapped)), function(i) {
      finemap_study <- finemapped[i, , drop = F]
      finemap_p <- as.numeric(finemap_study[["min_p"]])

      more_than_one_finemap_file <- nrow(dplyr::filter(finemapped, study == finemap_study[["study"]])) > 1
      if (!more_than_one_finemap_file) {
        # message(paste('skipping', finemap_study[["unique_study_id"]]))
        return(finemap_study)
      }

      imputed_file <- dplyr::filter(imputed, study == finemap_study[["study"]])$file
      # message("Checking ", finemap_study[["unique_study_id"]], " ", imputed_file)
      if (length(imputed_file) == 0 || is.null(imputed_file) || is.na(imputed_file)) {
        return(finemap_study)
      }

      imputed_gwas <- vroom::vroom(imputed_file, show_col_types = F) |> dplyr::slice(which.min(P))
      imputed_p <- as.numeric(imputed_gwas$P)
      if (finemap_p < imputed_p) {
        message(paste(
          "updating study:",
          finemap_study[["unique_study_id"]],
          "from finemap:",
          finemap_p,
          "to imputed:",
          imputed_p
        ))
        finemap_study$min_p <- imputed_p
      }

      return(finemap_study)
    })

    finemapped <- do.call(rbind, updated_p_vals)
    message(nrow(finemapped), " finemapped studies after update")
    vroom::vroom_write(finemapped, paste0(block$ld_block_data, "/finemapped_studies.tsv"))
    return()
  })
  return()
}

add_missing_svgs_to_finemapped_studies <- function() {
  source("../pipeline_steps/svg_helpers.R")

  ld_blocks <- vroom::vroom("../pipeline_steps/data/ld_blocks.tsv")
  ld_info <- construct_ld_block(ld_blocks$ancestry, ld_blocks$chr, ld_blocks$start, ld_blocks$stop)
  ld_info <- ld_info[dir.exists(ld_info$ld_block_data), ]
  ld_info <- dplyr::filter(ld_info, block == "EUR/22/35695831-37164130")

  results <- lapply(seq_len(nrow(ld_info)), function(i) {
    block <- ld_info[i, ]
    print(block$block)
    finemapped_studies_file <- paste0(block$ld_block_data, "/finemapped_studies.tsv")
    finemapped_studies <- vroom::vroom(finemapped_studies_file, show_col_types = F)

    updated_studies <- lapply(seq_len(nrow(finemapped_studies)), function(j) {
      study <- finemapped_studies[j, ]
      if (!is.na(study$svg_file)) {
        return(study)
      }
      message("creating svg for", study$study)

      svg_file <- sub("finemapped", "svgs/extractions", study$file)
      svg_file <- sub("tsv.gz", "svg", svg_file)
      study$svg_file <- svg_file
      gwas <- vroom::vroom(study$file, show_col_types = F)
      create_svg_for_ld_block(gwas, study$study, svg_file)
      return(study)
    })

    updated_studies <- do.call(rbind, updated_studies)
    if (nrow(updated_studies) == nrow(finemapped_studies)) {
      vroom::vroom_write(updated_studies, finemapped_studies_file)
    } else {
      print(paste("ERROR: number of studies changed for", block$block))
    }
    return()
  })
  return()
}

big_update_to_finemapped_results <- function() {
  # Control data.table threading to prevent spawning new processes
  data.table::setDTthreads(1)

  # Control ggplot2 to prevent excessive process spawning
  options(ggplot2.parallel = FALSE)
  Sys.setenv(OMP_THREAD_LIMIT = 1)

  source("../pipeline_steps/svg_helpers.R")

  ld_blocks <- data.table::fread("../pipeline_steps/data/ld_blocks.tsv")
  ld_info <- construct_ld_block(ld_blocks$ancestry, ld_blocks$chr, ld_blocks$start, ld_blocks$stop)
  ld_info <- ld_info[dir.exists(ld_info$ld_block_data), ]

  ld_info <- ld_info[4:nrow(ld_info), ]
  ld_info <- dplyr::filter(ld_info, block == "EUR/22/35695831-37164130")
  print(ld_info)

  # STEP 1: Update unique_study_id to use ld block, as opposed to chr_bp
  # For debugging, use mc.cores = 1 to see errors more clearly
  # results <- parallel::mclapply(seq_len(nrow(ld_info)), mc.cores = 1, function(i) {
  results <- parallel::mclapply(seq_len(nrow(ld_info)), mc.cores = 30, function(i) {
    gc()
    return(tryCatch(
      {
        options(error = function() {
          message("Error in worker process:")
          traceback(20)
          quit(status = 1)
        })
        block_name <- ld_info[i, ]$block
        block_data <- ld_info[i, ]$ld_block_data

        message(paste("Processing", block_name))
        finemapped_studies_file <- paste0(block_data, "/finemapped_studies.tsv")
        finemapped_studies <- data.table::fread(finemapped_studies_file)

        if ("svg_file" %in% names(finemapped_studies)) {
          message(paste(block_name, "has already been updated"))
          return()
        }

        finemap_id <- as.numeric(sub(".*_", "", finemapped_studies$unique_study_id))
        finemapped_studies[, unique_study_id := paste0(study, "_", ld_block, "_", finemap_id)]

        imputed_studies_file <- paste0(block_data, "/imputed_studies.tsv")
        imputed_studies <- data.table::fread(imputed_studies_file)

        old_files <- finemapped_studies$file
        flattened_block_name <- gsub("[/-]", "_", block_name)
        new_files <- paste0(sub("(.*/).*", "\\1", old_files), flattened_block_name, "_", finemap_id, ".tsv.gz")

        svg_files <- c()
        new_lbf_data <- c()
        finemap_full_gwas_with_lbf_files <- paste0(
          extracted_study_dir,
          "/",
          finemapped_studies$study,
          "/finemapped/",
          flattened_block_name,
          "_with_lbf.tsv.gz"
        )
        imputed_full_gwas_with_lbf_files <- paste0(
          extracted_study_dir,
          "/",
          imputed_studies$study,
          "/finemapped/",
          flattened_block_name,
          "_with_lbf.tsv.gz"
        )

        for (i in seq_along(imputed_studies$study)) {
          study <- imputed_studies[i, ]
          imputed_gwas <- data.table::fread(study$file,
            select = c("SNP", "CHR", "BP", "EA", "OA", "EAF", "Z", "BETA", "SE", "IMPUTED")
          )
          new_lbf_data[[study$study]] <- imputed_gwas
        }

        message(paste(
          block_name,
          "processing",
          length(old_files),
          "finemapped studies, creating lbf, svgs, and updating coloc and finemapped results"
        ))
        for (file_index in seq_along(old_files)) {
          specific_study <- finemapped_studies[file_index, ]
          finemapped_gwas <- data.table::fread(old_files[file_index],
            select = c("SNP", "CHR", "BP", "Z", "SE", "EAF", "IMPUTED")
          )
          if (!"Z" %in% names(finemapped_gwas)) {
            message(paste("Z missing for", old_files[file_index]))
          }
          finemapped_gwas[, LBF := convert_z_to_lbf(Z, SE, EAF, specific_study$sample_size, specific_study$category)]

          finemap_id <- as.numeric(sub(".*_", "", specific_study$unique_study_id))
          lbf_id <- paste0("LBF_", finemap_id)
          lbf_data <- finemapped_gwas[, .(SNP, LBF)]
          data.table::setnames(lbf_data, "LBF", lbf_id)
          new_lbf_data[[specific_study$study]] <- new_lbf_data[[specific_study$study]][lbf_data, on = "SNP"]

          data.table::fwrite(finemapped_gwas, new_files[file_index], sep = "\t")

          finemapped_studies[file_index, file := new_files[file_index]]
          # svg_files[[file_index]] <- create_svg_for_ld_block(finemapped_studies[file_index, ])
        }

        message(paste(block_name, "updating finemapped and coloc studies files"))
        finemapped_studies[, svg_file := as.character(svg_files)]
        finemapped_studies[, file_with_lbfs := as.character(finemap_full_gwas_with_lbf_files)]

        # coloc_pairwise_results_file <- paste0(block_data, '/coloc_pairwise_results.tsv.gz')
        # coloc_pairwise_results <- data.table::fread(coloc_pairwise_results_file)
        # finemap_id_a <- as.numeric(sub(".*_", "", coloc_pairwise_results$unique_study_a))
        # finemap_id_b <- as.numeric(sub(".*_", "", coloc_pairwise_results$unique_study_b))
        # coloc_pairwise_results[, unique_study_a := paste0(study_a, '_', block_name, '_', finemap_id_a)]
        # coloc_pairwise_results[, unique_study_b := paste0(study_b, '_', block_name, '_', finemap_id_b)]

        message(paste(block_name, "updating", nrow(imputed_studies), "imputed studies with LBF data and saving"))
        results <- lapply(seq_along(imputed_studies$study), function(i) {
          study <- imputed_studies[i, ]
          data.table::fwrite(new_lbf_data[[study$study]], imputed_full_gwas_with_lbf_files[i], sep = "\t")
          return()
        })

        data.table::fwrite(finemapped_studies, paste0(finemapped_studies_file), sep = "\t")
        # data.table::fwrite(coloc_pairwise_results, paste0(coloc_pairwise_results_file), sep = '\t')
        return(paste("Success:", block_name))
      },
      error = function(e) {
        message(paste("Error in worker process for block", i, ":", e$message))
        message("Full error:")
        print(e)
        return(paste("Error:", e$message))
      }
    ))
  })

  # Check results for errors
  error_results <- results[sapply(results, function(x) grepl("^Error:", x))]
  if (length(error_results) > 0) {
    message("Errors found in worker processes:")
    print(error_results)
  }
  return()
  # Rprof(NULL)
  # prof_result <- summaryRprof()
  # saveRDS(prof_result, 'prof_result.rds')
}


create_svgs_for_all_phenotypes <- function() {
  source("../pipeline_steps/svg_helpers.R")

  studies <- vroom::vroom(glue::glue(results_dir, "latest/studies_processed.tsv.gz"), show_col_types = F) |>
    dplyr::filter(data_type == "phenotype" & variant_type == "common")
  # dplyr::filter(study_name == "ebi-a-GCST90002384")

  lapply(seq_len(nrow(studies)), function(i) {
    study <- studies[i, ]
    # if (file.exists(glue::glue('{study$extracted_location}/svgs/full.zip'))) return()

    print(study$study_name)
    dir.create(glue::glue("{study$extracted_location}/svgs"), showWarnings = F, recursive = T)

    vcf_file <- glue::glue("{study$extracted_location}/vcf/hg38.vcf.gz")
    bcf_query <- glue::glue(
      "/home/bcftools/bcftools query ",
      '--format "[%ID]\t[%CHROM]\t[%POS]\t[%LP]" ',
      "{vcf_file}"
    )
    gwas <- system(bcf_query, wait = T, intern = T)
    gwas <- data.table::fread(text = gwas)

    colnames(gwas) <- c("RSID", "CHR", "BP", "LP")
    gwas <- gwas |> dplyr::mutate(LP = as.numeric(LP))

    create_svgs_from_gwas(study, gwas)
    return(invisible(NULL))
  })
  return(invisible(NULL))
}

rejoin_imputed_studies <- function() {
  ld_blocks <- vroom::vroom("../pipeline_steps/data/ld_blocks.tsv")
  ld_info <- construct_ld_block(ld_blocks$ancestry, ld_blocks$chr, ld_blocks$start, ld_blocks$stop)
  ld_info <- dplyr::filter(ld_info, dir.exists(ld_block_data))

  blocks <- ld_info$block
  for (block in blocks) {
    backup_ld_block <- glue::glue("/local-scratch/projects/genotype-phenotype-map/backup/data/ld_blocks/{block}")
    ld_block <- glue::glue("/local-scratch/projects/genotype-phenotype-map/data/ld_blocks/{block}")

    imputed_studies_file <- glue::glue("{ld_block}/imputed_studies.tsv")
    file.copy(imputed_studies_file, glue::glue("{imputed_studies_file}.backup"))
    backup_imputed_studies_file <- glue::glue("{backup_ld_block}/imputed_studies.tsv")

    print(file.info(imputed_studies_file)$mtime)
    update_me <- file.info(imputed_studies_file)$mtime > "2025-06-20"
    print(update_me)
    if (!update_me) {
      print(paste("NO UPDATE NEEDED:", block))
      next
    }

    backup_imputed_studies <- vroom::vroom(backup_imputed_studies_file, show_col_types = F) |>
      dplyr::mutate(time_taken = as.character(time_taken))
    imputed_studies <- vroom::vroom(imputed_studies_file, show_col_types = F) |>
      dplyr::mutate(time_taken = as.character(time_taken))

    imputed_studies <- dplyr::bind_rows(imputed_studies, backup_imputed_studies) |>
      dplyr::distinct(study, .keep_all = T)

    message(glue::glue("Updated imputed {block}: {nrow(imputed_studies)}"))

    vroom::vroom_write(imputed_studies, imputed_studies_file)
  }
  return(invisible(NULL))
}

delete_still_bad_finemapped_studies <- function() {
  ld_blocks <- vroom::vroom("../pipeline_steps/data/ld_blocks.tsv")
  ld_info <- construct_ld_block(ld_blocks$ancestry, ld_blocks$chr, ld_blocks$start, ld_blocks$stop)
  ld_info <- dplyr::filter(ld_info, dir.exists(ld_block_data))

  blocks <- ld_info$ld_block_data
  for (ld_block in blocks) {
    print(ld_block)
    studies_to_remove <- vroom::vroom(glue::glue("{ld_block}/imputed_studies.tsv.backup"), show_col_types = F)
    finemapped_studies_file <- glue::glue("{ld_block}/finemapped_studies.tsv")
    finemapped_studies <- vroom::vroom(finemapped_studies_file, show_col_types = F)

    good_finemapped_studies <- dplyr::filter(finemapped_studies, !study %in% studies_to_remove$study)
    bad_finemapped_studies <- dplyr::filter(finemapped_studies, study %in% studies_to_remove$study)
    message(paste("removing", nrow(bad_finemapped_studies), "finemapped studies from", ld_block))
    vroom::vroom_write(good_finemapped_studies, finemapped_studies_file)

    pairwise_coloc_file <- glue::glue("{ld_block}/coloc_pairwise_results.tsv.gz")
    pairwise_coloc <- vroom::vroom(pairwise_coloc_file, show_col_types = F)
    good_pairwise_coloc <- dplyr::filter(
      pairwise_coloc, !study_a %in% studies_to_remove$study & !study_b %in% studies_to_remove$study
    )
    bad_pairwise_coloc <- dplyr::filter(
      pairwise_coloc, study_a %in% studies_to_remove$study | study_b %in% studies_to_remove$study
    )
    message(paste("removing", nrow(bad_pairwise_coloc), "pairwise coloc studies from", ld_block))
    vroom::vroom_write(good_pairwise_coloc, pairwise_coloc_file)
  }
  return(invisible(NULL))
}

delete_bad_imputations <- function() {
  ld_blocks <- vroom::vroom("../pipeline_steps/data/ld_blocks.tsv")
  ld_info <- construct_ld_block(ld_blocks$ancestry, ld_blocks$chr, ld_blocks$start, ld_blocks$stop)
  ld_info <- dplyr::filter(ld_info, dir.exists(ld_block_data))

  block_data <- c("/local-scratch/projects/genotype-phenotype-map/data/ld_blocks/EUR/1/101384499-103762931/")

  parallel::mclapply(block_data, mc.cores = 10, function(ld_block) {
    print(ld_block)
    imputed_studies_file <- glue::glue("{ld_block}/imputed_studies.tsv")
    imputed_studies <- vroom::vroom(imputed_studies_file, show_col_types = F)

    bad_studies <- apply(imputed_studies, 1, function(imputed_study) {
      file <- imputed_study[["file"]]
      if (!file.exists(file)) {
        return(data.frame(ld_block = ld_block, study = imputed_study[["study"]], file = file))
      }
      imputed_gwas <- vroom::vroom(file, show_col_types = F)
      if (sum(is.na(imputed_gwas$SE)) > 0 || any(imputed_gwas$SE <= 0)) {
        return(data.frame(ld_block = ld_block, study = imputed_study[["study"]], file = file))
      }
    }) |> dplyr::bind_rows()

    message(paste("removing", nrow(bad_studies), "imputed studies from", ld_block))
    good_imputed_studies <- dplyr::filter(imputed_studies, !file %in% bad_studies$file)
    vroom::vroom_write(good_imputed_studies, imputed_studies_file)

    finemapped_studies_file <- glue::glue("{ld_block}/finemapped_studies.tsv")
    finemapped_studies <- vroom::vroom(finemapped_studies_file, show_col_types = F)

    good_finemapped_studies <- dplyr::filter(finemapped_studies, !study %in% bad_studies$study)
    bad_finemapped_studies <- dplyr::filter(finemapped_studies, study %in% bad_studies$study)
    message(paste("removing", nrow(bad_finemapped_studies), "finemapped studies from", ld_block))
    vroom::vroom_write(good_finemapped_studies, finemapped_studies_file)

    pairwise_coloc_file <- glue::glue("{ld_block}/coloc_pairwise_results.tsv.gz")
    pairwise_coloc <- vroom::vroom(pairwise_coloc_file, show_col_types = F)
    good_pairwise_coloc <- dplyr::filter(
      pairwise_coloc,
      !study_a %in% bad_studies$study & !study_b %in% bad_studies$study
    )
    bad_pairwise_coloc <- dplyr::filter(pairwise_coloc, study_a %in% bad_studies$study | study_b %in% bad_studies$study)
    message(paste("removing", nrow(bad_pairwise_coloc), "pairwise coloc studies from", ld_block))
    vroom::vroom_write(good_pairwise_coloc, pairwise_coloc_file)
    return()
  })
  return()
}

update_existing_eafs <- function() {
  ld_block_matrices <- list()
  study_to_remove <- "/local-scratch/projects/genotype-phenotype-map/data/study/ebi-a-GCST000758/"

  all_studies <- Sys.glob(paste0(extracted_study_dir, "[a-zC-Z]*/"))
  print(length(all_studies))
  remove_to <- match(study_to_remove, all_studies)
  all_studies <- all_studies[-seq(remove_to)]
  print(length(all_studies))

  for (study in all_studies) {
    extracted_snps_file <- paste0(study, "/extracted_snps.tsv")
    if (!file.exists(extracted_snps_file)) {
      print(paste("EXTRACTION MISSING: ", extracted_snps_file))
      next
    }

    extracted_snps <- vroom::vroom(extracted_snps_file, show_col_types = F)
    if (nrow(extracted_snps) == 0) next

    apply(extracted_snps, 1, function(extraction) {
      ld_region <- extraction[["ld_region"]]
      if (ld_region == "//_") {
        return()
      }

      if (is.null(ld_block_matrices[[ld_region]])) {
        ld_block_matrices[[ld_region]] <- vroom::vroom(
          paste0(ld_reference_panel_dir, ld_region, ".tsv"),
          show_col_types = F
        )
      }

      ld_block_matrix <- ld_block_matrices[[ld_region]]

      original_gwas <- vroom::vroom(extraction[["file"]], show_col_types = F)

      if (!all(is.na(as.numeric(original_gwas$EAF)))) {
        print(paste("EAF already included: ", extraction[["file"]]))
        return()
      }

      imputed_study <- sub("original", "imputed", extraction[["file"]])
      if (!file.exists(imputed_study)) {
        print(paste("IMPUTATION MISSING: ", imputed_study))
        return()
      }
      imputed_gwas <- vroom::vroom(imputed_study, show_col_types = F)
      print(imputed_study)
      updated_imputed_gwas <- flippy_flippy(ld_block_matrix, imputed_gwas)
      vroom::vroom_write(updated_imputed_gwas, imputed_study)

      finemapped_studies <- sub("original", "finemapped", extraction[["file"]])
      finemap_file_prefix <- sub("\\..*", "", finemapped_studies)
      all_finemaps <- Sys.glob(paste0(finemap_file_prefix, "*"))

      if (length(all_finemaps) == 0) {
        print(paste("FINEMAPPING MISSING: ", finemap_file_prefix))
        return()
      }

      for (finemap_study in all_finemaps) {
        finemap_gwas <- vroom::vroom(finemap_study, show_col_types = F)
        print(finemap_study)
        updated_finemap_gwas <- flippy_flippy(ld_block_matrix, finemap_gwas)
        vroom::vroom_write(updated_finemap_gwas, finemap_study)
      }
      return()
    })
  }
  return()
}


standardise_z_scores <- function() {
  library(parallel)
  all_studies <- Sys.glob(paste0(extracted_study_dir, "*/"))
  study_to_remove <- "/local-scratch/projects/genotype-phenotype-map/data/study/GTEx-sqtl-cis-Cells-Cultured-fibroblasts-chr1-chr1:100921841:100961828:clu-44102:ENSG00000162695-11/" # nolint: line_length_linter.
  remove_to <- match(study_to_remove, all_studies)
  all_studies <- all_studies[-seq(remove_to)]
  print(length(all_studies))

  results <- mclapply(X = all_studies, mc.cores = 100, FUN = function(study) {
    print(paste("standardising", study))
    extracted_snps_file <- paste0(study, "/extracted_snps.tsv")
    if (!file.exists(extracted_snps_file)) {
      print(paste("EXTRACTION MISSING: ", extracted_snps_file))
      return()
    }
    extracted_snps <- vroom::vroom(extracted_snps_file, show_col_types = F)
    if (nrow(extracted_snps) == 0) {
      return()
    }

    apply(extracted_snps, 1, function(extraction) {
      standardised_imputed_study <- sub("original", "imputed", extraction[["file"]])
      old_imputed_study <- sub("study/", "study_fix/study/", standardised_imputed_study)
      if (!file.exists(standardised_imputed_study)) {
        print(paste("IMPUTATION MISSING: ", standardised_imputed_study))
        return()
      }
      standardised_imputed_gwas <- vroom::vroom(standardised_imputed_study, show_col_types = F)
      old_imputed_gwas <- vroom::vroom(old_imputed_study, show_col_types = F)
      updated_imputed_gwas <- fix_z_score_direction(old_imputed_gwas, standardised_imputed_gwas)

      if (any((updated_imputed_gwas$Z * updated_imputed_gwas$BETA) < 0, na.rm = T)) {
        quit(paste("ERROR: ", standardised_imputed_study, "BETA and Z still dont match"))
      }

      vroom::vroom_write(updated_imputed_gwas, standardised_imputed_study)

      finemapped_studies <- sub("original", "finemapped", extraction[["file"]])
      finemap_file_prefix <- sub("\\..*", "", finemapped_studies)
      all_finemaps <- Sys.glob(paste0(finemap_file_prefix, "*"))

      if (length(all_finemaps) == 0) {
        print(paste("FINEMAPPING MISSING: ", finemap_file_prefix))
        return()
      }

      for (standardised_finemap_study in all_finemaps) {
        old_finemapped_study <- sub("study/", "study_fix/study/", standardised_finemap_study)

        standardised_finemap_gwas <- vroom::vroom(standardised_finemap_study, show_col_types = F)
        old_finemap_gwas <- vroom::vroom(old_finemapped_study, show_col_types = F)

        updated_finemap_gwas <- fix_z_score_direction(old_finemap_gwas, standardised_finemap_gwas)
        if (any((updated_finemap_gwas$Z * updated_finemap_gwas$BETA) < 0, na.rm = T)) {
          quit(paste("ERROR: ", standardised_finemap_gwas, "BETA and Z still dont match"))
        }
        vroom::vroom_write(updated_finemap_gwas, standardised_finemap_study)
      }
      return()
    })
    return()
  })
  return()
}

fix_z_score_direction <- function(old_gwas, standardised_gwas) {
  to_flip <- (old_gwas$EA > old_gwas$OA) & (!old_gwas$EA %in% c("D", "I"))
  print(paste("flipping", sum(to_flip), "z scores"))
  if (any(to_flip)) {
    standardised_gwas$Z[to_flip] <- -1 * standardised_gwas$Z[to_flip]
  }

  return(standardised_gwas)
}


update_studies_processed <- function() {
  studies_processed_file <- paste0(paste0(results_dir, "studies_processed.tsv"))
  studies_processed <- vroom::vroom(studies_processed_file, delim = "\t", show_col_types = F)
  studies_processed$study_name <- sub("\\.", "-", studies_processed$study_name)

  studies_processed$study_name <- sub(
    "BrainMeta-cis-eqtl-BrainMeta-cis-eQTL",
    "BrainMeta-cis-eQTL",
    studies_processed$study_name
  )
  studies_processed$study_name <- sub("Brain-eMeta-Brain-eMeta", "Brain-eMeta-full", studies_processed$study_name)

  vroom::vroom_write(studies_processed, studies_processed_file)
  return()
}

update_extracted_snps <- function() {
  brain_studies <- Sys.glob(paste0(extracted_study_dir, "BrainMeta*/"))

  dont_print <- lapply(brain_studies, function(study) {
    extracted_snps_file <- paste0(study, "/extracted_snps.tsv")
    extracted_snps <- vroom::vroom(extracted_snps_file, show_col_types = F)
    if (nrow(extracted_snps) > 0) print(extracted_snps_file)

    extracted_snps$file <- sub("BrainMeta-cis-eqtl-BrainMeta-cis-eQTL", "BrainMeta-cis-eQTL", extracted_snps$file)
    extracted_snps$file <- sub("Brain-eMeta-Brain-eMeta", "Brain-eMeta-full", extracted_snps$file)
    extracted_snps$file[grepl(
      "BrainMeta-cis-eQTL",
      extracted_snps$file
    )] <- sub("\\.", "-", extracted_snps$file[grepl("BrainMeta-cis-eQTL", extracted_snps$file)])

    vroom::vroom_write(extracted_snps, extracted_snps_file)
    return()
  })
  return()
}


brain_meta_updates <- function() {
  ld_regions <- vroom::vroom("../pipeline_steps/data/ld_regions.tsv")
  ld_info <- construct_ld_block(ld_regions$ancestry, ld_regions$chr, ld_regions$start, ld_regions$stop)

  lapply(ld_info$ld_block_data, function(ld_block) {
    print(ld_block)
    extracted_studies_file <- paste0(ld_block, "/extracted_studies.tsv")
    if (!file.exists(extracted_studies_file)) {
      return()
    }

    extracted_studies <- vroom::vroom(extracted_studies_file, show_col_types = F)

    extracted_studies$ancestry <- "EUR"
    extracted_studies$study <- sub(
      "BrainMeta-cis-eqtl-BrainMeta-cis-eQTL",
      "BrainMeta-cis-eQTL",
      extracted_studies$study
    )
    extracted_studies$study <- sub("Brain-eMeta-Brain-eMeta", "Brain-eMeta-full", extracted_studies$study)
    extracted_studies$study <- sub("\\.", "-", extracted_studies$study)
    extracted_studies$file <- sub("BrainMeta-cis-eqtl-BrainMeta-cis-eQTL", "BrainMeta-cis-eQTL", extracted_studies$file)
    extracted_studies$file <- sub("Brain-eMeta-Brain-eMeta", "Brain-eMeta-full", extracted_studies$file)
    extracted_studies$file[grepl(
      "BrainMeta-cis-eQTL",
      extracted_studies$file
    )] <- sub("\\.", "-", extracted_studies$file[grepl("BrainMeta-cis-eQTL", extracted_studies$file)])

    vroom::vroom_write(extracted_studies, extracted_studies_file)


    imputed_studies_file <- paste0(ld_block, "/imputed_studies.tsv")
    if (!file.exists(imputed_studies_file)) {
      return()
    }

    imputed_studies <- vroom::vroom(imputed_studies_file, show_col_types = F)

    imputed_studies$ancestry <- "EUR"
    imputed_studies$study <- sub("BrainMeta-cis-eqtl-BrainMeta-cis-eQTL", "BrainMeta-cis-eQTL", imputed_studies$study)
    imputed_studies$study <- sub("Brain-eMeta-Brain-eMeta", "Brain-eMeta-full", imputed_studies$study)
    imputed_studies$study <- sub("\\.", "-", imputed_studies$study)
    imputed_studies$file <- sub("BrainMeta-cis-eqtl-BrainMeta-cis-eQTL", "BrainMeta-cis-eQTL", imputed_studies$file)
    imputed_studies$file <- sub("Brain-eMeta-Brain-eMeta", "Brain-eMeta-full", imputed_studies$file)
    imputed_studies$file[grepl(
      "BrainMeta-cis-eQTL",
      imputed_studies$file
    )] <- sub("\\.", "-", imputed_studies$file[grepl("BrainMeta-cis-eQTL", imputed_studies$file)])

    vroom::vroom_write(imputed_studies, imputed_studies_file)

    finemapped_studies_file <- paste0(ld_block, "/finemapped_studies.tsv")
    if (!file.exists(finemapped_studies_file)) {
      return()
    }

    finemapped_studies <- vroom::vroom(finemapped_studies_file, show_col_types = F)

    finemapped_studies$ancestry <- "EUR"
    finemapped_studies$study <- sub(
      "BrainMeta-cis-eqtl-BrainMeta-cis-eQTL",
      "BrainMeta-cis-eQTL",
      finemapped_studies$study
    )
    finemapped_studies$study <- sub("Brain-eMeta-Brain-eMeta", "Brain-eMeta-full", finemapped_studies$study)
    finemapped_studies$study <- sub("\\.", "-", finemapped_studies$study)
    finemapped_studies$unique_study_id <- sub(
      "BrainMeta-cis-eqtl-BrainMeta-cis-eQTL",
      "BrainMeta-cis-eQTL",
      finemapped_studies$unique_study_id
    )
    finemapped_studies$unique_study_id <- sub(
      "Brain-eMeta-Brain-eMeta",
      "Brain-eMeta-full",
      finemapped_studies$unique_study_id
    )
    finemapped_studies$unique_study_id <- sub("\\.", "-", finemapped_studies$unique_study_id)
    finemapped_studies$unique_study_id <- gsub(" ", "", finemapped_studies$unique_study_id)
    fill_study_name <- grepl("^EUR", finemapped_studies$unique_study_id)
    finemapped_studies$unique_study_id[fill_study_name] <- paste0(
      finemapped_studies$study[fill_study_name],
      "_",
      finemapped_studies$unique_study_id[fill_study_name]
    )

    finemapped_studies$file <- sub(
      "BrainMeta-cis-eqtl-BrainMeta-cis-eQTL",
      "BrainMeta-cis-eQTL",
      finemapped_studies$file
    )
    finemapped_studies$file <- sub("Brain-eMeta-Brain-eMeta", "Brain-eMeta-full", finemapped_studies$file)

    chr_bp_version <- sub(
      "^[^_]*_",
      "",
      finemapped_studies$unique_study_id[grepl("BrainMeta-cis-eQTL", finemapped_studies$unique_study_id)]
    )
    move_to <- paste0(
      extracted_study_dir,
      finemapped_studies$study[grepl("BrainMeta-cis-eQTL", finemapped_studies$unique_study_id)],
      "/finemapped/EUR_",
      chr_bp_version,
      ".tsv.gz"
    )
    file.copy(finemapped_studies$file[grepl("BrainMeta-cis-eQTL", finemapped_studies$file)], move_to)

    finemapped_studies$file[grepl("BrainMeta-cis-eQTL", finemapped_studies$file)] <- move_to

    fill_eur <- !grepl("EUR", finemapped_studies$unique_study_id)
    finemapped_studies$unique_study_id[fill_eur] <- sub("_", "_EUR_", finemapped_studies$unique_study_id[fill_eur])

    no_min_p_finemapped <- finemapped_studies[is.na(finemapped_studies$min_p), ]
    actual_min_ps <- lapply(no_min_p_finemapped$file, function(file) {
      return(min(vroom::vroom(file, show_col_types = F)$P))
    })
    finemapped_studies$min_p[is.na(finemapped_studies$min_p)] <- unlist(actual_min_ps)

    vroom::vroom_write(finemapped_studies, finemapped_studies_file)
    return()
  })
  return()
}

update_imputed_studies <- function() {
  ld_regions <- vroom::vroom("../pipeline_steps/data/ld_regions.tsv")
  ld_info <- construct_ld_block(ld_regions$ancestry, ld_regions$chr, ld_regions$start, ld_regions$stop)

  lapply(ld_info$ld_block_data, function(ld_block) {
    extracted_studies_file <- paste0(ld_block, "/extracted_studies.tsv")
    imputed_studies_file <- paste0(ld_block, "/imputed_studies.tsv")
    if (!file.exists(extracted_studies_file)) {
      return()
    }

    extracted_studies <- vroom::vroom(extracted_studies_file, show_col_types = F)
    imputed_studies <- vroom::vroom(imputed_studies_file, col_types = vroom::cols(
      chr = vroom::col_character(),
      bp = vroom::col_number(),
      p_value_threshold = vroom::col_number(),
      eaf_from_reference_panel = vroom::col_logical(),
      sample_size = vroom::col_number()
    ), show_col_types = F)
    if (nrow(extracted_studies) == nrow(imputed_studies)) {
      return()
    }
    print(paste(ld_block, "sizes are different...", nrow(extracted_studies) - nrow(imputed_studies)))

    add_to_imputed <- apply(extracted_studies, 1, function(extracted_study) {
      already_in_imputed <- dplyr::filter(imputed_studies, study == extracted_study["study"])
      if (nrow(already_in_imputed) > 0) {
        return()
      }

      imputed_study <- sub("original", "imputed", extracted_study[["file"]])
      if (!file.exists(imputed_study)) {
        return()
      }

      eaf_from_reference_panel <- ifelse(grepl("GTEx", extracted_study["study"]), TRUE, NA)

      return(data.frame(
        study = extracted_study["study"],
        file = extracted_study["file"],
        chr = extracted_study["chr"],
        bp = as.numeric(extracted_study["bp"]),
        p_value_threshold = as.numeric(extracted_study["p_value_threshold"]),
        sample_size = as.numeric(extracted_study["sample_size"]),
        cis_trans = extracted_study["cis_trans"],
        rows_imputed = NA, eaf_from_reference_panel = eaf_from_reference_panel
      ))
    }) |> dplyr::bind_rows()

    imputed_studies <- dplyr::bind_rows(imputed_studies, add_to_imputed) |> dplyr::distinct()
    vroom::vroom_write(imputed_studies, imputed_studies_file)
    return()
  })
  return()
}

update_finemapped_studies <- function() {
  ld_regions <- vroom::vroom("../pipeline_steps/data/ld_regions.tsv")
  ld_regions <- ld_regions[ld_regions$chr >= 12, ]
  ld_info <- construct_ld_block(ld_regions$ancestry, ld_regions$chr, ld_regions$start, ld_regions$stop)

  lapply(ld_info$ld_block_data, function(ld_block) {
    finemapped_studies_file <- paste0(ld_block, "/finemapped_studies.tsv")
    imputed_studies_file <- paste0(ld_block, "/imputed_studies.tsv")
    if (!file.exists(finemapped_studies_file)) {
      return()
    }
    print(ld_block)

    imputed_studies <- vroom::vroom(imputed_studies_file, show_col_types = F)
    finemapped_studies <- vroom::vroom(finemapped_studies_file, col_types = vroom::cols(
      chr = vroom::col_character(),
      bp = vroom::col_number(),
      min_p = vroom::col_number(),
      p_value_threshold = vroom::col_number(),
      sample_size = vroom::col_number()
    ), show_col_types = F)

    if (nrow(imputed_studies) == 0) {
      return()
    }
    add_to_finemapped <- apply(imputed_studies, 1, function(imputed_study) {
      finemap_file_prefix <- sub("imputed", "finemapped", imputed_study[["file"]])
      finemap_file_prefix <- sub("\\..*", "", finemap_file_prefix)
      finemap_files <- Sys.glob(paste0(finemap_file_prefix, "*"))

      if (length(finemap_files) == 0) {
        return()
      }
      if (nrow(dplyr::filter(finemapped_studies, grepl(finemap_file_prefix, file))) > 0) {
        return()
      }

      unique_study_ids <- file_prefix(finemap_files)
      message <- if (length(finemap_files) > 1) "success" else NA

      bps <- lapply(finemap_files, function(file) {
        finemap <- vroom::vroom(file, show_col_types = F)
        min <- finemap[which.min(finemap$P), ]
        return(data.frame(bp = min$BP, min_p = min$P))
      }) |> dplyr::bind_rows()

      print(unique_study_ids)
      print(finemap_files)
      print(bps$bp)
      print(bps$min_p)

      return(data.frame(
        study = imputed_study["study"],
        unique_study_id = unique_study_ids,
        file = finemap_files,
        chr = imputed_study["chr"],
        bp = as.numeric(bps$bp),
        p_value_threshold = as.numeric(imputed_study["p_value_threshold"]),
        sample_size = as.numeric(imputed_study["sample_size"]),
        cis_trans = imputed_study["cis_trans"],
        message = message, min_p = bps$min_p
      ))
    }) |>
      dplyr::bind_rows()

    if (nrow(add_to_finemapped) == 0) {
      return()
    }
    print(paste("missing imputed from finemapped", nrow(add_to_finemapped)))

    finemapped_studies <- dplyr::bind_rows(finemapped_studies, add_to_finemapped) |> dplyr::distinct()
    vroom::vroom_write(finemapped_studies, finemapped_studies_file)
    return()
  })
  return()
}

populate_json_for_besd <- function() {
  gtex <- vroom::vroom("/local-scratch/projects/genotype-phenotype-map/GTEx_portal.csv")

  gtex$Tissue <- gsub("[ -]", "_", gtex$Tissue)
  gtex$Tissue <- gsub("_+", "_", gtex$Tissue)
  gtex$Tissue <- gsub("[()]", "", gtex$Tissue)

  apply(gtex, 1, function(tissue) {
    metadata <- list(
      sample_size = as.numeric(tissue["# RNASeq and Genotyped samples"]),
      ancestry = "EUR",
      cis_trans = "cis",
      category = "continuous"
    )
    metadata_json <- jsonlite::toJSON(metadata, auto_unbox = T, pretty = T)
    file_name <- paste0("/local-scratch/data/GTEx-cis/", tissue["Tissue"], ".json")
    write(metadata_json, file_name)
    return()
  })
  return()
}

remove_imputed_and_finemapped_results_from_pipeline <- function(study_pattern) {
  NUM_PARALLEL_JOBS <- 100

  all_studies <- Sys.glob(paste0(extracted_study_dir, study_pattern, "*/"))
  # study_to_remove <- '/local-scratch/projects/genotype-phenotype-map/data/study/study-name'
  # remove_to <- match(study_to_remove, all_studies)
  # all_studies <- all_studies[-seq(remove_to)]

  results <- parallel::mclapply(X = all_studies, mc.cores = NUM_PARALLEL_JOBS, FUN = function(study_dir) {
    print(paste("changing:", study_dir))
    study_id <- basename(study_dir)
    extracted_snps_file <- paste0(study_dir, "/extracted_snps.tsv")
    if (!file.exists(extracted_snps_file)) {
      print(paste("EXTRACTION MISSING: ", extracted_snps_file))
      return()
    }
    extracted_snps <- vroom::vroom(extracted_snps_file, show_col_types = F)
    if (nrow(extracted_snps) == 0) {
      return()
    }

    apply(extracted_snps, 1, function(extraction) {
      ld_block <- paste0(ld_block_data_dir, extraction[["ld_region"]])

      imputed_studies_file <- paste0(ld_block, "/imputed_studies.tsv")
      if (!file.exists(imputed_studies_file)) {
        print(paste("IMPUTED STUDIES FILE MISSING:", imputed_studies_file))
        return()
      }

      imputed_studies <- vroom::vroom(imputed_studies_file, show_col_types = F)
      entries <- nrow(imputed_studies)
      imputed_studies <- dplyr::filter(imputed_studies, !study == study_id)
      print(paste("removed", entries - nrow(imputed_studies), "rows from imputed_studies"))
      vroom::vroom_write(imputed_studies, imputed_studies_file)
      imputed_files <- Sys.glob(paste0(study_dir, "/imputed/*"))
      file.remove(imputed_files)


      finemapped_studies_file <- paste0(ld_block, "/finemapped_studies.tsv")
      if (!file.exists(finemapped_studies_file)) {
        print(paste("FINEMAPPED STUDIES FILE MISSING:", finemapped_studies_file))
        return()
      }

      finemapped_studies <- vroom::vroom(finemapped_studies_file, show_col_types = F)
      entries <- nrow(finemapped_studies)
      finemapped_studies <- dplyr::filter(finemapped_studies, !study == study_id)
      print(paste("removed", entries - nrow(finemapped_studies), "rows from finemapped_studies"))
      vroom::vroom_write(finemapped_studies, finemapped_studies_file)
      return()
    })
    finemapped_files <- Sys.glob(paste0(study_dir, "/finemapped/*"))
    file.remove(finemapped_files)
    return()
  })
  return()
}

print_alleles_to_flip <- function() {
  flipped_dir <- paste0(thousand_genomes_dir, "/flipped/")
  bims <- Sys.glob(paste0(flipped_dir, "*.bim"))

  lapply(bims, function(bim_file) {
    bim <- vroom::vroom(bim_file, col_names = F)
    print(names(bim))
    bim <- bim[(bim$X5 > bim$X6), ]
    vroom::vroom_write(dplyr::select(bim, X2), paste0(bim_file, "_toflip"), col_names = F)
    return()
  })
  return()
}


populate_beta_in_finemapped_files()
