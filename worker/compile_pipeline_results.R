#' Compile results from pipeline
#'
#' This function compiles the results from the pipeline into a series of files
#' * compiled_coloc_pairwise_results.tsv
#' * compiled_study_extractions.tsv
#' * compiled_coloc_clustered_results.tsv
#' * compiled_associations.tsv
#' * gwas_with_lbfs.tsv.gz
#'
#' The files are written to the extracted_study_dir directory and also uploaded to the Oracle bucket.
#'
#' @param gwas_info A list containing the GWAS information.
#' @return A list containing the compiled results.
compile_results <- function(gwas_info) {
  ld_block_dirs <- list.dirs(ld_block_data_dir, recursive = TRUE, full.names = TRUE) |>
    (\(dirs) dirs[!dirs %in% dirname(dirs[-1])])()

  compiled_coloc_pairwise_results_file <- glue::glue("{extracted_study_dir}/compiled_coloc_pairwise_results.tsv")
  compiled_study_extractions_file <- glue::glue("{extracted_study_dir}/compiled_extracted_studies.tsv")
  compiled_coloc_clustered_results_file <- glue::glue("{extracted_study_dir}/compiled_coloc_clustered_results.tsv")
  compiled_associations_file <- glue::glue("{extracted_study_dir}/compiled_associations.tsv")
  lbfs_concatenated_file <- glue::glue("{extracted_study_dir}/gwas_with_lbfs.tsv.gz")

  finemapped_studies_files <- Filter(
    function(file) file.exists(file),
    glue::glue("{ld_block_dirs}/finemapped_studies.tsv")
  )
  compare_guids <- gwas_info$metadata$gwas_upload_ids_to_compare
  if (is.null(compare_guids)) compare_guids <- gwas_info$metadata$compare_with_upload_guids
  if (is.null(compare_guids)) compare_guids <- character(0)
  compare_guids <- as.character(unlist(compare_guids))
  compare_guids <- setdiff(compare_guids, gwas_info$metadata$guid)
  if (length(compare_guids) > 0) {
    ld_blocks <- sub(paste0("^", ld_block_data_dir), "", ld_block_dirs)
    compare_files <- as.character(outer(
      compare_guids,
      ld_blocks,
      function(guid, block) glue::glue("{gwas_upload_dir}ld_blocks/gwas_upload/{guid}/{block}/finemapped_studies.tsv")
    ))
    finemapped_studies_files <- c(finemapped_studies_files, Filter(file.exists, compare_files))
  }

  if (length(finemapped_studies_files) > 0) {
    all_finemapped_studies <- lapply(finemapped_studies_files, function(file) {
      se_result <- data.table::fread(
        file,
        showProgress = FALSE,
        colClasses = list(character = c("chr", "snp", "study", "unique_study_id"))
      ) |>
        dplyr::filter(min_p <= gwas_info$metadata$p_value_threshold)
      return(se_result)
    }) |>
      data.table::rbindlist(fill = TRUE)
    study_extractions <- all_finemapped_studies |> dplyr::filter(study == gwas_info$metadata$guid)
  } else {
    flog.warn(paste(gwas_info$metadata$guid, "No finemapped study files found"))
    study_extractions <- data.table::data.table()
    all_finemapped_studies <- data.table::data.table()
  }

  all_snps_in_ld_blocks <- study_extractions |>
    dplyr::distinct(snp) |>
    dplyr::pull(snp)

  snp_annotations <- data.table::fread(
    file.path(variant_annotation_dir, "vep_annotations_hg38.tsv.gz"),
    select = c("snp", "rsid", "display_snp"),
    showProgress = FALSE,
    colClasses = list(character = c("snp", "rsid", "display_snp"))
  ) |>
    dplyr::mutate(snp = trimws(snp)) |>
    dplyr::filter(snp %in% all_snps_in_ld_blocks)

  if (nrow(study_extractions) > 0) {
    study_extractions <- merge(study_extractions, snp_annotations, by = "snp", all.x = TRUE)
  }

  study_extractions <- study_extractions |>
    dplyr::select(study, unique_study_id, snp, chr, bp, min_p, ld_block, file, file_with_lbfs)
  vroom::vroom_write(study_extractions, compiled_study_extractions_file)

  coloc_clustered_results_files <- Filter(
    function(file) file.exists(file),
    glue::glue("{ld_block_dirs}/coloc_clustered_results.tsv.gz")
  )
  if (length(coloc_clustered_results_files) > 0) {
    coloc_clustered_results <- lapply(coloc_clustered_results_files, function(file) {
      return(data.table::fread(file, showProgress = FALSE))
    }) |>
      data.table::rbindlist(fill = TRUE)

    if (nrow(coloc_clustered_results) > 0) {
      worker_unique_study_ids <- study_extractions$unique_study_id
      components_with_worker <- coloc_clustered_results |>
        dplyr::filter(unique_study_id %in% worker_unique_study_ids) |>
        dplyr::distinct(ld_block, component)
      coloc_clustered_results <- coloc_clustered_results |>
        dplyr::inner_join(components_with_worker, by = c("ld_block", "component")) |>
        dplyr::group_by(ld_block, component) |>
        dplyr::mutate(coloc_group_id = dplyr::cur_group_id()) |>
        dplyr::ungroup() |>
        dplyr::arrange(coloc_group_id) |>
        dplyr::select(unique_study_id, coloc_group_id, snp, ld_block, h4_connectedness, h3_connectedness) |>
        data.table::as.data.table()
    }
  } else {
    coloc_clustered_results <- data.table::data.table()
  }
  vroom::vroom_write(coloc_clustered_results, compiled_coloc_clustered_results_file)

  coloc_pairwise_results_files <- Filter(
    function(file) file.exists(file),
    glue::glue("{ld_block_dirs}/coloc_pairwise_results.tsv.gz")
  )
  if (length(compare_guids) > 0) {
    ld_blocks <- sub(paste0("^", ld_block_data_dir), "", ld_block_dirs)
    compare_pairwise_files <- as.character(outer(
      compare_guids,
      ld_blocks,
      function(guid, block) {
        return(paste0(gwas_upload_dir, "ld_blocks/gwas_upload/", guid, "/", block, "/coloc_pairwise_results.tsv.gz"))
      }
    ))
    coloc_pairwise_results_files <- c(coloc_pairwise_results_files, Filter(file.exists, compare_pairwise_files))
  }

  if (length(coloc_pairwise_results_files) > 0) {
    coloc_pairwise_results <- lapply(coloc_pairwise_results_files, function(file) {
      cp_result <- data.table::fread(file, showProgress = FALSE) |>
        dplyr::filter(study_a == gwas_info$metadata$guid | study_b == gwas_info$metadata$guid) |>
        dplyr::filter(PP.H4.abf >= posterior_prob_threshold_minimum & ignore == FALSE) |>
        dplyr::select(
          unique_study_a,
          unique_study_b,
          PP.H3.abf,
          PP.H4.abf,
          ld_block,
          false_positive,
          false_negative,
          ignore
        )
      return(cp_result)
    }) |>
      data.table::rbindlist(fill = TRUE)

    if (nrow(coloc_pairwise_results) > 0) {
      data.table::setnames(coloc_pairwise_results,
        old = c("unique_study_a", "unique_study_b", "PP.H3.abf", "PP.H4.abf"),
        new = c("unique_study_id_a", "unique_study_id_b", "h3", "h4")
      )
    }
  } else {
    coloc_pairwise_results <- data.table::data.table()
  }
  vroom::vroom_write(coloc_pairwise_results, compiled_coloc_pairwise_results_file)

  concatenate_file_with_lbfs(gwas_info, study_extractions)

  associations <- find_associations_for_coloc_clustered_snps(
    gwas_info,
    coloc_clustered_results,
    all_finemapped_studies,
    snp_annotations
  )

  vroom::vroom_write(associations, compiled_associations_file)

  return(list(
    coloc_pairwise_results = coloc_pairwise_results,
    study_extractions = study_extractions,
    coloc_clustered_results = coloc_clustered_results,
    associations = associations
  ))
}

#' Find associations for coloc clustered SNPs
#'
#' This function finds the associations for the coloc clustered SNPs.
#'
#' @param gwas_info A list containing the GWAS information.
#' @param coloc_clustered_results A data frame containing the coloc clustered results.
#' @param all_finemapped_studies A data frame containing the all finemapped studies.
#' @param snp_annotations A data frame containing the SNP annotations.
#' @return A data frame containing the associations.
find_associations_for_coloc_clustered_snps <- function(
  gwas_info,
  coloc_clustered_results,
  all_finemapped_studies,
  snp_annotations
) {
  flog.info(paste(gwas_info$metadata$guid, "Finding associations for coloc clustered SNPs"))
  if (nrow(coloc_clustered_results) == 0 || nrow(all_finemapped_studies) == 0) {
    return(data.frame())
  }

  finemapped_studies <- all_finemapped_studies |>
    dplyr::select(unique_study_id, study, file, ld_block) |>
    dplyr::filter(unique_study_id %in% coloc_clustered_results$unique_study_id)

  coloc_snps <- coloc_clustered_results |>
    dplyr::filter(unique_study_id %in% finemapped_studies$unique_study_id) |>
    dplyr::select(unique_study_id, snp, ld_block)

  if (nrow(coloc_snps) == 0 || nrow(finemapped_studies) == 0) {
    return(data.frame())
  }

  coloc_snps_with_files <- coloc_snps |>
    dplyr::left_join(finemapped_studies, by = c("unique_study_id", "ld_block"))

  files_to_read <- coloc_snps_with_files |>
    dplyr::group_by(file) |>
    dplyr::group_split()

  snp_mapping <- data.table::as.data.table(
    dplyr::distinct(coloc_snps_with_files, snp, study)
  )
  data.table::setkey(snp_mapping, snp)

  process_one_assoc_file <- function(file_group) {
    file_path <- file_group$file[1]
    study_for_file <- file_group$study[1]
    target_snps <- unique(file_group$snp)

    if (!grepl("^/", file_path) && !grepl("^study", file_path)) {
      file_path <- file.path(data_dir, file_path)
    } else if (grepl("^study", file_path)) {
      file_path <- file.path(data_dir, file_path)
    }

    if (!file.exists(file_path)) {
      return(NULL)
    }

    avail_cols <- names(data.table::fread(file_path, nrows = 0, showProgress = FALSE))
    required <- c("SNP", "BETA", "SE", "EAF")
    if (!all(required %in% avail_cols)) {
      return(NULL)
    }
    cols_to_read <- c(required, if ("IMPUTED" %in% avail_cols) "IMPUTED")

    gwas <- data.table::fread(
      file_path,
      select = cols_to_read,
      showProgress = FALSE,
      nThread = 1
    )

    gwas <- gwas[SNP %in% target_snps]
    if (nrow(gwas) == 0) {
      return(NULL)
    }

    gwas[, p := beta_se_to_p(BETA, SE)]
    if (!"IMPUTED" %in% colnames(gwas)) gwas[, IMPUTED := FALSE]
    snp_mapping_this_study <- snp_mapping[study == study_for_file]
    gwas <- gwas[snp_mapping_this_study, on = c(SNP = "snp"), nomatch = 0]

    if (nrow(gwas) == 0) {
      return(NULL)
    }
    return(gwas[, .(SNP, study, BETA, SE, p, EAF, IMPUTED)])
  }

  assoc_chunks <- parallel::mclapply(
    files_to_read,
    process_one_assoc_file,
    mc.cores = parallel_block_processing
  )
  valid_chunks <- Filter(Negate(is.null), assoc_chunks)
  if (length(valid_chunks) == 0) {
    return(data.frame())
  }
  associations <- data.table::rbindlist(valid_chunks)

  if (nrow(associations) == 0) {
    return(data.frame())
  }

  associations <- unique(associations[, .(
    snp = SNP,
    study_name = study,
    beta = BETA,
    se = SE,
    p,
    eaf = EAF,
    imputed = IMPUTED
  )], by = c("study_name", "snp"))

  associations <- stats::na.omit(associations)

  flog.info(paste(
    gwas_info$metadata$guid,
    "Extracted",
    nrow(associations),
    "associations for coloc clustered SNPs"
  ))
  return(as.data.frame(associations))
}

#' Concatenate file with LBFs
#'
#' This function concatenates the file with LBFs and writes to the extracted_study_dir directory.
#'
#' @param gwas_info A list containing the GWAS information.
#' @param study_extractions A data frame containing the study extractions.
#' @return NULL
concatenate_file_with_lbfs <- function(gwas_info, study_extractions) {
  lbfs_concatenated_file <- glue::glue("{extracted_study_dir}/gwas_with_lbfs.tsv.gz")
  if (nrow(study_extractions) > 0) {
    flog.info(paste(gwas_info$metadata$guid, "Concatenating file_with_lbfs files"))

    lbf_files <- unique(study_extractions$file_with_lbfs)
    lbf_files <- lbf_files[!is.na(lbf_files)]

    if (length(lbf_files) > 0) {
      if (file.exists(lbfs_concatenated_file)) unlink(lbfs_concatenated_file)

      for (i in seq_along(lbf_files)) {
        file_path <- glue::glue("{data_dir}/{lbf_files[i]}")
        if (!file.exists(file_path)) next

        dt <- data.table::fread(file_path, showProgress = FALSE)
        if (nrow(dt) > 0) {
          data.table::fwrite(dt, lbfs_concatenated_file,
            append = TRUE,
            compress = "gzip", sep = "\t"
          )
        }
        rm(dt)
        gc(verbose = FALSE)
      }
    }
  }
  return(NULL)
}
