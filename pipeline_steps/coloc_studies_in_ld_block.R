source("constants.R")
library(dplyr)
dplyr.summarise.inform <- FALSE

bp_range <- 50000

# Wallace 2020-style defaults: interpret priors at a reference locus size,
# then rescale to each harmonised locus so locus-level priors stay comparable.
wallace_reference_nsnps <- 1000
coloc_base_p1 <- 1e-4
coloc_base_p2 <- 1e-4
coloc_base_p12 <- 5e-6

parser <- argparser::arg_parser("Colocalise studies per LD block")
parser <- argparser::add_argument(
  parser,
  "--ld_block",
  help = "LD block that the studies are in",
  type = "character"
)
parser <- argparser::add_argument(
  parser,
  "--completed_output_file",
  help = "Sentinel file to write once coloc is complete",
  type = "character"
)
parser <- argparser::add_argument(
  parser,
  "--block_list",
  help = "CSV of studies to exclude from colocalisation (columns: id_pattern, cis_trans)",
  type = "character",
  default = NA
)
parser <- argparser::add_argument(
  parser,
  "--worker_guid",
  help = "Worker GUID",
  type = "character",
  default = NA
)
parser <- argparser::add_argument(
  parser,
  "--worker_p_value_threshold",
  help = "Worker p-value threshold",
  type = "numeric",
  default = NA
)
parser <- argparser::add_argument(
  parser,
  "--gwas_upload_ids_to_compare",
  help = "Comma-separated list of GWAS upload GUIDs to include in coloc comparison",
  type = "character",
  default = NA
)

args <- argparser::parse_args(parser)

main <- function() {
  start_time <- Sys.time()
  if (!is.na(args$worker_guid)) {
    update_directories_for_worker(args$worker_guid)
  }

  ld_info <- ld_block_dirs(args$ld_block)
  coloc_results_file <- glue::glue("{ld_info$ld_block_data}/coloc_pairwise_results.tsv.gz")
  if (file.exists(coloc_results_file)) {
    coloc_results <- vroom::vroom(coloc_results_file, delim = "\t", show_col_types = F)
  } else {
    coloc_results <- data.frame()
  }

  block <- vroom::vroom(
    glue::glue("{pipeline_metadata_dir}updated_ld_blocks_to_colocalise.tsv"),
    show_col_types = F
  ) |>
    dplyr::filter(data_dir == ld_info$ld_block_data)

  finemapped_file <- glue::glue("{ld_info$ld_block_data}/finemapped_studies.tsv")

  gwas_upload_ids_to_compare <- character(0)
  raw_compare_arg <- args$gwas_upload_ids_to_compare
  if (!is.na(raw_compare_arg) && nchar(trimws(raw_compare_arg)) > 0) {
    gwas_upload_ids_to_compare <- strsplit(trimws(raw_compare_arg), "\\s*,\\s*")[[1]]
    gwas_upload_ids_to_compare <- gwas_upload_ids_to_compare[nchar(gwas_upload_ids_to_compare) > 0]
  }

  if (!is.na(args$worker_guid)) {
    finemapped_studies <- load_worker_finemapped_studies(
      args$ld_block,
      args$worker_guid,
      args$worker_p_value_threshold,
      gwas_upload_ids_to_compare
    )
  } else if (file.exists(finemapped_file)) {
    finemapped_studies <- vroom::vroom(finemapped_file, col_types = finemapped_column_types, show_col_types = F) |>
      dplyr::arrange(unique_study_id)
  } else {
    finemapped_studies <- data.frame()
  }

  if ((is.na(args$worker_guid) && !file.exists(finemapped_file)) || nrow(block) == 0 || nrow(finemapped_studies) == 0) {
    message(glue::glue("{args$ld_block}: Nothing to coloc, skipping."))
    vroom::vroom_write(data.frame(), args$completed_output_file)
    return()
  }

  block_list <- NULL
  if (!is.null(args$block_list) && !is.na(args$block_list) && file.exists(args$block_list)) {
    block_list <- vroom::vroom(args$block_list, show_col_types = FALSE)
  }

  if (!is.null(block_list) && nrow(block_list) > 0) {
    blocked <- is_study_blocked(block_list, finemapped_studies$study, finemapped_studies$cis_trans)
    if (sum(blocked) > 0) {
      finemapped_studies <- finemapped_studies |> dplyr::filter(!blocked)
      message(glue::glue("{args$ld_block}: Excluded {sum(blocked)} blocked studies from colocalisation"))
    }
  }

  finemapped_studies <- dplyr::arrange(finemapped_studies, unique_study_id)

  study_pairs <- get_study_pairs_to_coloc(
    finemapped_studies,
    coloc_results,
    args$worker_guid,
    gwas_upload_ids_to_compare
  )
  message(
    glue::glue("{args$ld_block}: Found {nrow(study_pairs)} study pairs to coloc in {diff_time_taken(start_time)}")
  )
  message(glue::glue("Existing pairs: {nrow(coloc_results)}"))

  if (nrow(study_pairs) == 0) {
    message(glue::glue("{args$ld_block}: No new study pairs to coloc."))
    if (nrow(coloc_results) > 0) {
      coloc_results <- coloc_results |> dplyr::mutate(ld_block = args$ld_block)
      vroom::vroom_write(coloc_results, coloc_results_file)
    }
    if (!is.na(args$worker_guid) && nrow(finemapped_studies) > 0) {
      vroom::vroom_write(resolve_worker_file_paths(finemapped_studies), finemapped_file)
    }
    vroom::vroom_write(data.frame(), args$completed_output_file)
    return()
  }

  if (!is.na(args$worker_guid)) {
    finemapped_subset <- finemapped_studies |>
      dplyr::filter(
        min_p <= lowest_p_value_threshold &
          (unique_study_id %in% study_pairs$unique_study_a | unique_study_id %in% study_pairs$unique_study_b) &
          !ignore
      )
    finemapped_subset <- resolve_worker_file_paths(finemapped_subset)
  } else {
    finemapped_subset <- finemapped_studies
  }

  studies_to_colocalise <- load_studies_to_colocalise(finemapped_subset)
  message(
    glue::glue("{args$ld_block}: Loaded {length(studies_to_colocalise)} studies in {diff_time_taken(start_time)}")
  )

  gc(verbose = FALSE)

  coloc_results <- run_coloc_for_study_pairs(
    study_pairs = study_pairs,
    studies_to_colocalise = studies_to_colocalise,
    coloc_results = coloc_results,
    worker_guid = args$worker_guid,
    ld_block = args$ld_block,
    start_time = start_time
  )

  coloc_results <- coloc_results |> dplyr::mutate(ld_block = args$ld_block)
  vroom::vroom_write(coloc_results, coloc_results_file)
  message(glue::glue("{args$ld_block}: Wrote {nrow(coloc_results)} pairwise coloc results"))

  if (!is.na(args$worker_guid)) {
    vroom::vroom_write(resolve_worker_file_paths(finemapped_studies), finemapped_file)
  }

  vroom::vroom_write(data.frame(), args$completed_output_file)
  return()
}

load_studies_to_colocalise <- function(finemapped_subset) {
  studies_to_colocalise <- lapply(seq_len(nrow(finemapped_subset)), function(i) {
    file <- finemapped_subset$file[i]

    if (!file.exists(file)) {
      message(glue::glue("{file} does not exist, skipping."))
      return(NULL)
    }

    gwas <- data.table::fread(
      file,
      sep = "\t",
      select = c("SNP", "LBF"),
      colClasses = c(SNP = "character", LBF = "numeric"),
      verbose = FALSE,
      showProgress = FALSE,
      nThread = 1
    )

    lbf_vector <- gwas$LBF
    names(lbf_vector) <- gwas$SNP
    rm(gwas)

    return(lbf_vector)
  })
  names(studies_to_colocalise) <- finemapped_subset$unique_study_id
  studies_to_colocalise <- studies_to_colocalise[!sapply(studies_to_colocalise, is.null)]
  return(studies_to_colocalise)
}

run_coloc_for_study_pairs <- function(study_pairs, studies_to_colocalise, coloc_results,
                                      worker_guid, ld_block, start_time) {
  results <- lapply(seq_len(nrow(study_pairs)), function(i) {
    pair <- study_pairs[i, ]
    first_gwas <- studies_to_colocalise[[pair$unique_study_a]]
    second_gwas <- studies_to_colocalise[[pair$unique_study_b]]

    tryCatch(
      {
        result <- pairwise_coloc_analysis(first_gwas, second_gwas)
      },
      error = function(e) {
        message(glue::glue("Error colocating {pair$unique_study_a} and {pair$unique_study_b}: {e}"))
        stop(glue::glue("Error colocating {pair$unique_study_a} and {pair$unique_study_b}: {e}"))
      }
    )
    if (is.null(result)) {
      result <- data.frame(
        unique_study_a = pair$unique_study_a,
        study_a = NA,
        unique_study_b = pair$unique_study_b,
        study_b = NA,
        bp_distance = NA,
        ignore = T,
        false_positive = F,
        false_negative = F,
        nsnps = NA,
        hit1 = NA,
        hit2 = NA,
        PP.H0.abf = NA,
        PP.H1.abf = NA,
        PP.H2.abf = NA,
        PP.H3.abf = NA,
        PP.H4.abf = NA,
        idx1 = NA,
        idx2 = NA,
        h4 = NA,
        ld_block = ld_block
      )
      return(result)
    }

    result <- dplyr::bind_cols(pair, result)
    return(result)
  })
  new_coloc_results <- dplyr::bind_rows(results[!sapply(results, is.null)])

  message(
    glue::glue("{ld_block}: Colocated {nrow(new_coloc_results)} study pairs in {diff_time_taken(start_time)}")
  )

  coloc_results <- dplyr::bind_rows(coloc_results, new_coloc_results)

  if (!is.na(worker_guid) && nrow(study_pairs) > 0) {
    main_pipeline_coloc_file <- glue::glue("{data_dir}/ld_blocks/{ld_block}/coloc_pairwise_results.tsv.gz")
    if (file.exists(main_pipeline_coloc_file)) {
      main_pipeline_coloc <- vroom::vroom(
        main_pipeline_coloc_file,
        delim = "\t",
        show_col_types = FALSE
      )

      studies_in_pairs <- unique(c(study_pairs$unique_study_a, study_pairs$unique_study_b))
      main_pipeline_coloc <- main_pipeline_coloc |>
        dplyr::filter(
          unique_study_a %in% studies_in_pairs & unique_study_b %in% studies_in_pairs
        )

      coloc_results <- dplyr::bind_rows(coloc_results, main_pipeline_coloc)
      message(glue::glue(
        "{ld_block}: Added {nrow(main_pipeline_coloc)} main pipeline coloc pairs for clustering"
      ))
    }
  }

  return(coloc_results)
}

get_study_pairs_to_coloc <- function(studies, existing_results, worker_guid, compare_guids) {
  studies <- dplyr::mutate(studies, id = dplyr::row_number()) |>
    dplyr::filter(min_p <= lowest_p_value_threshold | !ignore)
  studies <- data.table::as.data.table(studies)

  pairs_filtered <- studies[studies, on = .(id < id), allow.cartesian = TRUE][
    , bp_distance := abs(i.bp - bp)
  ][
    bp_distance <= bp_range | i.study == study
  ][
    , .(
      unique_study_a = unique_study_id,
      study_a = study,
      unique_study_b = i.unique_study_id,
      study_b = i.study,
      bp_distance = bp_distance,
      ignore = ignore,
      false_positive = F,
      false_negative = F
    )
  ] |>
    tibble::as_tibble()

  if (!is.na(worker_guid) && length(compare_guids) == 0) {
    pairs_filtered <- pairs_filtered |> dplyr::filter(study_a == worker_guid | study_b == worker_guid)
  } else if (!is.na(worker_guid) && length(compare_guids) > 0) {
    pairs_filtered <- pairs_filtered |>
      dplyr::filter(
        study_a == worker_guid | study_b == worker_guid | study_a %in% compare_guids | study_b %in% compare_guids
      )
  }

  existing_results <- data.table::as.data.table(existing_results)
  if (nrow(existing_results) > 0) {
    pairs_filtered <- dplyr::anti_join(pairs_filtered, existing_results, by = c("unique_study_a", "unique_study_b"))
  }

  return(pairs_filtered)
}

pairwise_coloc_analysis <- function(first_gwas, second_gwas) {
  harmonised_gwases <- harmonise_gwases(first_gwas, second_gwas)
  if (length(harmonised_gwases) == 0) {
    return(NULL)
  }

  first_lbf <- harmonised_gwases[[1]]
  second_lbf <- harmonised_gwases[[2]]
  nsnps <- length(first_lbf)

  if (nsnps < 50 || length(second_lbf) < 50) {
    return(NULL)
  }

  dynamic_priors <- calculate_dynamic_coloc_priors(nsnps)
  result <- coloc::coloc.bf_bf(
    bf1 = first_lbf,
    bf2 = second_lbf,
    p1 = coloc_base_p1,
    p2 = coloc_base_p2,
    p12 = coloc_base_p12 
  )
  result$summary$h4 <- result$summary$PP.H4.abf
  result$summary$dynamic_p1 <- coloc_base_p1
  result$summary$dynamic_p2 <- coloc_base_p2
  result$summary$dynamic_p12 <- coloc_base_p12
  result$summary$nsnps <- nsnps
  coloc_results <- result$summary
  return(coloc_results)
}

harmonise_gwases <- function(...) {
  gwases <- list(...)

  snpids <- Reduce(intersect, lapply(gwases, function(gwas) names(gwas)))
  if (length(snpids) <= 1) {
    return(list())
  }
  snpids <- sort(snpids)

  gwases <- lapply(gwases, function(gwas) {
    gwas <- gwas[!duplicated(names(gwas))]
    gwas <- gwas[names(gwas) %in% snpids]
    gwas <- gwas[match(snpids, names(gwas))]
    gwas <- gwas[!is.na(gwas)]
    return(gwas)
  })

  final_snpids <- Reduce(intersect, lapply(gwases, function(gwas) names(gwas)))
  if (length(final_snpids) <= 1) {
    return(list())
  }
  final_snpids <- sort(final_snpids)

  gwases <- lapply(gwases, function(gwas) {
    gwas <- gwas[match(final_snpids, names(gwas))]
    return(gwas)
  })

  stopifnot(identical(names(gwases[[1]]), names(gwases[[2]])))
  return(gwases)
}

calculate_dynamic_coloc_priors <- function(nsnps) {
  nsnps <- max(1, as.numeric(nsnps))
  target_h1 <- coloc_base_p1 * wallace_reference_nsnps
  target_h2 <- coloc_base_p2 * wallace_reference_nsnps
  target_h4 <- coloc_base_p12 * wallace_reference_nsnps

  max_per_snp <- 0.99 / nsnps
  p1 <- min(target_h1 / nsnps, max_per_snp)
  p2 <- min(target_h2 / nsnps, max_per_snp)
  p12 <- min(target_h4 / nsnps, p1, p2, 0.99 / nsnps)
  p12 <- max(p12, .Machine$double.eps)

  # return(list(p1 = p1, p2 = p2, p12 = p12))
  return(list(p1 = coloc_base_p1, p2 = coloc_base_p2, p12 = coloc_base_p12))
}

resolve_worker_file_paths <- function(finemapped_studies) {
  finemapped_studies |>
    dplyr::mutate(
      file = dplyr::case_when(
        grepl("^gwas_upload", file) ~ glue::glue("{data_dir}/{file}"),
        grepl("^study", file) ~ glue::glue("{data_dir}/{file}"),
        TRUE ~ file
      ),
      file_with_lbfs = dplyr::case_when(
        is.na(file_with_lbfs) ~ NA_character_,
        grepl("^gwas_upload", file_with_lbfs) ~ glue::glue("{data_dir}/{file_with_lbfs}"),
        grepl("^study", file_with_lbfs) ~ glue::glue("{data_dir}/{file_with_lbfs}"),
        TRUE ~ file_with_lbfs
      ),
      svg_file = dplyr::case_when(
        is.na(svg_file) ~ NA_character_,
        grepl("^gwas_upload", svg_file) ~ glue::glue("{data_dir}/{svg_file}"),
        grepl("^study", svg_file) ~ glue::glue("{data_dir}/{svg_file}"),
        TRUE ~ svg_file
      )
    )
}

load_worker_finemapped_studies <- function(ld_block, worker_guid, worker_p_value_threshold, compare_guids) {
  finemapped_file <- glue::glue("{ld_block_dirs(ld_block)$ld_block_data}/finemapped_studies.tsv")
  if (!file.exists(finemapped_file)) {
    return(data.frame())
  }

  finemapped_studies <- vroom::vroom(finemapped_file, col_types = finemapped_column_types, show_col_types = F) |>
    dplyr::filter(study == worker_guid) |>
    dplyr::filter(min_p <= worker_p_value_threshold)
  worker_guid_bp <- finemapped_studies$bp[finemapped_studies$study == worker_guid]

  existing_finemapped_studies_file <- glue::glue("{data_dir}/ld_blocks/{ld_block}/finemapped_studies.tsv")
  if (file.exists(existing_finemapped_studies_file)) {
    existing_finemapped_studies <- vroom::vroom(
      existing_finemapped_studies_file,
      col_types = finemapped_column_types,
      show_col_types = F
    ) |>
      dplyr::filter(min_p <= min_p_allowed_for_worker & abs(min(worker_guid_bp) - bp) <= bp_range)

    finemapped_studies <- dplyr::bind_rows(finemapped_studies, existing_finemapped_studies) |>
      dplyr::arrange(unique_study_id)
  }

  for (compare_guid in compare_guids) {
    if (compare_guid == worker_guid) next
    compare_finemapped_file <- glue::glue(
      "{gwas_upload_dir}/ld_blocks/gwas_upload/{compare_guid}/{ld_block}/finemapped_studies.tsv"
    )
    if (!file.exists(compare_finemapped_file)) {
      message(glue::glue("Compare file {compare_finemapped_file} does not exist"))
      next
    }

    compare_finemapped <- vroom::vroom(
      compare_finemapped_file,
      col_types = finemapped_column_types,
      show_col_types = F
    ) |>
      dplyr::filter(study == compare_guid & min_p <= p_value_threshold)

    message(glue::glue("Found {nrow(compare_finemapped)} finemapped studies in {compare_finemapped_file}"))
    if (nrow(compare_finemapped) > 0) {
      finemapped_studies <- dplyr::bind_rows(finemapped_studies, compare_finemapped) |>
        dplyr::distinct(unique_study_id, .keep_all = TRUE) |>
        dplyr::arrange(unique_study_id)
      message(glue::glue(
        "{ld_block}: Added {nrow(compare_finemapped)} finemapped studies from compare upload {compare_guid}"
      ))
    }
  }

  return(finemapped_studies)
}

invisible(main())
