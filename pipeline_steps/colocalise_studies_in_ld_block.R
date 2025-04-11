source('constants.R')
bp_range <- 10000

parser <- argparser::arg_parser('Colocalise studies per region')
parser <- argparser::add_argument(parser, '--ld_block', help = 'LD block that the ', type = 'character')
parser <- argparser::add_argument(parser, '--completed_output_file', help = 'Coloc result file to save', type = 'character')
parser <- argparser::add_argument(parser, '--worker_guid', help = 'Worker GUID', type = 'character', default = NA)
args <- argparser::parse_args(parser)

main <- function() {
  if (!is.na(args$worker_guid)) {
    update_directories_for_worker(args$worker_guid)
  }
  ld_info <- ld_block_dirs(args$ld_block)
  block <- vroom::vroom(glue::glue('{pipeline_metadata_dir}updated_ld_blocks_to_colocalise.tsv'), show_col_types=F) |>
    dplyr::filter(data_dir == ld_info$ld_block_data)

  finemapped_file <- glue::glue('{ld_info$ld_block_data}/finemapped_studies.tsv')
  if (file.exists(finemapped_file)) {
    finemapped_studies <- vroom::vroom(finemapped_file, col_types = finemapped_column_types, show_col_types = F) |>
      dplyr::filter(variant_type == variant_types$common & min_p <= p_value_threshold)
  }

  if (!is.na(args$worker_guid)) {
    existing_finemapped_studies_file <- glue::glue('{data_dir}/ld_blocks/{args$ld_block}/finemapped_studies.tsv')
    existing_finemapped_studies <- vroom::vroom(existing_finemapped_studies_file, col_types = finemapped_column_types, show_col_types=F)
    finemapped_studies <- dplyr::bind_rows(finemapped_studies, existing_finemapped_studies) |>
      dplyr::filter(variant_type == variant_types$common & min_p <= p_value_threshold) |>
      dplyr::mutate(file = sub('/local-scratch/projects/genotype-phenotype-map/data/', data_dir, file))
  }

  cached_coloc_group_file <- glue::glue('{ld_info$ld_block_data}/coloc_cached_groups.rds')
  cached_studies_file <- glue::glue('{ld_info$ld_block_data}/coloc_cached_studies.tsv')
  cached_data <- load_cached_data(finemapped_studies, cached_coloc_group_file, cached_studies_file)
  if (cached_data$skip_coloc_block) {
    message('No new common studies to compare, skipping.')
    vroom::vroom_write(data.frame(), args$completed_output_file)
    return()
  }

  if (!file.exists(finemapped_file) || nrow(block) == 0 || nrow(finemapped_studies) == 0) {
    message(glue::glue('Nothing to coloc in LD region {ld_info$ld_block_data}, skipping.'))
    vroom::vroom_write(data.frame(), args$completed_output_file)
    return()
  }

  studies_to_colocalise <- lapply(finemapped_studies$file, function(file) vroom::vroom(file, delim = '\t', show_col_types = F))
  names(studies_to_colocalise) <- finemapped_studies$unique_study_id

  start_time <- Sys.time()
  grouped_studies <- group_studies_in_same_bp_range(finemapped_studies, args$worker_guid)
  if (length(grouped_studies) == 0) {
    message(glue::glue('No grouped studies for LD region {ld_info$ld_block_data}, skipping.'))
    vroom::vroom_write(data.frame(), args$completed_output_file)
    return()
  }
  message(glue::glue('Grouping studies took {as.character(hms::as_hms(difftime(Sys.time(), start_time)))}'))

  start_time <- Sys.time()
  hyprcoloc_results <- colocalise_based_on_group(studies_to_colocalise, grouped_studies, finemapped_studies, cached_data$cached_coloc_groups)
  message(glue::glue('Colocalisation took {as.character(hms::as_hms(difftime(Sys.time(), start_time)))}'))

  start_time <- Sys.time()
  filtered_hyprcoloc_results <- post_coloc_filtering(hyprcoloc_results)
  message(glue::glue('Post-colocalisation filtering took {as.character(hms::as_hms(difftime(Sys.time(), start_time)))}'))

  coloc_results_file <- glue::glue('{ld_info$ld_block_data}/coloc_results.tsv')
  vroom::vroom_write(filtered_hyprcoloc_results, coloc_results_file)

  if (is.na(args$worker_guid)) {
    studies_to_cache <- data.frame(study=finemapped_studies$study)
    vroom::vroom_write(studies_to_cache, cached_studies_file)
    saveRDS(hyprcoloc_results, cached_coloc_group_file)
  }

  vroom::vroom_write(data.frame(), args$completed_output_file)
}

load_cached_data <- function(finemapped_studies, cached_coloc_group_file, cached_studies_file) {
  skip_coloc_block <- FALSE
  if (file.exists(cached_coloc_group_file)) {
    cached_coloc_groups <- readRDS(cached_coloc_group_file)
  } else {
    cached_coloc_groups <- list()
  }

  if (file.exists(cached_studies_file)) {
    cached_studies <- vroom::vroom(cached_studies_file, delim = '\t', show_col_types = F)
    already_run <- setdiff(finemapped_studies$study, cached_studies$study)
    if (length(already_run) == 0) {
      skip_coloc_block <- TRUE
    }
  }

  return(list(cached_coloc_groups = cached_coloc_groups,
              skip_coloc_block = skip_coloc_block
             )
  )
}

post_coloc_filtering <- function(hyprcoloc_results) {
  if (is.null(hyprcoloc_results)) hyprcoloc_results <- data.frame()
  hyprcoloc_results <- lapply(hyprcoloc_results, function(result) {
    if (!is.null(result)) return(result$results)
  }) |> dplyr::bind_rows()

  hyprcoloc_results <- hyprcoloc_results[!duplicated(hyprcoloc_results$traits), ]
  traits <- strsplit(hyprcoloc_results$traits, ', ')

  #we have to remove subsets again because hyprcoloc breaks the larger groups up and reports on different iterations
  subsets <- find_subsets_of_grouped_studies(traits)
  if (length(subsets) > 0) hyprcoloc_results <- hyprcoloc_results[-subsets, ]

  return(hyprcoloc_results)
}

#' group_studies_in_same_bp_range, in 2 steps
#'   1. created 'grouped_studies' list with all studies whose top snp is within bp_range, and only one study per group
#'   2. Filter out all grouped studies that are of length 1, or that are subsets of other groups
group_studies_in_same_bp_range <- function(studies, worker_guid) {
  grouped_studies <- list()
  for(i in seq_len(nrow(studies))) {
    study <- studies[i,]$unique_study_id
    study_name <- studies[i,]$study
    if (!is.na(worker_guid) && study_name != worker_guid) next

    study_bp <- studies[i,]$bp
    filtered_studies <- dplyr::filter(studies, (bp-bp_range) < study_bp & study_bp < (bp+bp_range))
    filtered_studies <- filtered_studies[!duplicated(filtered_studies$study), ]

    found_grouped_studies <- filtered_studies$unique_study_id

    if (length(found_grouped_studies) > 1) {
      grouped_studies[[study]] <- found_grouped_studies
    }
  }

  len <- sapply(grouped_studies, length)
  grouped_studies <- grouped_studies[order(len)]

  if (length(grouped_studies) == 0) {
    return(list())
  }

  subsets <- find_subsets_of_grouped_studies(grouped_studies)
  if (length(subsets) > 0) grouped_studies <- grouped_studies[-subsets]
  return(grouped_studies)
}

find_subsets_of_grouped_studies <- function(grouped_studies) {
  subsets <- c()
  if (length(grouped_studies) < 2) return(subsets)
  for (i in 1:(length(grouped_studies)-1)) {
    for (j in (i+1):length(grouped_studies)) {
      if (all(grouped_studies[[i]] %in% grouped_studies[[j]])) {
        subsets <- c(subsets, i)
      }
    }
  }

  subsets <- subsets[!duplicated(subsets)]
  return(subsets)
}

colocalise_based_on_group <- function(studies, groupings, metadata, cached_coloc_groups) {
  all_results <- lapply(groupings, function(group) {
    specific_group <- studies[group]

    unique_group_id <- paste(unlist(group), collapse = '')
    cached_result <- find_cached_coloc_groups(cached_coloc_groups, unique_group_id)

    if (!is.null(cached_result)) {
      return(cached_result)
    } else {
      specific_group <- do.call(harmonise_gwases, specific_group)
      if (length(specific_group) == 0 || nrow(specific_group[[1]]) == 0) return()

      snps <- specific_group[[1]]$SNP
      trait_names <- names(specific_group)
      categories <- dplyr::filter(metadata, study %in% trait_names)$category
      binary_outcomes <- lapply(categories, function(category) {
        return (category == study_categories$binary)
      })

      beta_matrix <- lapply(specific_group, function(study) study$BETA) |>
        dplyr::bind_cols() |>
        as.matrix()
      se_matrix <- lapply(specific_group, function(study) study$SE) |>
        dplyr::bind_cols() |>
        as.matrix()

      results <- hyprcoloc::hyprcoloc(effect.est = beta_matrix,
                                      effect.se = se_matrix,
                                      trait.names = trait_names,
                                      binary.outcomes = binary_outcomes,
                                      snp.id = snps,
                                      snpscores = T
      )
      results$unique_group_id <- unique_group_id
      return(results)
    }
  })
  return(all_results)
}


find_cached_coloc_groups <- function(cached_coloc_groups, unique_group_id) {
  cached_result <- Filter(function(cached_result) cached_result$unique_group_id == unique_group_id, cached_coloc_groups)

  if (length(cached_result) == 0) return(NULL)
  else {
    return(cached_result[[1]])
  }
}


harmonise_gwases <- function(...) {
  gwases <- list(...)
  snpids <- Reduce(intersect, lapply(gwases, function(gwas) gwas$SNP))
  if (length(snpids) <= 1) return(list())

  gwases <- lapply(gwases, function(gwas) {
    dplyr::filter(gwas, SNP %in% snpids & !duplicated(SNP)) |>
      dplyr::arrange(SNP)
  })

  return(gwases)
}

main()
