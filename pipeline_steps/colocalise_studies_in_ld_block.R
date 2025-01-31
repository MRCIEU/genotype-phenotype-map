source('constants.R')
bp_range <- 10000

parser <- argparser::arg_parser('Colocalise studies per region')
parser <- argparser::add_argument(parser, '--ld_block', help = 'LD block that the ', type = 'character')
parser <- argparser::add_argument(parser, '--completed_output_file', help = 'Coloc result file to save', type = 'character')
args <- argparser::parse_args(parser)

main <- function() {
  ld_info <- ld_block_dirs(args$ld_block)
  block <- vroom::vroom(glue::glue('{pipeline_metadata_dir}updated_ld_blocks_to_colocalise.tsv'), show_col_types=F) |>
    dplyr::filter(data_dir == ld_info$ld_block_data)

  finemapped_file <- glue::glue('{ld_info$ld_block_data}/finemapped_studies.tsv')
  if (file.exists(finemapped_file)) {
    finemapped_studies <- vroom::vroom(finemapped_file, delim = '\t', show_col_types = F) |>
      dplyr::filter(variant_type == variant_types$common & min_p <= p_value_threshold)
  }

  if (!file.exists(finemapped_file) || nrow(block) == 0 || nrow(finemapped_studies) == 0) {
    message(glue::glue('Nothing to coloc in LD region {ld_info$ld_block_data}, skipping.'))
    vroom::vroom_write(data.frame(), args$completed_output_file)
    return()
  }

  studies_to_colocalise <- lapply(finemapped_studies$file, function(file) vroom::vroom(file, show_col_types = F))
  names(studies_to_colocalise) <- finemapped_studies$unique_study_id

  grouped_studies <- group_studies_in_same_bp_range(finemapped_studies)
  if (length(grouped_studies) == 0) {
    message(glue::glue('No grouped studies for LD region {ld_info$ld_block_data}, skipping.'))
    vroom::vroom_write(data.frame(), args$completed_output_file)
    return()
  }

  hyprcoloc_results <- colocalise_based_on_group(studies_to_colocalise, grouped_studies, finemapped_studies)
  hyprcoloc_results <- post_coloc_filtering(hyprcoloc_results)

  coloc_results_file <- glue::glue('{ld_info$ld_block_data}/coloc_results.tsv')
  vroom::vroom_write(hyprcoloc_results, coloc_results_file)
  vroom::vroom_write(data.frame(), args$completed_output_file)
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
group_studies_in_same_bp_range <- function(studies) {
  grouped_studies <- list()
  for(i in seq_len(nrow(studies))) {
    study_bp <- studies[i,]$bp
    study <- studies[i,]$unique_study_id

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

colocalise_based_on_group <- function(studies, groupings, metadata) {
  results <- lapply(groupings, function(group) {
    specific_group <- studies[group]
    specific_group <- do.call(harmonise_gwases, specific_group)

    if (length(specific_group) == 0 || nrow(specific_group[[1]])==0) return()

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
      return(results)
  })
  return(results)
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
