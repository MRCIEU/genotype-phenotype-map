source('constants.R')
library(argparser, quietly = TRUE)
bp_range <- 10000

parser <- argparser::arg_parser('Colocalise studies per region')
parser <- argparser::add_argument(parser, '--ld_block', help = 'LD block that the ', type = 'character')
parser <- argparser::add_argument(parser, '--coloc_result_file', help = 'Coloc result file to save', type = 'character')
args <- argparser::parse_args(parser)

main <- function(args) {
  ld_info <- ld_block_dirs(args$ld_block)
  block <- vroom::vroom(paste0(pipeline_metadata_dir, 'updated_ld_blocks_to_colocalise.tsv'), show_col_types=F) |>
    dplyr::filter(data_dir == ld_info$ld_block_data)

  finemapped_file <- paste0(ld_info$ld_block_data, '/finemapped_studies.tsv')
  if (file.exists(finemapped_file)) {
    finemapped_studies <- vroom::vroom(finemapped_file, show_col_types = F)
    finemapped_studies$unique_study_id <- paste0(finemapped_studies$study, "_", file_prefix(finemapped_studies$file))
  }

  if (!file.exists(finemapped_file) || nrow(block) == 0 || nrow(finemapped_studies) == 0) {
    coloc_result_dir <- dirname(args$coloc_result_file)
    coloc_files <- Sys.glob(paste0(coloc_result_dir, '/hyprcoloc_results*'))
    coloc_files <- sort(coloc_files, decreasing=T)
    if (length(coloc_files) == 0) {
      vroom::vroom_write(data.frame(), args$coloc_result_file)
    } else {
      file.symlink(coloc_files[1], args$coloc_result_file)
    }

    message(paste0('Nothing to process for LD region ', ld_info$ld_block_data ,', skipping.'))
    return()
  }

  studies_to_colocalise <- lapply(finemapped_studies$file, function(file) {
    return(vroom::vroom(file, show_col_types = F))
  })
  names(studies_to_colocalise) <- finemapped_studies$unique_study_id

  grouped_studies <- group_studies_in_same_bp_range(finemapped_studies)
  if (length(grouped_studies) == 0) {
    vroom::vroom_write(data.frame(), args$coloc_result_file)
    return()
  }

  hyprcoloc_results <- colocalise_based_on_group(studies_to_colocalise, grouped_studies, finemapped_studies)

  all_results <- lapply(hyprcoloc_results, function(result) {
    if (!is.null(result)) return(result$results)
  }) |> dplyr::bind_rows()

  if (is.null(all_results)) all_results <- data.frame()

  vroom::vroom_write(all_results, args$coloc_result_file)
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

  to_remove <- c()
  for (i in 1:(length(grouped_studies)-1)) {
    for (j in (i+1):length(grouped_studies)) {
      if (all(grouped_studies[[i]] %in% grouped_studies[[j]])) {
        to_remove <- c(to_remove, i)
      }
    }
  }
  grouped_studies <- grouped_studies[-to_remove]

  return(grouped_studies)
}

colocalise_based_on_group <- function(studies, groupings, metadata) {
  results <- lapply(groupings, function(group) {
    specific_group <- studies[group]
    specific_group <- do.call(harmonise_gwases, specific_group)
    if (nrow(specific_group[[1]])==0) return()

    snps <- specific_group[[1]]$RSID
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
  snpids <- Reduce(intersect, lapply(gwases, function(gwas) gwas$RSID))
  message(paste('Number of shared SNPs after harmonisation:', length(snpids)))

  gwases <- lapply(gwases, function(gwas) {
    dplyr::filter(gwas, RSID %in% snpids & !duplicated(RSID)) |>
      dplyr::arrange(RSID)
  })

  return(gwases)
}

main(args)
