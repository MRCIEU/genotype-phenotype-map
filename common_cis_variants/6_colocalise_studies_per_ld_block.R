source('constants.R')
library(argparser, quietly = TRUE)
bp_range <- 10000

parser <- arg_parser("Colocalising ld blocks")
parser <- add_argument(parser,
                       "--chr",
                       help = "Chromosome",
                       type = "numeric"
)
args <- parse_args(parser)

ld_regions <- vroom::vroom("data/ld_regions.tsv")
ld_blocks_to_colocalise <- vroom::vroom(paste0(pipeline_metadata_dir, "updated_ld_blocks_to_colocalise.tsv")) |>
  dplyr::filter(chr == args$chr)

harmonise_gwases <- function(...) {
  gwases <- list(...)

  snpids <- Reduce(intersect, lapply(gwases, function(gwas) gwas$rsid))
  message(paste("Number of shared SNPs after harmonisation:", length(snpids)))

  gwases <- lapply(gwases, function(gwas) {
    dplyr::filter(gwas, rsid %in% snpids & !duplicated(rsid)) |>
      dplyr::arrange(rsid)
  })

  return(gwases)
}

harmonise_gwases_missing_rows <- function(gwases=list()) {
  all_snps <- unique(do.call(c, lapply(gwases, function(gwas) gwas$rsid)))
  message(paste("Number of SNPs across all gwases:", length(all_snps)))

  gwases <- lapply(gwases, function(gwas) {
    snp_diff <- setdiff(all_snps, gwas$rsid)
    message(paste('adding', length(snp_diff), 'rows to gwas'))
    tidyr::complete(gwas, rsid=snp_diff, beta=0, se=1, lp=0.01)
  })

  return(gwases)
}

all_results <- apply(ld_blocks_to_colocalise, 1, function(block) {
  specific_ld_block_data_dir <- block[['data_dir']]
  ld_block_results_dir <- block[['results_dir']]
  extracted_studies <- vroom::vroom(paste0(specific_ld_block_data_dir, "/finemapped_results.tsv"), show_col_types = F)

  studies_to_colocalise <- lapply(extracted_studies$file, function(file) vroom::vroom(file, show_col_types = F))
  names(studies_to_colocalise) <- file_prefix(extracted_studies$study)

  grouped_studies <- list()
  for(i in seq_len(nrow(extracted_studies))) {
    study_bp <- extracted_studies[i,]$bp
    study <- extracted_studies[i,]$unique_study_id
    found_grouped_studies <- dplyr::filter(extracted_studies, (bp-bp_range) < study_bp & study_bp < (bp+bp_range))$unique_study_id

    #TODO: change this to maybe ad all at the beginning, then sort by length and remove subsets
    if (!list(found_grouped_studies) %in% grouped_studies) {
      grouped_studies[study] <- list(found_grouped_studies)
    }
  }
  grouped_studies <- Filter(function(s) length(s) > 1, grouped_studies)

  results <- lapply(grouped_studies, function(group) {
    specific_group <- studies_to_colocalise[group]
    specific_group <- do.call(harmonise_gwases, specific_group)
    if (nrow(specific_group[[1]])==0) return()

    snps <- specific_group[[1]]$rsid
    trait_names <- names(specific_group)
    beta_matrix <- lapply(specific_group, function(study) study$beta) |>
      dplyr::bind_cols() |>
      as.matrix()
    se_matrix <- lapply(specific_group, function(study) study$se) |>
      dplyr::bind_cols() |>
      as.matrix()
    binary_outcomes <- lapply(specific_group, function(group) {
      return(0)
      #return (category == data_categories$binary)
    })

    tryCatch(
      expr = {
        results <- hyprcoloc::hyprcoloc(beta_matrix, se_matrix,
                                        trait.names = trait_names,
                                        binary.outcomes = binary_outcomes,
                                        snp.id = snps,
                                        snpscores = T)

        return(results)
      },
      error = function(e) {
        print(e)
        return()
      }
    )
  })
  saveRDS(results, paste0(ld_block_results_dir, "/hyprcoloc_results.rds"))

  significant_results <- lapply(results, function(result) {
    if (is.null(result)) return()
    dplyr::filter(result$results, posterior_prob >= 0.8)
  })
  significant_results <- Filter(function(result) !is.null(result) && nrow(result) != 0, significant_results) |>
    dplyr::bind_rows()
  vroom::vroom_write(significant_results, paste0(ld_block_results_dir, "/significant_results.tsv"))
  print(paste0("LD Block Results: ", ld_block_results_dir, "/hyprcoloc_results.tsv"))

  return(significant_results)
})
