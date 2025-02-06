source('../pipeline_steps/constants.R')

# TODO: change this to include candidate snp, posterior prob and pp explained by snp for each study


additional_coloc_columns <- function(studies, coloc_results) {
  study_a_match <- match(coloc_results$unique_study_a, studies$unique_study_id)
  study_b_match <- match(coloc_results$unique_study_b, studies$unique_study_id)

  coloc_results$trait_a <- studies$trait[study_a_match]
  coloc_results$trait_b <- studies$trait[study_b_match]
  coloc_results$data_type_a <- studies$data_type[study_a_match]
  coloc_results$data_type_b <- studies$data_type[study_b_match]
  coloc_results$cis_trans_a <- studies$cis_trans[study_a_match]
  coloc_results$cis_trans_b <- studies$cis_trans[study_b_match]
  coloc_results$min_p_a <- studies$min_p[study_a_match]
  coloc_results$min_p_b <- studies$min_p[study_b_match]
  coloc_results$gene_a <- studies$known_gene[match(coloc_results$unique_study_a, studies$unique_study_id)]
  coloc_results$gene_b <- studies$known_gene[match(coloc_results$unique_study_b, studies$unique_study_id)]

  return(coloc_results)
}

results_dir <- '/local-scratch/projects/genotype-phenotype-map/results/'
last_pipeline_run <- head(sort(list.files(results_dir, pattern = "^20", include.dirs=T), decreasing=T), 1)[1]
last_result_dir <- glue::glue('{results_dir}/{last_pipeline_run}')

studies_processed <- vroom::vroom(glue::glue('{last_result_dir}/studies_processed.tsv'), show_col_types=F)
coloc_results <- vroom::vroom(glue::glue('{last_result_dir}/coloc_results.tsv'), show_col_types=F) |>
  dplyr::filter(posterior_prob > 0.5) |>
  tidyr::separate(candidate_snp, into=c("chr", "bp"), sep="[:_]", remove=F) |>
  dplyr::mutate(chr = as.numeric(chr), bp = as.numeric(bp))
all_study_blocks <- vroom::vroom(glue::glue('{last_result_dir}/all_study_blocks.tsv'), show_col_types=F)
all_study_blocks <- merge(all_study_blocks, studies_processed, by.x='study', by.y='study_name')

study_of_interest <- 'ebi-a-GCST90029025'
coloc_results_for_study <- dplyr::filter(coloc_results, grepl(study_of_interest, unique_study_a) | grepl(study_of_interest, unique_study_b))
associated_studies <- unique(c(coloc_results_for_study$unique_study_a, coloc_results_for_study$unique_study_b))

only_interesting_blocks <- dplyr::filter(all_study_blocks, unique_study_id %in% associated_studies) 
only_interesting_blocks <- only_interesting_blocks |>
  dplyr::select(trait, study, unique_study_id, chr, bp, min_p, cis_trans, gene, tissue, data_type, sample_size)

snp_a <- coloc_results_for_study$candidate_snp[match(only_interesting_blocks$unique_study_id, coloc_results_for_study$unique_study_a)]
snp_b <- coloc_results_for_study$candidate_snp[match(only_interesting_blocks$unique_study_id, coloc_results_for_study$unique_study_b)]

only_interesting_blocks$candidate_snp <- pmax(snp_a, snp_b, na.rm = T)

only_interesting_blocks_by_snp <- split(only_interesting_blocks, only_interesting_blocks$candidate_snp)
only_interesting_blocks_by_snp <- lapply(only_interesting_blocks_by_snp, function(group) {
  return(list(
    snp = as.character(group$candidate_snp[[1]]),
    posterior_prob = as.numeric(coloc_results_for_study$posterior_prob[coloc_results_for_study$candidate_snp == group$candidate_snp[1]][[1]]),
    posterior_explained_by_snp  = as.numeric(coloc_results_for_study$posterior_explained_by_snp[coloc_results_for_study$candidate_snp == group$candidate_snp[1]][[1]]),
    studies = group |>
      dplyr::select(-chr, -bp)
  ))
})
jsonlite::toJSON(only_interesting_blocks_by_snp[[1]], pretty = T, auto_unbox = T)

coloc_results_for_study <- coloc_results_for_study |>
  dplyr::filter(unique_study_a %in% only_interesting_blocks$unique_study_id & unique_study_b %in% only_interesting_blocks$unique_study_id)

coloc_results_for_study <- additional_coloc_columns(only_interesting_blocks, coloc_results_for_study)
coloc_results_for_study <- coloc_results_for_study |>
  dplyr::mutate(min_p = pmin(min_p_a, min_p_b),
               includes_trans = (!is.na(cis_trans_a) & cis_trans_a == 'trans') | (!is.na(cis_trans_b) & cis_trans_b == 'trans'),
               includes_qtl = data_type_a != 'phenotype' | data_type_b != 'phenotype') |>
  dplyr::select(-min_p_a, -min_p_b, -cis_trans_a, -cis_trans_b, -data_type_a, -data_type_b, -posterior_explained_by_snp, -unique_study_a, -unique_study_b)


freq_of_snp <- table(coloc_results_for_study$candidate_snp)
coloc_results_for_study$num_unique_studies <- as.numeric(freq_of_snp[coloc_results_for_study$candidate_snp])

combined_data <- list(
  studies = only_interesting_blocks_by_snp,
  colocs = coloc_results_for_study
)

coloc_result_json <- jsonlite::toJSON(combined_data, pretty = T, auto_unbox = T)
write(coloc_result_json, '~/coloc_result.json')


# Getting data for SNP

# specific_snp <- "21:46628861_A_G" #this one has 24 entries
specific_snp <- "2:210675783_A_C" #this one has 168 entries
coloc_for_snp <- coloc_results[coloc_results$candidate_snp == specific_snp,]
associated_studies <- unique(c(coloc_for_snp$unique_study_a, coloc_for_snp$unique_study_b))
blocks_for_snp <- dplyr::filter(all_study_blocks, unique_study_id %in% associated_studies) |>
  dplyr::select(trait, study, unique_study_id, chr, bp, min_p, cis_trans, gene, tissue, data_type, sample_size)

variant_annotations <- vroom::vroom(glue::glue('{variant_annotation_dir}/vep_annotations_hg38.tsv.gz'), show_col_types = F) |>
  dplyr::filter(SNP == specific_snp)

  snp_data_json <- list(
    annotation = variant_annotations[1,],
    studies = blocks_for_snp
  )
snp_result_json <- jsonlite::toJSON(snp_data_json, pretty = T, auto_unbox = T)
write(snp_result_json, '~/snp_result.json')