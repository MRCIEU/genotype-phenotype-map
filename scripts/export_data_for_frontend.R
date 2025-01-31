source('../pipeline_steps/constants.R')

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
print(length(associated_studies))

only_interesting_blocks <- dplyr::filter(all_study_blocks, unique_study_id %in% associated_studies) 
only_interesting_blocks <- only_interesting_blocks |>
  dplyr::select(trait, study, unique_study_id, chr, bp, min_p, cis_trans, gene, tissue, data_type, sample_size)

coloc_results_for_study <- dplyr::filter(coloc_results_for_study, unique_study_a %in% only_interesting_blocks$unique_study_id & unique_study_b %in% only_interesting_blocks$unique_study_id)
print(length(coloc_results_for_study))

freq_of_snp <- table(coloc_results_for_study$candidate_snp)
coloc_results_for_study$num_unique_studies <- as.numeric(freq_of_snp[coloc_results_for_study$candidate_snp])

combined_data <- list(
  studies = only_interesting_blocks,
  colocs = coloc_results_for_study
)

coloc_result_json <- jsonlite::toJSON(combined_data, pretty = T)
write(coloc_result_json, '~/coloc_result.json')
