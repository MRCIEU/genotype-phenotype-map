source('../pipeline_steps/constants.R')

results_dir <- '/local-scratch/projects/genotype-phenotype-map/results/'
last_pipeline_run <- head(sort(list.files(results_dir, pattern = "^20", include.dirs=T), decreasing=T), 1)[1]
last_result_dir <- glue::glue('{results_dir}/{last_pipeline_run}')

studies_processed <- vroom::vroom(glue::glue('{last_result_dir}/studies_processed.tsv'), show_col_types=F)
coloc_results <- vroom::vroom(glue::glue('{last_result_dir}/coloc_results.tsv'), show_col_types=F)
all_study_blocks <- vroom::vroom(glue::glue('{last_result_dir}/all_study_blocks.tsv'), show_col_types=F)
raw_coloc_results <- vroom::vroom(glue::glue('{last_result_dir}/raw_coloc_results.tsv'), show_col_types=F) |>
  dplyr::filter(posterior_prob > 0.5) |>
  dplyr::rename(unique_studies = traits) |> 
  dplyr::select(-iteration) |>
  tidyr::separate(candidate_snp, into=c("chr", "bp"), sep="[:_]", remove=F) |>
  dplyr::mutate(chr = as.numeric(chr), bp = as.numeric(bp))

raw_coloc_results$unique_studies <- strsplit(raw_coloc_results$unique_studies, ", ")
# filtered_studies <- sapply(raw_coloc_results$unique_studies, function(studies) {
#   for (study in studies) {
#     if (study %in% all_study_blocks$unique_study_id) return(study)
#     else return()
#   }
# })
# print(head(filtered_studies))
# q()

raw_coloc_results$studies <- lapply(raw_coloc_results$unique_studies, function(studies) {
  return(sapply(studies, function(study) sub('_.*', '', study)))
})

all_study_blocks[all_study_blocks$unique_study_id == 'UKB-PPP-european-CTSE:P14091:OID30628:v1_EUR_6_25790150_1', ]
all_study_blocks <- merge(all_study_blocks, studies_processed, by.x='study', by.y='study_name')
all_study_blocks[all_study_blocks$unique_study_id == 'UKB-PPP-european-CTSE:P14091:OID30628:v1_EUR_6_25790150_1', ]

study_of_interest <- 'ebi-a-GCST90029025'
# all_coloc_results_for_study <- dplyr::filter(coloc_results, grepl(study_of_interest, unique_study_a) | grepl(study_of_interest, unique_study_b))

is_in <- sapply(raw_coloc_results$studies, function(studies) study_of_interest %in% studies)
raw_coloc_results <- raw_coloc_results[is_in, ]
associated_studies <- unique(unlist(raw_coloc_results$unique_studies))

only_interesting_blocks <- dplyr::filter(all_study_blocks, unique_study_id %in% associated_studies) 
only_interesting_blocks <- only_interesting_blocks |>
  dplyr::select(trait, study, unique_study_id, chr, bp, min_p, cis_trans, gene, tissue, data_type, sample_size)

combined_data <- list(
  studies = only_interesting_blocks,
  # colocs = all_coloc_results_for_study,
  grouped_colocs = raw_coloc_results
)

coloc_result_json <- jsonlite::toJSON(combined_data, pretty = T)
write(coloc_result_json, '~/coloc_result.json')

