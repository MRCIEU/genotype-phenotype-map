studies_processed <- vroom::vroom(paste0('scratch/results/',timestamp,'/studies_processed.tsv'), show_col_types=F)
all_study_results <- vroom::vroom(paste0('scratch/results/',timestamp,'/all_study_regions.tsv'), show_col_types=F)
coloc_results <- vroom::vroom(paste0('scratch/results/',timestamp,'/coloc_results.tsv'), show_col_types=F)
studies <- merge(all_study_results, studies_processed, by.x='study', by.y='study_name')
ebi <- vroom::vroom(paste0('scratch/results/phenotype_categorisation/ebi_gwas_metadata.tsv'), col_select=c('accessionId', 'reportedTrait', 'efoTraits'), show_col_types=F) |>
  tidyr::separate(sep=',', col = "efoTraits", into = c("efo_trait_1", "efo_trait_2"), remove = F) |>
  dplyr::filter(accessionId %in% studies$study)


specific_gene <- 'RAD9A'

study_a_match <- match(coloc_results$unique_study_a, studies$unique_study_id)
study_b_match <- match(coloc_results$unique_study_b, studies$unique_study_id)

pheno_study_a_match <- match(coloc_results$unique_study_a, phenotype_studies$unique_study_id)
pheno_study_b_match <- match(coloc_results$unique_study_b, phenotype_studies$unique_study_id)

coloc_results$gene_a <- studies$known_gene[study_a_match]
coloc_results$gene_b <- studies$known_gene[study_b_match]
coloc_results$data_type_a <- studies$data_type[study_a_match]
coloc_results$data_type_b <- studies$data_type[study_b_match]
coloc_results$efo_trait_1a <- phenotype_studies$efo_trait_1[pheno_study_a_match]
coloc_results$efo_trait_1b <- phenotype_studies$efo_trait_1[pheno_study_b_match]

known_gene_results <- dplyr::filter(coloc_results, !is.na(gene_a) & gene_a == specific_gene) |>
dplyr::select(-efo_trait_2a, -efo_trait_3a, -efo_trait_4a, -efo_trait_2b, -efo_trait_3b, -efo_trait_4b, -data_type_a, -data_type_b, -efo_traits_b)
implied_gene_results <- dplyr::filter(coloc_results, unique_study_a %in% known_gene_results$unique_study_b & data_type_a == 'phenotype')
