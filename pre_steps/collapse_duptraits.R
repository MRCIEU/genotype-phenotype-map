## Date: 20-10-2024
## A. Hanson
## Collapsing duplicate study traits to derive index and duplicate study IDs

library(data.table)
library(dplyr)
library(here)

here()

# Open GWAS ebi-a studies from the EBI GWAS Catalog
studies <- fread("/local-scratch/projects/genotype-phenotype-map/data/pipeline_metadata/studies_to_process.tsv")

# EBI experimental factor ontology (EFO) trait ID
efos <- fread("/local-scratch/projects/genotype-phenotype-map/results/phenotype_categorisation/study_trait_map.tsv")

study_efos <- left_join(
  studies[,.(study_name, data_type, data_format, ancestry, sample_size, trait)],
  efos[,.(study_name, efoTraits, main_efo_trait)])

# How many clumped top hits per study?
file_names <- studies$study_name
file_paths <- file.path("/local-scratch/data/opengwas/igd", file_names, "clump.txt")

study_efos$num_clumped <- unlist(lapply(file_paths, function(x){
  as.numeric(sub(" .*","",system(paste("wc -l", x), intern = TRUE)))
  }))

# Collapse by identical trait
# Annotate if study also has greatest number of clustered hits
study_efos_unique <- study_efos[,{
  study_name_collapsed <- if(.N == 1){
    ""
  } else {
    paste(study_name[-which.max(sample_size)], collapse = ",")
  }
  top_num_clumped <- if(.N == 1){
    TRUE
  } else{
    which.max(sample_size) == which.max(num_clumped)
  }
  .SD[sample_size == max(sample_size), .(study_name, data_type, data_format, ancestry, sample_size, num_clumped, efoTraits, study_name_collapsed, top_num_clumped)]
}, by = trait]

print("Writing output...")

# Write out info, index study and duplicate study IDs
fwrite(study_efos_unique, 
       here("pipeline_steps/data/collapsed_studies_info.tsv"), sep = "\t")

fwrite(as.data.frame(study_efos_unique[,study_name]), 
       here("pipeline_steps/data/index_studies.tsv"), col.names = F)

fwrite(as.data.frame(
  unlist(strsplit(paste(study_efos_unique[nchar(study_name_collapsed)>0,study_name_collapsed], collapse = ","),","))),
  here("pipeline_steps/data/ignore_studies.tsv"), col.names = F)
