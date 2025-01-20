data_dir <- "/local-scratch/data/ukb-seq/downloads"

# Individual metadata file locations
meta_bm <- file.path(data_dir,"backmanexwas/ukb-wes-bm-studies.tsv")
meta_gb <- file.path(data_dir,"genebass/rare/ukb-wes-gb-studies.tsv")
meta_az <- file.path(data_dir,"azexwas/ukb-wes-az-studies.tsv")

studies_list <- do.call("rbind", lapply(list(meta_bm, meta_gb, meta_az),
    function(x){vroom::vroom(x, show_col_types = F)}))

# Format studies list
studies_list_out <- data.frame(
    data_type = NA,
    data_format = studies_list$source_study,
    study_name = studies_list$study_id,
    trait = studies_list$trait,
    ancestry= studies_list$ancestry,
    sample_size = studies_list$sample_size,
    category = NA,
    study_location = studies_list$file_path,
    extracted_location = paste0("/local-scratch/projects/genotype-phenotype-map/data/study_sequencing", studies_list$study_id, "/"),
    reference_build = studies_list$build,
    p_value_threshold = 1.5e-4,
    gene = NA,
    probe = NA,
    tissue = NA
)

# Annotate data_type (protein or phenotype) for studies from in AZ PheWAS Portal data
# Read in AZ study class infomation
az_studies <- vroom::vroom(file.path(data_dir,"azexwas/azphewas-unique-phenocodes-classes.csv"), show_col_types = F)
az_studies <- az_studies |> dplyr::filter(trait %in% studies_list_out$trait)

# Populate data_type and category
studies_list_out[match(az_studies$trait,studies_list_out$trait),"data_type"] <- az_studies$data_type
studies_list_out[match(az_studies$trait,studies_list_out$trait),"category"] <- az_studies$category

# All remaining data_type = "phenotype"
studies_list_out$data_type[is.na(studies_list_out$data_type)] <- "phenotype"

# Annotate category (binary, categorical, continuous) for Genebass studies
gb_studies <- vroom::vroom(file.path(data_dir,"genebass/genebass_wes_studymetadata.tsv"), show_col_types = F)

# Populate category
studies_list_out[match(gb_studies$file_name, basename(studies_list_out$study_location)), "category"] <- gb_studies$trait_type
studies_list_out$category[studies_list_out$category %in% c("icd10","icd_first_occurrence")] <- "binary"

## Need to 
# add categories for backman datasets
# fix trait names