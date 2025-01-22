# Update of rare variant study metadata files with data_type, category, gene, tissue and trait information

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
    extracted_location = paste0("/local-scratch/projects/genotype-phenotype-map/data/study_sequencing/", studies_list$study_id),
    reference_build = studies_list$build,
    p_value_threshold = 1.5e-4,
    gene = NA,
    probe = NA,
    tissue = NA,
    variant_type = "rare_exome"
)

### ------- Add data_type and category ------- ###

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
# Read in Genebass study metadata
gb_studies <- vroom::vroom(file.path(data_dir,"genebass/genebass_wes_studymetadata.tsv"), show_col_types = F)

# Populate category
studies_list_out[match(gb_studies$file_name, basename(studies_list_out$study_location)), "category"] <- gb_studies$trait_type
studies_list_out$category[studies_list_out$category %in% c("icd10","icd_first_occurrence")] <- "binary"

# Annotate category (binary, continuous) for Backman studies
# This information is provided in the trait_string in the study metadata file 
bm_studies <- vroom::vroom(meta_bm, show_col_types = F)
bm_studies$category <- ifelse(grepl("BIN", bm_studies$trait_string), "binary", "continuous")

# Populate category
studies_list_out[match(bm_studies$study_id, studies_list_out$study_name), "category"] <- bm_studies$category

## ------- Neaten up trait names (predominantly for az studies) -------- ##
# Move original traits to column trait_original
studies_list_out$trait_original <- studies_list_out$trait

studies_list_out <- studies_list_out |>
    dplyr::mutate(trait = ifelse(data_format == "azphewas" & data_type == "phenotype", sub("^[0-9].*#","",trait_original), trait)) |>
    dplyr::mutate(trait = ifelse(data_format == "azphewas" & data_type == "phenotype" & grepl("^union", trait), paste(sub("union.*#","",trait_original), "(union)"), trait)) |>
    dplyr::mutate(trait = ifelse(data_format == "azphewas" & data_type == "protein", gsub("#"," ",sub("^[0-9]+#[0-9]+#","",trait_original)), trait)) |>
    dplyr::mutate(gene = ifelse(data_format == "azphewas" & data_type == "protein", sub(" .*","",trait), NA)) |>
    dplyr::mutate(tissue = ifelse(data_format == "azphewas" & data_type == "protein", "plasma", NA))

# Write out metadata file
vroom::vroom_write(studies_list_out |> dplyr::filter(data_format == "azphewas"), file.path(data_dir,"azexwas/ukb-wes-az-metadata.tsv"))
vroom::vroom_write(studies_list_out |> dplyr::filter(data_format == "backman"), file.path(data_dir,"backmanexwas/ukb-wes-bm-metadata.tsv"))
vroom::vroom_write(studies_list_out |> dplyr::filter(data_format == "genebass"), file.path(data_dir,"genebass/rare/ukb-wes-gb-metadata.tsv"))