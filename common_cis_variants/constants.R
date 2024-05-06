opengwas_dir <- "/mnt/storage/private/mrcieu/data/IGD/data/public/"
data_dir <- "./"
results_dir <- "./"
pipeline_metadata_dir <- paste0(data_dir, "pipeline_metadata/")

data_types <- list(phenotypes="phenotypes", proteomics="proteomics", metabolomics="metabolomics")
data_sources <- list(opengwas="opengwas")
databases<- list(ukb="UK Biobank", gtex="GTEx", finngen="Finn Gen", bbj="Biobank Japan", eqtl_gen="eQTL Gen", ukb_ppp="UK Biobank Proteomic...")
ancestry_map <- list(EUR="European", EAS="East Asian", AFR="African")