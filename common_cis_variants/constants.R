options(error = function() traceback(20))

data_dir <- Sys.getenv('DATA_DIR')
results_dir <- Sys.getenv('RESULTS_DIR')
pipeline_metadata_dir <- paste0(data_dir, "pipeline_metadata/")
ld_block_dir <- paste0(data_dir, "ld_blocks/")
extracted_study_dir <- paste0(data_dir, "study/")

data_types <- list(phenotype="phenotype", proteomics="proteomics", metabolomics="metabolomics", eqtl="eqtl")
databases <- list(opengwas="opengwas")
data_source <- list(ukb="UK Biobank", gtex="GTEx", finngen="Finn Gen", bbj="Biobank Japan", eqtl_gen="eQTL Gen", ukb_ppp="UK Biobank Proteomic...")
ancestry_map <- list(EUR="European", EAS="East Asian", AFR="African")
reverse_ancestry_map <- setNames(names(ancestry_map), ancestry_map)
