opengwas_dir <- "/mnt/storage/private/mrcieu/data/IGD/data/public/"
data_dir <- "/user/work/wt23152/test/data"
results_dir <- "/user/work/wt23152/test/results"
pipeline_metadata_dir <- paste0(data_dir, "pipeline_metadata/")

data_types <- list(phenotype="phenotype", proteomics="proteomics", metabolomics="metabolomics")
databases <- list(opengwas="opengwas")
databases<- list(ukb="UK Biobank", gtex="GTEx", finngen="Finn Gen", bbj="Biobank Japan", eqtl_gen="eQTL Gen", ukb_ppp="UK Biobank Proteomic...")
ancestry_map <- list(EUR="European", EAS="East Asian", AFR="African")
