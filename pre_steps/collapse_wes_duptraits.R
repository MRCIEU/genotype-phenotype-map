## Date: 21-02-2025
## A.L.Hanson

### De-duplicating binary traits within az, bm and gb UKB WES-derived rare variant summary stats

# Read metadata
az_meta <- "/local-scratch/data/ukb-seq/downloads/azexwas/ukb-wes-az-metadata.tsv"
bm_meta <- "/local-scratch/data/ukb-seq/downloads/backmanexwas/ukb-wes-bm-metadata.tsv"
gb_meta <- "/local-scratch/data/ukb-seq/downloads/genebass/rare/ukb-wes-gb-metadata.tsv"

meta <- lapply(list(az_meta, bm_meta, gb_meta), data.table::fread)
names(meta) <- c("az", "bm", "gb")

# For azexwas studies:
# "union" traits for binary phenotypes, these combine diagnosis calls across:
# #41202 ICD10 Diagnoses - main
# #40001 Primary cause of death
# #40002 Secondary cause of death
# #40006 Type of cancer
# #20002 Non-cancer illness, self reported
# But have substantially fewer included controls (presumably due to missing data?)

# Keep #41202 ICD10 diagnoses and exclude those traits for which a corresponding union, self reported, cancer registry, or death phenotype exists
# Remove all self reported phenoytpes, primary and secondary cause of death calls
# Keep union and cancer calls only if a corresponding ICD10 codes doesn't exist

#az_cancer <- grep("^40006", meta$az$trait_original , value = T) #313 traits
az_ICD10 <- grep("^41202", meta$az$trait_original , value = T) #3253 traits
#az_ICD10 <- az_ICD10[!(az_ICD10 %in% sub("40006","41202",az_cancer))] # 2966 ICD10 (not in cancer type)

# Search for equivalent traits in other diagnosis call sources:
# Primary cause of death
az_search_40001_1 <- paste0("^\\Q",sub("41202","40001",az_ICD10[1:round(length(az_ICD10)/2)]),"\\E$", collapse = "|")
az_search_40001_2 <- paste0("^\\Q",sub("41202","40001",az_ICD10[round(length(az_ICD10)/2+1):length(az_ICD10)]),"\\E$", collapse = "|")
az_match_40001 <- c(grep(az_search_40001_1, meta$az$trait_original, value = T), grep(az_search_40001_2, meta$az$trait_original, value = T))
# 260
az_40001 <- grep("^40001", meta$az$trait_original , value = T) #263 traits (remove all)

# Secondary cause of death
az_search_40002_1 <- paste0("^\\Q",sub("41202","40002",az_ICD10[1:round(length(az_ICD10)/2)]),"\\E$", collapse = "|")
az_search_40002_2 <- paste0("^\\Q",sub("41202","40002",az_ICD10[round(length(az_ICD10)/2+1):length(az_ICD10)]),"\\E$", collapse = "|")
az_match_40002 <- c(grep(az_search_40002_1, meta$az$trait_original, value = T), grep(az_search_40002_2, meta$az$trait_original, value = T))
# 350
az_40002 <- grep("^40002", meta$az$trait_original , value = T) #361 traits (remove all)

# Type of cancer
az_search_40006_1 <- paste0("^\\Q",sub("41202","40006",az_ICD10[1:round(length(az_ICD10)/2)]),"\\E$", collapse = "|")
az_search_40006_2 <- paste0("^\\Q",sub("41202","40006",az_ICD10[round(length(az_ICD10)/2+1):length(az_ICD10)]),"\\E$", collapse = "|")
az_match_40006 <- c(grep(az_search_40006_1, meta$az$trait_original, value = T), grep(az_search_40006_2, meta$az$trait_original, value = T))
# 287
az_40006 <- grep("^40006", meta$az$trait_original , value = T) #313 traits (26 to retain)

# Self report
az_search_20002_1 <- paste0("^\\Q",sub("41202","20002",az_ICD10[1:round(length(az_ICD10)/2)]),"\\E$", collapse = "|")
az_search_20002_2 <- paste0("^\\Q",sub("41202","20002",az_ICD10[round(length(az_ICD10)/2+1):length(az_ICD10)]),"\\E$", collapse = "|")
az_match_20002 <- c(grep(az_search_20002_1, meta$az$trait_original, value = T), grep(az_search_20002_2, meta$az$trait_original, value = T))
# none here
az_20002 <- grep("^20002", meta$az$trait_original , value = T) #353 (remove all)

# Union
az_search_union_1 <- paste0("^\\Q",sub("41202","union",az_ICD10[1:round(length(az_ICD10)/2)]),"\\E$", collapse = "|")
az_search_union_2 <- paste0("^\\Q",sub("41202","union",az_ICD10[round(length(az_ICD10)/2+1):length(az_ICD10)]),"\\E$", collapse = "|")
az_match_union <- c(grep(az_search_union_1, meta$az$trait_original, value = T), grep(az_search_union_2, meta$az$trait_original, value = T))
# 3233
az_union <- grep("^union", meta$az$trait_original , value = T) #4501 (1268 to retain)

## az traits to remove:
az_out <- c(az_40001,az_40002,az_20002,az_match_40006,az_match_union)
az_out_studies <- meta$az |> dplyr::filter(trait_original %in% az_out) |> dplyr::pull(study_name) #4497

# For backman studies:
# Keep only ICD10 phenotypes (remove self reported and composite codes)

bm_ICD10 <- grep("^ICD10", meta$bm$trait_original , value = T)
bm_selfreported <- grep("^Non cancer illness code self reported", meta$bm$trait_original , value = T) #347

bm_out_studies <- meta$bm |> dplyr::filter(trait_original %in% bm_selfreported) |> dplyr::pull(study_name) #347

# For genebass studies:
# Remove self reported cancer and non cancer codes
gb_selfreported <- grep("self-reported", meta$gb$trait_original , value = T)

gb_out_studies <- meta$gb |> dplyr::filter(trait_original %in% gb_selfreported) |> dplyr::pull(study_name) #257

# Check
meta_all <- do.call("rbind", meta)
meta_filt <- meta_all |> dplyr::filter(!(study_name %in% c(az_out_studies, bm_out_studies, gb_out_studies))) #5101 removed

#meta_filt |> dplyr::filter(grepl("psoriasis", trait))
#meta_filt |> dplyr::filter(grepl("eczema", trait))

write.table(c(az_out_studies, bm_out_studies, gb_out_studies), "pipeline_steps/data/ignore_studies_rare.tsv", row.names = F, col.names = F, quote = F)