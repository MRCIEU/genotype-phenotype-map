## Date: 21-02-2025
## A.L.Hanson

library(dplyr)

### De-duplicating binary traits within and between az, bm and gb UKB WES-derived rare variant summary stats

# Read metadata
az_meta <- "/local-scratch/data/ukb-seq/downloads/azexwas/ukb-wes-az-metadata.tsv"
bm_meta <- "/local-scratch/data/ukb-seq/downloads/backmanexwas/ukb-wes-bm-metadata.tsv"
gb_meta <- "/local-scratch/data/ukb-seq/downloads/genebass/rare/ukb-wes-gb-metadata.tsv"

meta <- lapply(list(az_meta, bm_meta, gb_meta), data.table::fread)
names(meta) <- c("az", "bm", "gb")

#### ----------- STEP 1: De-duplicate within each study set ------------ #####

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

# Source of report 
az_search_sor_1 <- paste(paste("Source of report of", unlist(lapply(strsplit(az_ICD10[1:round(length(az_ICD10)/2)],"#"), function(x){x[2]}))), collapse = "|")
az_search_sor_2 <- paste(paste("Source of report of", unlist(lapply(strsplit(az_ICD10[round(length(az_ICD10)/2+1):length(az_ICD10)],"#"), function(x){x[2]}))), collapse = "|")

## Remove "Source of report" if "union" trait exists
az_search_sor_3 <- paste(paste("Source of report of", unlist(lapply(strsplit(az_union[1:round(length(az_union)/2)],"#"), function(x){x[2]}))), collapse = "|")
az_search_sor_4 <- paste(paste("Source of report of", unlist(lapply(strsplit(az_union[round(length(az_union)/2+1):length(az_union)],"#"), function(x){x[2]}))), collapse = "|")

az_match_sor <- c(grep(az_search_sor_1, meta$az$trait_original, value = T), grep(az_search_sor_2, meta$az$trait_original, value = T),
    grep(az_search_sor_3, meta$az$trait_original, value = T), grep(az_search_sor_4, meta$az$trait_original, value = T)) |> unique()  
# 709
az_sor <- grep("Source of report of", meta$az$trait_original , value = T) #743 (34 to retain)

## az traits to remove:
az_out <- c(az_40001,az_40002,az_20002,az_match_40006,az_match_union, az_match_sor)
az_out_studies <- meta$az |> dplyr::filter(trait_original %in% az_out) |> dplyr::pull(study_name) #5206

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
out_studies <- c(az_out_studies, bm_out_studies, gb_out_studies)

meta_all <- do.call("rbind", meta)
meta_filt <- meta_all |> dplyr::filter(!(study_name %in% out_studies)) #5810 removed

#meta_filt |> dplyr::filter(grepl("psoriasis", trait))
#meta_filt |> dplyr::filter(grepl("eczema", trait))

#### ----------- STEP 2: De-duplicate across study sets ------------ #####

# Remove duplicate ICD10 codes (n=1904 from bm) from az and gb
# Search terms
bm_ICD10_code <- paste(paste0("\\b",sub("\\.", "\\\\.", sub(":","",unlist(lapply(strsplit(bm_ICD10, " "), function(x){x[2]})))),"\\b"), collapse = "|")

# az: identify "Source or report", #41202 and #union with matched codes --> REMOVE
# az: keep "Blocks" of traits (these are unique phenotype sets)
# az studies do not have the "." between codes and subcodes (e.g. L40.5 is L405)
bm_ICD10_code_az <- paste(paste0("\\b",sub("\\.", "", sub(":","",unlist(lapply(strsplit(bm_ICD10, " "), function(x){x[2]})))),"\\b"), collapse = "|")

az_check <- meta_filt |> filter(grepl(bm_ICD10_code_az, trait_original)) |> filter(data_format == "azphewas") |> pull(trait_original)
az_out_r2 <- grep("Source|41202#[A-Z][0-9]|union", az_check, value = T)
az_out_studies_r2 <- meta$az |> dplyr::filter(trait_original %in% az_out_r2) |> dplyr::pull(study_name) #1280

# gb: identify "Date first reported" and numerical codes matched with codes --> REMOVE
# gb: keep "Operartive procedures"
gb_check <- meta_filt |> filter(grepl(bm_ICD10_code, trait_original)) |> filter(data_format == "genebass") |> pull(trait_original)
gb_out_r2 <- grep("Date|^[A-Z][0-9]", gb_check, value = T)
gb_out_studies_r2 <- meta$gb |> dplyr::filter(trait_original %in% gb_out_r2) |> dplyr::pull(study_name) #450

# Remove 
meta_filt_r2 <- meta_filt |> dplyr::filter(!(study_name %in% c(az_out_studies_r2, gb_out_studies_r2))) #1730 removed

#meta_filt_r2 |> dplyr::filter(grepl("psoriasis", trait, ignore.case = T))
#meta_filt_r2 |> dplyr::filter(grepl("multiple sclerosis", trait, ignore.case = T))

# Remove duplicate gb "Date first reported" if az "Source of report" or "union" exist
# Search terms
az_sor_code <- paste(paste("Date", unlist(lapply(strsplit(grep("Source of report", meta_filt_r2$trait_original, value = T), " "), function(x){x[5]})), "first reported"), collapse = "|")
az_union_code <- paste(paste("Date", unlist(lapply(strsplit(grep("^union", meta_filt_r2$trait_original, value = T), "#"), function(x){x[2]})), "first reported"), collapse = "|")

# gb: identify "Date first reported" and numerical codes matched with az "Source of report" or "union" codes --> REMOVE
gb_check_again <- meta_filt_r2 |> filter(grepl(az_sor_code, trait_original) | grepl(az_union_code, trait_original) ) |> filter(data_format == "genebass") |> pull(trait_original)
gb_out_studies_r3 <- meta$gb |> dplyr::filter(trait_original %in% gb_check_again) |> dplyr::pull(study_name) #63

# Remove 
meta_filt_r3 <- meta_filt_r2 |> dplyr::filter(!(study_name %in% c(gb_out_studies_r3))) # 63 removed

#meta_filt_r3 |> dplyr::filter(grepl("psoriasis", trait, ignore.case = T))
#meta_filt_r3 |> dplyr::filter(grepl("type 1 diabetes", trait, ignore.case = T))
#meta_filt_r3 |> dplyr::filter(grepl("haemoglobin", trait, ignore.case = T))

#### ----------- STEP 3: De-duplicate across study sets ------------ #####

# Use efo trait mappings from Gib (turns out these are actually to coarse, ignoring for now)

trait_mapping <- read.csv("https://raw.githubusercontent.com/MRCIEU/finetune-efo/refs/heads/main/manual-mapping/trait-efo.csv", header = T, row.names = 1)

mapping_rare <- meta_filt_r3 |> 
    left_join(trait_mapping, by = c("study_name" = "trait_id"))

# data.table::fwrite(mapping_rare, "/local-scratch/projects/genotype-phenotype-map/data/trait_formatting/trait_mapping_rare.tsv", sep = "\t")  

# Ignore proteins and NMR metabolites, brain MRI measures (az only), opperative procedures (gb only)
rare_proteins <- mapping_rare |> filter(data_type == "protein") |> pull(study_name) #109
rare_nmr <- mapping_rare |> filter(grepl("^NMR ", trait.x)) |> pull(study_name) #325
rare_opp <- mapping_rare |> filter(grepl("Operative procedure", trait.x)) |> pull(study_name) #1522

# Remove smallest set of Backman continuous traits (n = 289)
bm_out_cont <- mapping_rare |> filter(data_format == "backman" & category == "continuous") |> pull(study_name) #289
# Remove self reported cancer codes
code20001_out <- mapping_rare |> filter(grepl("^20001|\\(20001\\)", trait_original)) |> pull(study_name) #77

mapping_rare <- mapping_rare |> filter(!study_name %in% c(rare_proteins, rare_nmr, rare_opp, bm_out_cont, code20001_out))

# First group identical trait strings and keep study with greatest sample size
mapping_rare <- mapping_rare |> arrange(trait.x, desc(sample_size))
mapping_rare <- mapping_rare |> filter(duplicated(trait.x) == FALSE) # 9672 remaining

# Remove duplicate union traits
mapping_rare <- mapping_rare |> filter(duplicated(sub(" \\(union\\)","",trait.x)) == FALSE) # 9631

#### ----------- STEP 4: Write out list of de-duplicated traits and efo terms ------------ #####
meta_filt_r4 <- meta_filt_r3 |> dplyr::filter(study_name %in% c(mapping_rare$study_name,rare_proteins,rare_nmr,rare_opp)) # 14419 rare studies remaining

out <- meta_filt_r4 |> 
    left_join(trait_mapping, by = c("study_name" = "trait_id", "trait"))

data.table::fwrite(out, "/local-scratch/projects/genotype-phenotype-map/data/trait_formatting/deduplicated_rarestudies.tsv", sep = "\t")

### Studies to ignore
ignore <- meta_all$study_name[!(meta_all$study_name %in% out$study_name)]

write.table(c(az_out_studies, bm_out_studies, gb_out_studies), "pipeline_steps/data/ignore_studies_rare.tsv", row.names = F, col.names = F, quote = F)