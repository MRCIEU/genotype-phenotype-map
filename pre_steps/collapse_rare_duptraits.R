## Date: 21-02-2025
## A.L.Hanson
### De-duplicating traits within and between az, bm and gb UKB WES-derived rare variant summary stats

library(dplyr)

# Read metadata
az_meta <- "/local-scratch/data/ukb-seq/downloads/azexwas/ukb-wes-az-metadata.tsv"
bm_meta <- "/local-scratch/data/ukb-seq/downloads/backmanexwas/ukb-wes-bm-metadata.tsv"
gb_meta <- "/local-scratch/data/ukb-seq/downloads/genebass/rare/ukb-wes-gb-metadata.tsv"

meta <- lapply(list(az_meta, bm_meta, gb_meta), data.table::fread)
names(meta) <- c("az", "bm", "gb")

## Add UKB field IDs if available
meta$bm$fieldID <- ifelse(grepl("\\([0-9]+\\)$",meta$bm$trait_original), sub(".*\\(([0-9]+)\\)$","\\1",meta$bm$trait_original), NA)

gb_fields <- data.table::fread("/local-scratch/data/ukb-seq/downloads/genebass/genebass_wes_studymetadata.tsv")
meta$gb$fieldID <- gb_fields[match(gsub(".*/","",meta$gb$study_location), gb_fields$file_name),"phenocode"]

meta$az$fieldID <- ifelse(grepl("^[0-9]+#[A-Z].*", meta$az$trait_original), sub("#.*","", meta$az$trait_original), NA)

#### ----------- STEP 1: Identify unique traits across study sets ------------ #####

# For backman studies:
# Keep binary traits (including ICD10 codes #41202 and type of cancer #40006)
# Remove #20002 and #20001 self reported traits
# Remove treatment/medication codes (larger set in gb)
# Keep continuous traits (then exclude duplicates from larger set in gb)

bm_selfreported <- grep("^Non cancer illness code self reported|Self Report", meta$bm$trait_original , value = T) #365 OUT
bm_continuous <- meta$bm |> dplyr::filter(category == "continuous") |> pull(trait_original) # 289 KEEP
bm_medication <- grep("Treatment or Medication|medication for", meta$bm$trait_original, ignore.case = T, value = T) #257 OUT

bm_out_studies <- meta$bm |> dplyr::filter(trait_original %in% c(bm_selfreported, bm_medication)) |> dplyr::pull(study_name) #622 to remove

# For genebass studies:
# Remove binary traits (mostly "Date first reported" and cancer codes, larger set in bm)
# Keep operative procedures (main and secondary) and operation code, treatment/medication codes
# Remove remaining categorical traits (including self-reported traits)
# Keep continuous traits (including imaging measures)

gb_binary <- meta$gb |> dplyr::filter(category == "binary") |> pull(trait_original) #725 OUT
gb_opps <- grep("operat|^Treatment|^Medication", meta$gb$trait_original , value = T, ignore.case = T) #2139 KEEP
gb_categorical <- meta$gb |> dplyr::filter(category == "categorical" & !(trait_original %in% gb_opps)) |> pull(trait_original) #432 OUT

gb_continuous <- meta$gb |> dplyr::filter(category == "continuous") |> pull(trait_original) # 1233 

# Remove continuous variables alread in bm (based on field ID)
bm_continuous_match <- meta$bm |> dplyr::filter(trait_original %in% c(bm_continuous)) |> dplyr::pull(fieldID) |> na.exclude()
gb_continuous_match <- meta$gb |> dplyr::filter(category == "continuous" & fieldID %in% bm_continuous_match) |> pull(trait_original) #134 OUT

gb_out_studies <- meta$gb |> dplyr::filter(trait_original %in% c(gb_binary, gb_categorical, gb_continuous_match)) |> dplyr::pull(study_name) #1300 to remove

# For azexwas studies:
# Remove binary disease traits (keeping only bm ICD10 codes as master studies, there is too much duplication here otherwise)

# Keep "union" traits for binary phenotypes, these combine diagnosis calls across:
# #41202 ICD10 Diagnoses - main
# #40001 Primary cause of death
# #40002 Secondary cause of death
# #40006 Type of cancer
# #20002 Non-cancer illness, self reported)

# Keep ICD10 #41202 "blocks" which capture subtraits across chapters
az_union <- grep("^union", meta$az$trait_original , value = T, ignore.case = T) #4516 KEEP
az_block <- grep("^41202#Block", meta$az$trait_original , value = T, ignore.case = T) #157 KEEP

# Keep proteins and NMR metabolites
az_proteins <- meta$az |> filter(data_type == "protein") |> pull(trait_original) #2941 KEEP
az_nmr <- meta$az |> filter(grepl("^NMR ", trait)) |> pull(trait_original) #325 KEEP

# Remove remaining continuous traits (duplicated with gb)
az_continuous <- meta$az |> dplyr::filter(category == "continuous", !(trait_original %in% c(az_proteins,az_nmr))) |> pull(trait_original) #1602 OUT

az_out_studies <- meta$az |> dplyr::filter(!(trait_original %in% c(az_union, az_block, az_proteins, az_nmr))) |> dplyr::pull(study_name) #7031 to remove

#### ----------- STEP 2: Remove duplicate sets of studies ------------ #####

out_studies <- c(az_out_studies, bm_out_studies, gb_out_studies) #8953 to remove

meta_all <- do.call("rbind", meta)
meta_filt <- meta_all |> dplyr::filter(!(study_name %in% out_studies)) # 14507 rare variant studies remaining

### Studies to ignore
ignore <- meta_all$study_name[!(meta_all$study_name %in% meta_filt$study_name)] #8953

#### ----------- STEP 3: Merge with EFO terms from trait mappings ------------ #####

# Use efo trait mappings from Gib:
trait_mapping <- read.csv("https://raw.githubusercontent.com/MRCIEU/finetune-efo/refs/heads/main/manual-mapping/trait-efo.csv", header = T, row.names = 1)

mapping_rare <- meta_filt |> 
    left_join(trait_mapping, by = c("study_name" = "trait_id", "trait"))

# Check a few...
# mapping_rare |> dplyr::filter(grepl("psoriasis", trait))
# mapping_rare |> dplyr::filter(grepl("eczema", trait))
# mapping_rare |> dplyr::filter(grepl("haemoglobin", trait))
# mapping_rare |> dplyr::filter(grepl("diabetes", trait))

# Write
data.table::fwrite(mapping_rare, "/local-scratch/projects/genotype-phenotype-map/data/trait_cleaning/deduplicated_rarestudies.tsv", sep = "\t")
write.table(ignore, "pipeline_steps/data/ignore_studies_rare.tsv", row.names = F, col.names = F, quote = F)

# --------------------- Previous attempts... ----------------------- #

#az_cancer <- grep("^40006", meta$az$trait_original , value = T) #313 traits
#az_ICD10 <- grep("^41202", meta$az$trait_original , value = T) #3253 traits
#az_ICD10 <- az_ICD10[!(az_ICD10 %in% sub("40006","41202",az_cancer))] # 2966 ICD10 (not in cancer type)

# Search for equivalent traits in other diagnosis call sources:
# Primary cause of death
# az_search_40001_1 <- paste0("^\\Q",sub("41202","40001",az_ICD10[1:round(length(az_ICD10)/2)]),"\\E$", collapse = "|")
# az_search_40001_2 <- paste0("^\\Q",sub("41202","40001",az_ICD10[round(length(az_ICD10)/2+1):length(az_ICD10)]),"\\E$", collapse = "|")
# az_match_40001 <- c(grep(az_search_40001_1, meta$az$trait_original, value = T), grep(az_search_40001_2, meta$az$trait_original, value = T))
# # 260
# az_40001 <- grep("^40001", meta$az$trait_original , value = T) #263 traits (remove all)

# # Secondary cause of death
# az_search_40002_1 <- paste0("^\\Q",sub("41202","40002",az_ICD10[1:round(length(az_ICD10)/2)]),"\\E$", collapse = "|")
# az_search_40002_2 <- paste0("^\\Q",sub("41202","40002",az_ICD10[round(length(az_ICD10)/2+1):length(az_ICD10)]),"\\E$", collapse = "|")
# az_match_40002 <- c(grep(az_search_40002_1, meta$az$trait_original, value = T), grep(az_search_40002_2, meta$az$trait_original, value = T))
# # 350
# az_40002 <- grep("^40002", meta$az$trait_original , value = T) #361 traits (remove all)

# # Type of cancer
# az_search_40006_1 <- paste0("^\\Q",sub("41202","40006",az_ICD10[1:round(length(az_ICD10)/2)]),"\\E$", collapse = "|")
# az_search_40006_2 <- paste0("^\\Q",sub("41202","40006",az_ICD10[round(length(az_ICD10)/2+1):length(az_ICD10)]),"\\E$", collapse = "|")
# az_match_40006 <- c(grep(az_search_40006_1, meta$az$trait_original, value = T), grep(az_search_40006_2, meta$az$trait_original, value = T))
# # 287
# az_40006 <- grep("^40006", meta$az$trait_original , value = T) #313 traits (26 to retain)

# # Self report
# az_search_20002_1 <- paste0("^\\Q",sub("41202","20002",az_ICD10[1:round(length(az_ICD10)/2)]),"\\E$", collapse = "|")
# az_search_20002_2 <- paste0("^\\Q",sub("41202","20002",az_ICD10[round(length(az_ICD10)/2+1):length(az_ICD10)]),"\\E$", collapse = "|")
# az_match_20002 <- c(grep(az_search_20002_1, meta$az$trait_original, value = T), grep(az_search_20002_2, meta$az$trait_original, value = T))
# # none here
# az_20002 <- grep("^20002", meta$az$trait_original , value = T) #353 (remove all)

# # Union
# az_search_union_1 <- paste0("^\\Q",sub("41202","union",az_ICD10[1:round(length(az_ICD10)/2)]),"\\E$", collapse = "|")
# az_search_union_2 <- paste0("^\\Q",sub("41202","union",az_ICD10[round(length(az_ICD10)/2+1):length(az_ICD10)]),"\\E$", collapse = "|")
# az_match_union <- c(grep(az_search_union_1, meta$az$trait_original, value = T), grep(az_search_union_2, meta$az$trait_original, value = T))
# # 3233
# az_union <- grep("^union", meta$az$trait_original , value = T) #4501 (1268 to retain)

# # Source of report 
# az_search_sor_1 <- paste(paste("Source of report of", unlist(lapply(strsplit(az_ICD10[1:round(length(az_ICD10)/2)],"#"), function(x){x[2]}))), collapse = "|")
# az_search_sor_2 <- paste(paste("Source of report of", unlist(lapply(strsplit(az_ICD10[round(length(az_ICD10)/2+1):length(az_ICD10)],"#"), function(x){x[2]}))), collapse = "|")

# ## Remove "Source of report" if "union" trait exists
# az_search_sor_3 <- paste(paste("Source of report of", unlist(lapply(strsplit(az_union[1:round(length(az_union)/2)],"#"), function(x){x[2]}))), collapse = "|")
# az_search_sor_4 <- paste(paste("Source of report of", unlist(lapply(strsplit(az_union[round(length(az_union)/2+1):length(az_union)],"#"), function(x){x[2]}))), collapse = "|")

# az_match_sor <- c(grep(az_search_sor_1, meta$az$trait_original, value = T), grep(az_search_sor_2, meta$az$trait_original, value = T),
#     grep(az_search_sor_3, meta$az$trait_original, value = T), grep(az_search_sor_4, meta$az$trait_original, value = T)) |> unique()  
# # 709
# az_sor <- grep("Source of report of", meta$az$trait_original , value = T) #743 (34 to retain)

# ## az traits to remove:
# az_out <- c(az_40001,az_40002,az_20002,az_match_40006,az_match_union, az_match_sor)
# az_out_studies <- meta$az |> dplyr::filter(trait_original %in% az_out) |> dplyr::pull(study_name) #5206

# # For backman studies:
# # Keep only ICD10 phenotypes (remove self reported and composite codes)

# bm_ICD10 <- grep("^ICD10", meta$bm$trait_original , value = T)
# bm_selfreported <- grep("^Non cancer illness code self reported", meta$bm$trait_original , value = T) #347

# bm_out_studies <- meta$bm |> dplyr::filter(trait_original %in% bm_selfreported) |> dplyr::pull(study_name) #347

# # For genebass studies:
# # Remove self reported cancer and non cancer codes
# gb_selfreported <- grep("self-reported", meta$gb$trait_original , value = T)

# gb_out_studies <- meta$gb |> dplyr::filter(trait_original %in% gb_selfreported) |> dplyr::pull(study_name) #257

# Check
# out_studies <- c(az_out_studies, bm_out_studies, gb_out_studies)

# meta_all <- do.call("rbind", meta)
# meta_filt <- meta_all |> dplyr::filter(!(study_name %in% out_studies)) #5810 removed

#meta_filt |> dplyr::filter(grepl("psoriasis", trait))
#meta_filt |> dplyr::filter(grepl("eczema", trait))

#### ----------- STEP 2: De-duplicate across study sets ------------ #####

# Remove duplicate ICD10 codes (n=1904 from bm) from az and gb
# Search terms
# bm_ICD10_code <- paste(paste0("\\b",sub("\\.", "\\\\.", sub(":","",unlist(lapply(strsplit(bm_ICD10, " "), function(x){x[2]})))),"\\b"), collapse = "|")

# # az: identify "Source or report", #41202 and #union with matched codes --> REMOVE
# # az: keep "Blocks" of traits (these are unique phenotype sets)
# # az studies do not have the "." between codes and subcodes (e.g. L40.5 is L405 and R11.0 is R11)
# bm_ICD10_code_az <- paste(paste0("\\b",sub("\\.", "", sub(":","",unlist(lapply(strsplit(bm_ICD10, " "), function(x){x[2]})))),"\\b"), collapse = "|")

# az_check <- meta_filt |> filter(grepl(bm_ICD10_code_az, trait_original)) |> filter(data_format == "azphewas") |> pull(trait_original)
# az_out_r2 <- grep("Source|41202#[A-Z][0-9]|union", az_check, value = T)
# az_out_studies_r2 <- meta$az |> dplyr::filter(trait_original %in% az_out_r2) |> dplyr::pull(study_name) #1280

# # gb: identify "Date first reported" and numerical codes matched with codes --> REMOVE
# # gb: keep "Operartive procedures"
# gb_check <- meta_filt |> filter(grepl(bm_ICD10_code, trait_original)) |> filter(data_format == "genebass") |> pull(trait_original)
# gb_out_r2 <- grep("Date|^[A-Z][0-9]", gb_check, value = T)
# gb_out_studies_r2 <- meta$gb |> dplyr::filter(trait_original %in% gb_out_r2) |> dplyr::pull(study_name) #450

# # Remove 
# meta_filt_r2 <- meta_filt |> dplyr::filter(!(study_name %in% c(az_out_studies_r2, gb_out_studies_r2))) #1730 removed

# #meta_filt_r2 |> dplyr::filter(grepl("psoriasis", trait, ignore.case = T))
# #meta_filt_r2 |> dplyr::filter(grepl("multiple sclerosis", trait, ignore.case = T))

# # Remove duplicate gb "Date first reported" if az "Source of report" or "union" exist
# # Search terms
# az_sor_code <- paste(paste("Date", unlist(lapply(strsplit(grep("Source of report", meta_filt_r2$trait_original, value = T), " "), function(x){x[5]})), "first reported"), collapse = "|")
# az_union_code <- paste(paste("Date", unlist(lapply(strsplit(grep("^union", meta_filt_r2$trait_original, value = T), "#"), function(x){x[2]})), "first reported"), collapse = "|")

# # gb: identify "Date first reported" and numerical codes matched with az "Source of report" or "union" codes --> REMOVE
# gb_check_again <- meta_filt_r2 |> filter(grepl(az_sor_code, trait_original) | grepl(az_union_code, trait_original) ) |> filter(data_format == "genebass") |> pull(trait_original)
# gb_out_studies_r3 <- meta$gb |> dplyr::filter(trait_original %in% gb_check_again) |> dplyr::pull(study_name) #63

# # Remove 
# meta_filt_r3 <- meta_filt_r2 |> dplyr::filter(!(study_name %in% c(gb_out_studies_r3))) # 63 removed

#meta_filt_r3 |> dplyr::filter(grepl("psoriasis", trait, ignore.case = T))
#meta_filt_r3 |> dplyr::filter(grepl("type 1 diabetes", trait, ignore.case = T))
#meta_filt_r3 |> dplyr::filter(grepl("haemoglobin", trait, ignore.case = T))

#### ----------- STEP 3: Collapse like traits ------------ #####

# Use efo trait mappings from Gib (turns out these are actually to coarse, ignoring for now)

# trait_mapping <- read.csv("https://raw.githubusercontent.com/MRCIEU/finetune-efo/refs/heads/main/manual-mapping/trait-efo.csv", header = T, row.names = 1)

# mapping_rare <- meta_filt_r3 |> 
#     left_join(trait_mapping, by = c("study_name" = "trait_id"))

# # data.table::fwrite(mapping_rare, "/local-scratch/projects/genotype-phenotype-map/data/trait_formatting/trait_mapping_rare.tsv", sep = "\t")  

# # Ignore proteins and NMR metabolites, brain MRI measures (az only), opperative procedures (gb only)
# rare_proteins <- mapping_rare |> filter(data_type == "protein") |> pull(study_name) #109
# rare_nmr <- mapping_rare |> filter(grepl("^NMR ", trait.x)) |> pull(study_name) #325
# rare_opp <- mapping_rare |> filter(grepl("Operative procedure", trait.x)) |> pull(study_name) #1522

# # Remove smallest set of Backman continuous traits (n = 289)
# bm_out_cont <- mapping_rare |> filter(data_format == "backman" & category == "continuous") |> pull(study_name) #289
# # Remove smaller set of Treatment or Medication code traits (n=247), or "Medication for" (n=10) also in gb
# bm_out_med <- mapping_rare |> filter(grepl("Treatment or Medication|medication for", trait_original, ignore.case = T)) |> 
#     filter(data_format == "backman") |> 
#     pull(study_name)
# # Remove self reported cancer codes
# code20001_out <- mapping_rare |> filter(grepl("^20001|\\(20001\\)", trait_original)) |> pull(study_name) #77

# mapping_rare <- mapping_rare |> filter(!study_name %in% c(rare_proteins, rare_nmr, rare_opp, bm_out_cont, bm_out_med, code20001_out))

# # First group identical trait strings and keep study with greatest sample size
# mapping_rare <- mapping_rare |> arrange(trait.x, desc(sample_size))
# mapping_rare <- mapping_rare |> filter(duplicated(trait.x) == FALSE) # 9416 remaining

# # Remove duplicate union traits
# mapping_rare <- mapping_rare |> filter(duplicated(sub(" \\(union\\)","",trait.x)) == FALSE) # 9375

#### ----------- STEP 4: Write out list of de-duplicated traits and efo terms ------------ #####
# meta_filt_r4 <- meta_filt_r3 |> dplyr::filter(study_name %in% c(mapping_rare$study_name,rare_proteins,rare_nmr,rare_opp)) # 14173 rare studies remaining

# out <- meta_filt_r4 |> 
#     left_join(trait_mapping, by = c("study_name" = "trait_id", "trait"))

# out |> dplyr::filter(grepl("psoriasis", trait, ignore.case = T))
# out |> dplyr::filter(grepl("diabetes", trait, ignore.case = T))
#meta_filt_r3 |> dplyr::filter(grepl("haemoglobin", trait, ignore.case = T))

# data.table::fwrite(out, "/local-scratch/projects/genotype-phenotype-map/data/trait_formatting/deduplicated_rarestudies.tsv", sep = "\t")

# ### Studies to ignore
# ignore <- meta_all$study_name[!(meta_all$study_name %in% out$study_name)]

# write.table(ignore, "pipeline_steps/data/ignore_studies_rare.tsv", row.names = F, col.names = F, quote = F)