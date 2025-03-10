## Date: 28-02-2025
## A.L.Hanson
### Paring rare and common variant UKB studies

library(dplyr)

rare_studies <- data.table::fread("/local-scratch/projects/genotype-phenotype-map/data/trait_cleaning/deduplicated_rarestudies.tsv")
# 14507 studies

common_studies <- data.table::fread("/local-scratch/projects/genotype-phenotype-map/results/studies_processed.tsv") |> 
   filter(variant_type == "common" & grepl("ukb", study_name, ignore.case = T))
# 9056 studies (ukb and ebi)
# 6260 ukb only

#### ----------- STEP 1: Identify nonsense/non-heritable/'non-health' traits from common variant studies (not included in rare) ------------ #####
# See /local-scratch/projects/genotype-phenotype-map/data/trait_cleaning/showcase_outcategories.txt
# See /local-scratch/projects/genotype-phenotype-map/data/trait_cleaning/showcase_outfields.txt

outfields <- data.table::fread("/local-scratch/projects/genotype-phenotype-map/data/trait_cleaning/showcase_outfields.txt")

# Remove fields
out_studies <- common_studies |> 
    filter(sub(":.*","",common_studies$trait) %in% outfields$Field_name) |>
    dplyr::select(study_name, trait)

# Write out:
#data.table::fwrite(out_studies, "/local-scratch/projects/genotype-phenotype-map/data/trait_cleaning/commonstudies_nonsensetraits.txt", sep = "\t")

common_studies <- common_studies |> filter(!(study_name %in% out_studies$study_name)) # 905 removed (5355 remaining)

# ----------- Match UKB-PPP datasets ---------- #
common_proteins <- common_studies |> filter(data_type == "protein") |> 
    dplyr::select(data_type,study_name,trait,gene) |>
    mutate(OID = unlist(lapply(strsplit(study_name, ":"), function(x){x[3]}))) #2940

rare_proteins <- rare_studies |> filter(data_type == "protein") |> 
    dplyr::select(data_type,study_name,trait,gene,fieldID) |>
    mutate(OID = unlist(lapply(strsplit(trait, " "), function(x){x[3]}))) #2941

common_rare_proteins <- left_join(common_proteins, rare_proteins,
    by = c("data_type","gene","OID"), suffix = c("_common","_rare")) |> 
    select(-OID) |>
    relocate(data_type,study_name_common, study_name_rare, trait_common, trait_rare, gene) #2940 direct matches


### Need to source exact field IDs


# # ----------- Direct string matches ---------- #
# common_rare_match <- inner_join(common_studies[,c("data_type","study_name","trait")], rare_studies[,c("data_type","study_name","trait","fieldID")],
#     by = c("data_type","trait"), suffix = c("_common","_rare")) |> 
#     rename(trait_common = trait) |>
#     mutate(trait_rare = trait_common) |> 
#     relocate(data_type, study_name_common, study_name_rare, trait_common, trait_rare, fieldID) #54 direct matches


# # ----------- Match ICD10 codes ---------- #
# common_ICD10 <- common_studies |> filter(grepl("main ICD10", trait)) |> 
#     dplyr::select(data_type,study_name,trait) |>
#     mutate(ICD10 = sub(":","",sub(".*(ICD10: [A-Z][0-9.]+).*", "\\1", trait))) #329

# rare_ICD10 <- rare_studies |> filter(grepl("ICD10", trait)) |> 
#     dplyr::select(data_type,study_name,trait, fieldID) |>
#     mutate(ICD10 = sub(".*(ICD10 [A-Z][0-9.]+):.*", "\\1", trait)) #1904

# common_rare_ICD10 <- inner_join(common_ICD10, rare_ICD10,
#     by = c("data_type","ICD10"), suffix = c("_common","_rare")) |> 
#     select(-ICD10) |>
#     relocate(data_type,study_name_common, study_name_rare, trait_common, trait_rare) #320 direct matches

# # ----------- Direct string matches ---------- #
# common_rare_match <- inner_join(common_studies[,c("data_type","study_name","trait")], rare_studies[,c("data_type","study_name","trait","fieldID")],
#     by = c("data_type","trait"), suffix = c("_common","_rare")) |> 
#     rename(trait_common = trait) |>
#     mutate(trait_rare = trait_common) |> 
#     relocate(data_type, study_name_common, study_name_rare, trait_common, trait_rare, fieldID) #54 direct matches

# common_remaining <- common_studies |> filter(!(study_name %in% 
#     c(common_rare_proteins$study_name_common, 
#     common_rare_ICD10$study_name_common, 
#     common_rare_match$study_name_common,
#     common_rare_field$study_name_common)))

# # ----------- Modified string matches ---------- #




















# ----------- Direct string matches ---------- #
common_rare_match <- inner_join(common_studies[,c("data_type","study_name","trait")], rare_studies[,c("data_type","study_name","trait")],
    by = c("data_type","trait"), suffix = c("_common","_rare")) |> 
    relocate(data_type,study_name_common,study_name_rare) |>
    rename(trait_common = trait) |>
    mutate(trait_rare = trait_common) #287 direct matches





## What kinds of traits are still unmatched?
common_studies |> filter(!study_name %in% 
    c(common_rare_proteins$study_name_common,
    common_rare_ICD10$study_name_common,
    common_rare_field$study_name_common)) |> pull(trait)


