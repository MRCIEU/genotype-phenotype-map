## Date: 28-02-2025
## A.L.Hanson
### Paring rare and common variant UKB studies

library(dplyr)
library(here)

rare_studies <- data.table::fread("/local-scratch/projects/genotype-phenotype-map/data/trait_cleaning/deduplicated_rarestudies.tsv") # nolint: line_length_linter.
# 14507 studies

common_studies <- data.table::fread("/local-scratch/projects/genotype-phenotype-map/results/studies_processed.tsv") |>
  filter(variant_type == "common" & grepl("ukb", study_name, ignore.case = T))
# 9056 studies (ukb and ebi - there are UKB derived sum stats in ebi- studies
# but there is no clear way to find which ones these are?)
# 6260 ukb only

# UKB studies to source .json files for (ukb-b and ukb-d)
ukb_studies <- common_studies |>
  filter(grepl("ukb-*", study_name)) |>
  pull(study_name)

# data.table::fwrite(common_studies |> filter(grepl("ukb-*", study_name)) |> pull(study_name) |> as.data.frame(),
#   here("pre_steps/data/opengwas_ukbstudies.txt"),
#   row.names = F, quote = F)

#### ----------- Identify nonsense/non-heritable/'non-health'
# traits from common variant studies (not included in rare) ------------ #####
# See /local-scratch/projects/genotype-phenotype-map/data/trait_cleaning/showcase_outcategories.txt
# See /local-scratch/projects/genotype-phenotype-map/data/trait_cleaning/showcase_outfields.txt

outfields <- data.table::fread("/local-scratch/projects/genotype-phenotype-map/data/trait_cleaning/showcase_outfields.txt") # nolint: line_length_linter.

# Remove fields
out_studies <- common_studies |>
  filter(grepl("Job SOC coding", trait) | sub(":.*", "", common_studies$trait) %in% outfields$Field_name) |>
  dplyr::select(study_name, trait)

# Write out:
data.table::fwrite(
  out_studies,
  "/local-scratch/projects/genotype-phenotype-map/data/trait_cleaning/commonstudies_nonsensetraits.txt", sep = "\t"
)

common_studies <- common_studies |> filter(!(study_name %in% out_studies$study_name)) # 980 removed (5281 remaining)

# ----------- Match UKB-PPP datasets ---------- #
common_proteins <- common_studies |>
  filter(data_type == "protein") |>
  dplyr::select(data_type, study_name, trait, gene) |>
  mutate(OID = unlist(lapply(strsplit(study_name, ":"), function(x) {
    return(x[3])
  }))) # 2940

rare_proteins <- rare_studies |>
  filter(data_type == "protein") |>
  dplyr::select(data_type, study_name, trait, gene, fieldID) |>
  mutate(OID = unlist(lapply(strsplit(trait, " "), function(x) {
    return(x[3])
  }))) # 2941

common_rare_proteins <- left_join(common_proteins, rare_proteins,
  by = c("data_type", "gene", "OID"), suffix = c("_common", "_rare")
) |>
  select(-OID) |>
  relocate(data_type, study_name_common, study_name_rare, trait_common, trait_rare, gene) # 2940 direct matches

# ----------- Match ICD10 codes ---------- #
common_ICD10 <- common_studies |>
  filter(grepl("main ICD10", trait)) |>
  dplyr::select(data_type, study_name, trait) |>
  mutate(ICD10 = sub(":", "", sub(".*(ICD10: [A-Z][0-9.]+).*", "\\1", trait))) # 329

rare_ICD10 <- rare_studies |>
  filter(grepl("ICD10", trait)) |>
  dplyr::select(data_type, study_name, trait, fieldID) |>
  mutate(ICD10 = sub(".*(ICD10 [A-Z][0-9.]+):.*", "\\1", trait)) # 1904

common_rare_ICD10 <- inner_join(common_ICD10, rare_ICD10,
  by = c("data_type", "ICD10"), suffix = c("_common", "_rare")
) |>
  select(-ICD10) |>
  relocate(data_type, study_name_common, study_name_rare, trait_common, trait_rare) # 320 direct matches

# ----------- Match treatment/medication codes ---------- #
common_med <- common_studies |>
  filter(grepl("Treatment/medication code", trait)) |>
  dplyr::select(data_type, study_name, trait) |>
  mutate(med = sub(".*: ", "", trait)) # 144

rare_med <- rare_studies |>
  filter(grepl("Treatment/medication code", trait)) |>
  dplyr::select(data_type, study_name, trait, fieldID) |>
  mutate(med = sub(".*\\((.*)\\)", "\\1", trait)) # 422

common_rare_med <- inner_join(common_med, rare_med,
  by = c("data_type", "med"), suffix = c("_common", "_rare")
) |>
  select(-med) |>
  relocate(data_type, study_name_common, study_name_rare, trait_common, trait_rare) # 137 direct matches

# ----------- Operative procedures ---------- #
common_opp <- common_studies |>
  filter(grepl("Operative procedures", trait)) |>
  dplyr::select(data_type, study_name, trait) |>
  mutate(opp = trait) # 294

rare_opp <- rare_studies |>
  filter(grepl("Operative procedures", trait)) |>
  dplyr::select(data_type, study_name, trait, fieldID) |>
  mutate(opp = sub(")$", "", sub("OPCS4 \\(", "OPCS: ", trait))) # 1522

common_rare_opp <- inner_join(common_opp, rare_opp,
  by = c("data_type", "opp"), suffix = c("_common", "_rare")
) |>
  select(-opp) |>
  relocate(data_type, study_name_common, study_name_rare, trait_common, trait_rare) # 294 direct matches

common_opp_code <- common_studies |>
  filter(grepl("Operation code", trait)) |>
  dplyr::select(data_type, study_name, trait) |>
  mutate(code = sub(".*: ", "", trait)) # 109

rare_opp_code <- rare_studies |>
  filter(grepl("Operation code", trait)) |>
  dplyr::select(data_type, study_name, trait, fieldID) |>
  mutate(code = sub("Operation code \\((.*)\\)$", "\\1", trait)) # 182

common_rare_opp_code <- inner_join(common_opp_code, rare_opp_code,
  by = c("data_type", "code"), suffix = c("_common", "_rare")
) |>
  select(-code) |>
  relocate(data_type, study_name_common, study_name_rare, trait_common, trait_rare) # 106 direct matches

# ----------- Remaining ---------- #
common_studies_remaining <- common_studies |>
  filter(!(study_name %in% c(
    common_rare_proteins$study_name_common,
    common_rare_ICD10$study_name_common,
    common_rare_med$study_name_common,
    common_rare_opp$study_name_common,
    common_rare_opp_code$study_name_common
  ))) # 1485 remaining

rare_studies_remaining <- rare_studies |>
  filter(!(study_name %in% c(
    rare_proteins$study_name,
    rare_ICD10$study_name,
    rare_med$study_name,
    rare_opp$study_name,
    rare_opp_code$study_name
  ))) # 7536 remaining

# ----------- Amalgamated diseases ---------- #
rare_union <- rare_studies |>
  filter(grepl("union", trait)) |>
  dplyr::select(data_type, study_name, trait, fieldID) |>
  mutate(union = sub("Chapter [A-Z]+ ", "", sub(" \\(union\\)", "", trait))) # 4505

common_union <- common_studies_remaining |>
  dplyr::select(data_type, study_name, trait) |>
  mutate(union = trait)

common_rare_union <- inner_join(common_union, rare_union,
  by = c("data_type", "union"), suffix = c("_common", "_rare")
) |>
  select(-union) |>
  relocate(data_type, study_name_common, study_name_rare, trait_common, trait_rare) # 56 direct matches

# ----------- Direct string matches ---------- #
common_rare_match <- inner_join(common_studies_remaining[, c(
  "data_type",
  "study_name",
  "trait"
)], rare_studies[, c("data_type", "study_name", "trait", "fieldID")],
by = c("data_type", "trait"), suffix = c("_common", "_rare")
) |>
  rename(trait_common = trait) |>
  mutate(trait_rare = trait_common) |>
  relocate(data_type, study_name_common, study_name_rare, trait_common, trait_rare, fieldID) # 51 direct matches

# ----------- Modified string matches ---------- #
common_rare_modmatch <- inner_join(
  common_studies_remaining[, c("data_type", "study_name", "trait")] |>
    mutate(traitmod = sub("\\/", " ", trait)),
  rare_studies_remaining[, c("data_type", "study_name", "trait", "fieldID")] |>
    mutate(traitmod = sub(" \\([0-9].*\\)$", "", sub(" - ", ": ", trait))),
  by = c("data_type", "traitmod"), suffix = c("_common", "_rare")
) |>
  relocate(data_type, study_name_common, study_name_rare, trait_common, trait_rare, fieldID) # 166 matches

# ----------- Remaining ---------- #
common_studies_remaining <- common_studies_remaining |>
  filter(!(study_name %in% c(
    common_rare_union$study_name_common,
    common_rare_match$study_name_common,
    common_rare_modmatch$study_name_common
  ))) # 1263 remaining

rare_studies_remaining <- rare_studies_remaining |>
  filter(!(study_name %in% c(
    common_rare_union$study_name_rare,
    common_rare_match$study_name_rare,
    common_rare_modmatch$study_name_rare
  ))) # 7326 remaining

# ----------- Field IDs from .json files ---------- #
notes <- data.table::fread(here("pre_steps/data/opengwas_notes.txt")) |>
  mutate(
    study_name = sub("_", "-", sub("UKB", "ukb", id)),
    fieldID = sub(":.*", "", note)
  ) |>
  filter(study_name %in% common_studies_remaining$study_name & !(fieldID == ""))

common_rare_fields <- inner_join(
  notes,
  rare_studies_remaining[, c("data_type", "study_name", "trait", "fieldID")],
  by = "fieldID", suffix = c("_common", "_rare")
) |>
  group_by(fieldID) |>
  filter(n() == 1) |>
  ungroup() |>
  as.data.frame() |>
  dplyr::select(study_name_common, study_name_rare, trait_common, trait_rare, fieldID) |>
  mutate(data_type = "phenotype", .before = "study_name_common") # 200 matches

# ----------- Modified string match 2 ---------- #
common_rare_modmatch2 <- inner_join(
  common_studies_remaining[, c("data_type", "study_name", "trait")],
  rare_studies_remaining[, c("data_type", "study_name", "trait", "fieldID")] |>
    filter(grepl("^[A-Z][0-9].*-[A-Z][0-9].*\\(union\\)", trait)) |>
    mutate(traitmod = sub(" \\(union\\)", "", sub("^[A-Z][0-9]+-[A-Z][0-9]+ ", "", trait))),
  by = c("data_type", "trait" = "traitmod"), suffix = c("_common", "_rare")
) |>
  rename(trait_common = trait) |>
  relocate(data_type, study_name_common, study_name_rare, trait_common, trait_rare, fieldID)

## Remove self reported and secondary codes
common_studies_remaining <- common_studies_remaining |>
  filter(!(study_name %in% c(
    common_rare_fields$study_name_common, common_rare_modmatch2$study_name_common
  )) &
    !(grepl("self-reported|secondary", trait))) # 724 remaining

rare_studies_remaining <- rare_studies_remaining |>
  filter(!(study_name %in% c(
    common_rare_fields$study_name_rare, common_rare_modmatch2$study_name_rare
  ))) # 7101

# ----------- Manual EFO matching for remainder --------- #
efos <- readr::read_csv(url("https://raw.githubusercontent.com/MRCIEU/finetune-efo/refs/heads/main/manual-mapping/trait-efo.csv")) # nolint: line_length_linter.

# Pull out EFO terms for remaining unmatched common studies
efo_terms <- efos |>
  filter(trait_id %in% common_studies_remaining$study_name & confidence > 4) |>
  pull(efo) |>
  unique()

# Pull high confidence terms including common and rare studies remaining unmatched
efos_remainder <- efos |>
  arrange(efo) |>
  filter(efo %in% efo_terms & confidence > 4) |>
  filter(trait_id %in% c(common_studies_remaining$study_name, rare_studies_remaining$study_name))

# Write out and manually match
# data.table::fwrite(efos_remainder, here("pre_steps/data/tmp_unmatched_efos.txt"))

manual <- data.table::fread(here("pre_steps/data/manual_common_rare_pairing.csv")) |>
  unique() |>
  relocate(study_name_common, study_name_rare, trait_common, trait_rare) |>
  mutate(data_type = "phenotype", .before = study_name_common) # 199

# Combine all successful pairs

out_match <- rbind(
  common_rare_proteins[, 1:5],
  common_rare_ICD10[, 1:5],
  common_rare_med[, 1:5],
  common_rare_opp[, 1:5],
  common_rare_opp_code[, 1:5],
  common_rare_union[, 1:5],
  common_rare_match[, 1:5],
  common_rare_modmatch[, 1:5],
  common_rare_fields[, 1:5],
  common_rare_modmatch2[, 1:5],
  manual[, 1:5]
)

data.table::fwrite(out_match, here("pipeline_steps/data/common_rare_ukbstudypairs.csv"))
