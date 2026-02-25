# Matching T-I pairs from Minikel et al. to GPMAP molecular traits and phenotypes

library(here)
library(dotenv)

renv::load(project = here("genotype-phenotype-map-analysis"))

readRenviron(here("genotype-phenotype-map-analysis/.env"))
data_dir <- Sys.getenv("DATA_DIR")

# Data from Minikel et al: https://zenodo.org/records/10783210

## Supplementary Table1 (target indication pairs: n=1255 with supporting genetic evidence, n=24458 unsupported)
ti_pairs <- data.table::fread(file.path(data_dir, "ericminikel-genetic_support/Supplement_ST01.csv"))

## Genetically supported/unsupported T-I pairs (from Figure 1A, ST04)
table(ti_pairs$combined_max_phase, ti_pairs$target_status)

## All significant genetic associations
assocs <- data.table::fread(file.path(data_dir, "ericminikel-genetic_support/data/assoc.tsv.gz"))

### ----- STEP 1: -----
# Derive table of all association study trait --> MeSH term (+ MeSH ID) mappings

term_tab <- assocs |>
  dplyr::select(original_trait, mesh_id, mesh_term) |>
  unique()

# data.table::fwrite(term_tab, file = file.path(data_dir,"target-indication_data/ericminike_traittomesh.csv"), quote = TRUE) # nolint: line_length_linter.

### ----- STEP 2: -----
# Derive table of ingested GPMAP study traits to match
# Retain common variant studies (colocalisation has not been run for rare exome variant studies)

gpmap_traits <- gpmapr::traits_api()[[1]] |>
  dplyr::filter(variant_type == "Common") |>
  dplyr::filter(duplicated(trait_name) == FALSE)

colnames(gpmap_traits) <- paste0(colnames(gpmap_traits), "_gpmap")

### ----- STEP 3: -----
# Map GPMAP traits to the trait terms in term_tab so the matched MeshID can be identified
# The assumption here is that the phenotypic associations sourced in Minikel et al. are the more exhaustive set, and should # nolint: line_length_linter.
# contain the EBI and UKB studies ingested for the GPMAP

# Pass a .csv containing the original traits from Minikel et al and GPMAP traits
forAI <- cbind(
  gpmap_trait = gpmap_traits$trait_name_gpmap[seq_len(nrow(term_tab))],
  original_trait = term_tab$original_trait
)
# data.table::fwrite(forAI, file = file.path(data_dir,"target-indication_data/ericminike_gpmap_traitstomatch.csv"), quote = TRUE) # nolint: line_length_linter.

# ChatGPT5 was run with the following prompt:
# nolint start: line_length_linter.
# In the .csv file attached, the first column ("gpmap_trait") contain the trait names from a set of genetic studies. The second column
# ("original_trait") contains an independent set of trait names from a more comprehensive second set of studies. For each "gpmap_trait",
# can you select the "original_trait" you believe best matches from the set provided, and place these in column 3. Rate your confidence
# in the match (between 0 and 1) in column 4. For example, a direct string match would receive a score of 1, and a more partial match
# (e.g "diabetes mellitus" vs "type 2 diabetes") may receive a score of 0.5. You can recycle rows in original_trait if one serves as the
# best match for multiple gpmap_traits. There are more original_traits than gpmap_traits, so not all will be used. Please return the output .csv file
# nolint end

# Output: gpmap_matched_traits.csv

matched <- data.table::fread(file = file.path(data_dir, "target-indication_data/gpmap_matched_traits.csv")) |>
  dplyr::select(-original_trait) |>
  dplyr::arrange(desc(confidence)) |>
  dplyr::filter(nchar(gpmap_trait) > 0) |>
  dplyr::left_join(term_tab, by = c("best_match" = "original_trait"))

# Join to gpmap_traits
gpmap_traits <- gpmap_traits |> dplyr::left_join(matched, by = c("trait_name_gpmap" = "gpmap_trait"))

### ----- STEP 4: -----
# Retain only the GPMAP traits which can be matched to a drug-target indication in ti_pairs (based on MeSH term similarity of >=0.8) # nolint: line_length_linter.

# Read in MeSH ID similarity matrix from Minikel et al.
mesh_match <- data.table::fread(file.path(data_dir, "ericminikel-genetic_support/data/sim.tsv.gz"))
mesh_terms <- data.frame(
  mesh_id = c(ti_pairs$indication_mesh_id, term_tab$mesh_id),
  mesh_term = c(ti_pairs$indication_mesh_term, term_tab$mesh_term)
) |>
  unique()
mesh_match_wterms <- mesh_match |>
  dplyr::left_join(mesh_terms, by = c("meshcode_a" = "mesh_id")) |>
  dplyr::left_join(mesh_terms, by = c("meshcode_b" = "mesh_id"), suffix = c("_a", "_b")) |>
  na.exclude()

# Indication MeSH IDs
i_mesh <- unique(ti_pairs$indication_mesh_id) # 769 IDs
i_term_tab <- unique(ti_pairs[, c("indication_mesh_id", "indication_mesh_term")])

# Table of indicator MeSH terms and matched MeSH terms to pair by
mesh_match_highconf <- mesh_match |>
  dplyr::filter(
    meshcode_a %in% i_mesh | meshcode_b %in% i_mesh,
    comb_norm >= 0.8
  )

# Which MeSH code is the indicator code, and which is the matched code
i_idx <- apply(mesh_match_highconf, 1, function(x) {
  return(which(x %in% i_mesh))
})
i_idx_both <- which(lapply(i_idx, length) == 2) # Both codes are in the indicator set
i_idx <- sapply(i_idx, function(x) {
  return(x[1])
})

mesh_match_highconf <- mesh_match_highconf |>
  dplyr::mutate(
    indication_mesh_id = ifelse(i_idx == 1, meshcode_a, meshcode_b),
    matched_mesh_id = ifelse(i_idx == 1, meshcode_b, meshcode_a)
  ) |>
  dplyr::select(-c(meshcode_a, meshcode_b)) |>
  dplyr::relocate(comb_norm, .after = matched_mesh_id) |>
  rbind(data.frame(
    indication_mesh_id = mesh_match_highconf$meshcode_b[i_idx_both],
    matched_mesh_id = mesh_match_highconf$meshcode_a[i_idx_both],
    comb_norm = mesh_match_highconf$comb_norm[i_idx_both]
  )) |> # Add inverse for rows with two indicator Mesh IDs
  dplyr::arrange(indication_mesh_id, comb_norm) |>
  unique() |>
  dplyr::left_join(i_term_tab)

# Full set of MeSH IDs which directly or closely match indication MeSH IDs
mesh_set <- unique(c(mesh_match_highconf$indication_mesh_id, mesh_match_highconf$matched_mesh_id)) # 1622

# Extract these from the GPMAP studies
gpmap_matched <- gpmap_traits |> dplyr::filter(mesh_id %in% mesh_set) # 2171 studies

# Manually check the pairing of GPMAP traits to original traits conducted above (based on AI reported best_match and confidence) # nolint: line_length_linter.
# and exclude those where the GPMAP trait is not the same or sufficiently comparable to the original trait, or
# sufficiently captured by the MeSH term assigned to the original trait in Minikel et al.
# (see gpmap_matched_traits_filtered.xlsx for excluded traits)

# data.table::fwrite(gpmap_matched, file.path(data_dir,"target-indication_data/gpmap_matched_traits_tofilter.csv"))

gpmap_matched_filtered <- data.table::fread(file.path(
  data_dir,
  "target-indication_data/gpmap_matched_traits_filtered_corrected.csv"
))
# 1423 studies with high confidence mapping to MeSH terms for disease indications (343 unique MeSH terms)

# Match GPMAP studies to indicators
# If a study MeSH term matches to multiple indicator MeSH terms keep the match to both

gpmap_matched_filtered <- gpmap_matched_filtered |>
  dplyr::left_join(mesh_match_highconf, by = c("mesh_id" = "matched_mesh_id"), relationship = "many-to-many")

### ----- STEP 5: -----
# Identify GPMAP molecular studies available for drug targets, including:
# GWAS of target gene expression (identifying cis-eQTLs)
# GWAS of target isoform expression (identifying cis-sQTLs)
# GWAS of target protein expression (identifying cis and trans pQTLs)
# GWAS of cpg methylation across sites proximal to the target (identifying cis-mQTLs)

targets <- unique(ti_pairs$target) # 2480 drug targets

gpmap_molecular <- data.table::fread(file.path(data_dir, "studies_processed.tsv.gz")) |>
  dplyr::filter(variant_type == "common" & gene %in% targets)

matched_targets <- unique(gpmap_molecular$gene)
length(matched_targets) # 2350 drug targets with genetic associations tested across molecular assays

# What % of drug targets have data from each molecular assays available (in at least one tissue)?:
available_assays <- gpmap_molecular |>
  dplyr::group_by(gene, data_type) |>
  dplyr::summarise(n_studies = dplyr::n())
table(available_assays$data_type) / length(matched_targets) * 100 # 97.9% eQTL, 91.3% sQTL, 87.5% mQTL, 38.4% pQTL

### ----- STEP 6: -----

# Bind GPMAP association studies to ti_pairs (keeping pairs for which both the indication and the target are included in GPMAP studies) # nolint: line_length_linter.
ti_pairs_gpmap <- dplyr::left_join(
  ti_pairs[, c(
    "ti_uid", "target",
    "indication_mesh_id", "indication_mesh_term", "historical_max_phase",
    "active_max_phase", "combined_max_phase",
    "succ_p_1", "succ_1_2", "succ_2_3", "succ_3_a", "orphan", "year_launch", "target_status"
  )],
  gpmap_matched_filtered,
  by = c("indication_mesh_id", "indication_mesh_term"), relationship = "many-to-many"
) |>
  dplyr::select(-poor_match) |>
  dplyr::relocate("indication_association_similarity" = comb_norm, .after = indication_mesh_term) |>
  dplyr::filter(!is.na(data_type_gpmap) & target %in% matched_targets) # Keep t-i pairs with targets in GPMAP only

length(unique(ti_pairs_gpmap$ti_uid)) # 18932/25713 (74%) unique T-I pairs tested in Minikel et al. have genetic evidence in GPMAP # nolint: line_length_linter.

## Check which indications are not matching
missing <- i_term_tab[!(i_term_tab$indication_mesh_id %in% gpmap_matched_filtered$indication_mesh_id)]

# Write out final table of gpmap studies paired with phamaceutical indications
# data.table::fwrite(ti_pairs_gpmap, file.path(data_dir, "target-indication_data/target-indicationpairs_gpmapevidence.tsv"), sep = "\t") # nolint: line_length_linter.
