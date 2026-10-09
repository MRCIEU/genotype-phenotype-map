#!/usr/bin/env Rscript
# Deduplicate common phenotypic traits for EBMF trait decomposition.
#
# Rationale: the GPMap phenotype catalogue holds the same trait several times
# under slightly different names (Firth and SPA runs of one UKB GWAS,
# "(UKB data field N)" re-analyses of UKB fields, "Standing height" vs
# "Height", "Eosinophill count" vs "Eosinophil counts", ...). Near-identical
# rows in the traits x SNPs matrix soak up their own EBMF factors, so obvious
# duplicates are collapsed to one representative trait per group.
#
# Method: each trait name is reduced to a dedup key by an ordered set of named
# rules (`build_name_rules()`), followed by a curated spelling map and synonym
# map. Traits sharing (key, category) form a group, so continuous and
# categorical traits are never merged. The representative of each group is the
# trait with the most coloc groups, then the largest sample size, then the
# most study extractions, then the trait id string.
#
# Input: the raw results files (studies_processed, traits_processed,
# study_extractions, coloc_clustered_results) in `--results_dir`. Only traits
# with at least one study extraction are included, matching the traits table
# built by create_db_from_results.R, which loads this output into the
# `trait_duplicates` table.
#
# Output: one row per in-scope trait (singletons included), so `keep == TRUE`
# is the complete list of traits to use. `match_rules` records which rules
# changed each trait name, so merges made by a given rule (e.g. left_right,
# main_secondary) can be filtered out downstream.
#
# Usage (from scripts/): Rscript deduplicate_traits.R [--results_dir <dir>]
# After ingesting new traits: re-run, review the multi-member groups and the
# high-coverage singletons printed at the end, and extend `spelling_map` /
# `synonym_map` below for any obvious duplicates the rules miss.

source("../pipeline_steps/constants.R")
suppressPackageStartupMessages(library(dplyr))

parser <- argparser::arg_parser("Deduplicate common phenotypic traits")
parser <- argparser::add_argument(
  parser,
  "--results_dir",
  help = "Results directory to read raw results files from (default: latest results directory)",
  type = "character",
  default = as.character(latest_results_dir)
)
args <- argparser::parse_args(parser)

# ---- Settings ---------------------------------------------------------------
output_path <- trait_deduplication_file
collapse_left_right <- TRUE       # merge left and right side measurements
collapse_main_secondary <- TRUE   # merge main and secondary ICD10/OPCS codes
n_review_singletons <- 50L        # high-coverage singletons printed for review

# Whole-word replacements applied to normalised names
spelling_map <- c(
  eosinophill = "eosinophil",
  neutrophill = "neutrophil",
  basophill = "basophil",
  hematocrit = "haematocrit",
  hemoglobin = "haemoglobin",
  counts = "count",
  triglycerides = "triglyceride"
)

# Whole-key replacements applied last; names and values are normalised keys
wbc_types <- c("eosinophil", "basophil", "lymphocyte", "monocyte", "neutrophil")
wbc_count_types <- c("eosinophil", "basophil", "lymphocyte", "monocyte")
synonym_map <- c(
  "clotting disorder excessive bleeding" = "clotting disorder or excessive bleeding",
  "standing height" = "height",
  "haemoglobin concentration" = "haemoglobin",
  "haematocrit percentage" = "haematocrit",
  "platelet crit" = "plateletcrit",
  "mean spheric corpuscular volume" = "mean sphered cell volume",
  "red blood cell distribution width" = "red cell distribution width",
  "immature fraction of reticulocytes" = "immature reticulocyte fraction",
  "reticulocyte percentage" = "reticulocyte fraction of red cells",
  "high light scatter reticulocyte percentage" =
    "high light scatter reticulocyte percentage of red cells",
  "glycated haemoglobin hba1c" = "haemoglobin a1c",
  "glycated haemoglobin" = "haemoglobin a1c",
  "cholesterol" = "total cholesterol",
  "cholesterol total" = "total cholesterol",
  "high density lipoprotein cholesterol" = "hdl cholesterol",
  "low density lipoprotein cholesterol" = "ldl cholesterol",
  "direct low density lipoprotein" = "ldl cholesterol",
  "direct low density lipoprotein cholesterol" = "ldl cholesterol",
  "apolipoprotein a" = "apolipoprotein a1",
  "apolipoprotein a i" = "apolipoprotein a1",
  "igf 1" = "insulin like growth factor 1",
  "gamma glutamyltransferase" = "gamma glutamyl transferase",
  "gamma glutamyl transpeptidase" = "gamma glutamyl transferase",
  "fev1" = "forced expiratory volume in 1 second",
  "lung function fvc" = "forced vital capacity",
  stats::setNames(
    paste(wbc_count_types, "count"), paste("white blood cell count", wbc_count_types)
  ),
  stats::setNames(
    paste(wbc_types, "percentage of white cells"), paste(wbc_types, "percentage")
  )
)

# ---- Name normalisation rules ----------------------------------------------
# Each rule maps a character vector of names to normalised names. Rules run in
# order; a rule is recorded in `match_rules` for every name it changes.
build_name_rules <- function(collapse_left_right, collapse_main_secondary,
                             spelling_map, synonym_map) {
  rules <- list(
    lowercase = function(x) {
      return(tolower(gsub("\u00a0", " ", x)))
    },
    firth_spa = function(x) {
      return(gsub("\\s*\\((firth|spa) correction\\)", "", x))
    },
    ukb_field = function(x) {
      return(gsub("\\s*\\(ukb data field [0-9_]+\\)", "", x))
    },
    main_secondary = function(x) {
      return(sub(
        "^(diagnoses|operative procedures) - (main|secondary) (icd10|opcs):",
        "\\1 - \\3:", x
      ))
    },
    # Only the trailing side qualifier is dropped, so mid-name laterality
    # ("left ventricular failure", "left side of heart") is left untouched
    left_right = function(x) {
      x <- gsub("\\s*\\((left|right)\\)", "", x)
      x <- sub(":\\s*(left|right) eye$", ": single eye", x)
      return(sub("\\s+(left|right)$", "", x))
    },
    gloss = function(x) {
      x <- gsub("(blood cell|platelet) \\(?(erythrocyte|leukocyte|thrombocyte)\\)?", "\\1", x)
      return(gsub("\\s*\\(fev1\\)", "", x))
    },
    # "Body mass index (BMI)", "Peak expiratory flow (PEF)": drop the first
    # parenthetical when it only restates the initials of the words before it
    acronym = function(x) {
      parts <- regmatches(x, regexec("^(.*?)\\s*\\(([a-z]{2,6})\\)(.*)$", x))
      return(vapply(seq_along(x), function(i) {
        m <- parts[[i]]
        if (length(m) == 0) return(x[i])
        words <- strsplit(m[2], "[^a-z0-9]+")[[1]]
        words <- words[words != ""]
        n <- nchar(m[3])
        if (length(words) < n) return(x[i])
        initials <- paste(substr(utils::tail(words, n), 1, 1), collapse = "")
        if (initials != m[3]) return(x[i])
        return(paste0(m[2], m[4]))
      }, character(1)))
    },
    # "+"/"-" are kept as pos/neg so cell populations such as CD14+ CD16- and
    # CD14- CD16+ stay distinct. A hyphen is a negative marker when it trails
    # a token or directly follows a CD marker (CD4-CD8-); otherwise it joins
    # words (C-reactive, Interleukin-2) and becomes a space.
    punctuation = function(x) {
      x <- gsub("\\s+-\\s+", " ", x)
      x <- gsub("\\b(cd[0-9]+[a-z]*)-", "\\1 neg ", x, perl = TRUE)
      x <- gsub("-(?=\\s|$|[,;:/)+-])", " neg ", x, perl = TRUE)
      x <- gsub("-", " ", x)
      x <- gsub("\\+", " pos ", x)
      x <- gsub("%", " pct ", x)
      x <- gsub("[^a-z0-9]+", " ", x)
      return(trimws(gsub("\\s+", " ", x)))
    },
    # Specimen prefixes are only dropped from "... levels" names, so traits
    # such as "Blood pressure ..." or "Blood clot ..." keep their prefix
    affix = function(x) {
      x <- sub("^(serum|plasma|circulating|blood) (.+ levels?)$", "\\2", x)
      return(sub(" (levels?|measurement)$", "", x))
    },
    spelling = function(x) {
      for (from in names(spelling_map)) {
        x <- gsub(paste0("\\b", from, "\\b"), spelling_map[[from]], x)
      }
      return(x)
    },
    synonym = function(x) {
      hit <- x %in% names(synonym_map)
      x[hit] <- unname(synonym_map[x[hit]])
      return(x)
    }
  )
  if (!collapse_left_right) rules$left_right <- NULL
  if (!collapse_main_secondary) rules$main_secondary <- NULL
  return(rules)
}

normalise_trait_names <- function(trait_names, rules) {
  key <- trait_names
  applied <- vector("list", length(trait_names))
  for (rule in names(rules)) {
    new_key <- rules[[rule]](key)
    changed <- new_key != key
    # lowercase and punctuation fire on almost every name; only record the
    # rules that say something about why two names matched
    if (!rule %in% c("lowercase", "punctuation")) {
      applied[changed] <- lapply(applied[changed], c, rule)
    }
    key <- new_key
  }
  match_rules <- vapply(applied, paste, character(1), collapse = ";")
  return(list(key = key, match_rules = match_rules))
}

# ---- Fetch traits and build dedup groups -----------------------------------
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)

study_extractions <- vroom::vroom(
  file.path(args$results_dir, "study_extractions.tsv.gz"),
  col_select = c(study, unique_study_id, ignore),
  show_col_types = FALSE
) |>
  filter(is.na(ignore) | ignore == FALSE)

studies <- vroom::vroom(
  file.path(args$results_dir, "studies_processed.tsv.gz"),
  col_types = studies_processed_column_types,
  show_col_types = FALSE
) |>
  filter(
    data_type == data_types$phenotype,
    variant_type == variant_types$common,
    study_name %in% study_extractions$study
  ) |>
  select(study_name, trait, sample_size, category)

coloc_groups <- vroom::vroom(
  file.path(args$results_dir, "coloc_clustered_results.tsv.gz"),
  col_select = c(unique_study_id, coloc_group_id),
  show_col_types = FALSE
)

trait_counts <- study_extractions |>
  inner_join(select(studies, study_name, trait), by = c("study" = "study_name")) |>
  left_join(coloc_groups, by = "unique_study_id") |>
  group_by(trait) |>
  summarise(
    num_study_extractions = n_distinct(unique_study_id),
    num_coloc_groups = n_distinct(coloc_group_id, na.rm = TRUE),
    .groups = "drop"
  )

trait_studies <- studies |>
  group_by(trait) |>
  summarise(sample_size = max(sample_size), category = first(category), .groups = "drop")

traits <- vroom::vroom(
  file.path(args$results_dir, "traits_processed.tsv.gz"),
  show_col_types = FALSE
) |>
  select(study_name, trait, category) |>
  rename(trait_name = trait, trait = study_name, trait_category = category) |>
  inner_join(trait_studies, by = "trait") |>
  inner_join(trait_counts, by = "trait") |>
  select(
    trait, trait_name, trait_category, category, sample_size,
    num_study_extractions, num_coloc_groups
  )

rules <- build_name_rules(
  collapse_left_right, collapse_main_secondary, spelling_map, synonym_map
)
normalised <- normalise_trait_names(traits$trait_name, rules)
traits$dedup_key <- paste(normalised$key, traits$category, sep = " | ")
traits$match_rules <- normalised$match_rules

ranked <- traits |>
  group_by(dedup_key) |>
  arrange(
    desc(num_coloc_groups), desc(sample_size), desc(num_study_extractions), trait,
    .by_group = TRUE
  ) |>
  mutate(
    group_size = n(),
    group_rank = row_number(),
    keep = group_rank == 1L,
    representative_trait = first(trait),
    representative_trait_name = first(trait_name)
  ) |>
  ungroup()

group_order <- ranked |>
  filter(keep) |>
  arrange(desc(num_coloc_groups), trait) |>
  transmute(dedup_key, dedup_group = row_number())

dedup <- ranked |>
  inner_join(group_order, by = "dedup_key") |>
  arrange(dedup_group, group_rank) |>
  select(
    trait, trait_name, trait_category, category, sample_size,
    num_study_extractions, num_coloc_groups, dedup_group, dedup_key,
    group_size, group_rank, keep, representative_trait,
    representative_trait_name, match_rules
  )

# ---- Checks -----------------------------------------------------------------
group_checks <- dedup |>
  group_by(dedup_group) |>
  summarise(
    n_keep = sum(keep),
    rep_in_group = all(representative_trait %in% trait),
    n_category = n_distinct(category),
    .groups = "drop"
  )
stopifnot(
  !anyDuplicated(dedup$trait),
  nrow(dedup) == nrow(traits),
  all(group_checks$n_keep == 1L),
  all(group_checks$rep_in_group),
  all(group_checks$n_category == 1L)
)

readr::write_tsv(dedup, output_path, na = "")

# ---- Summary ----------------------------------------------------------------
multi <- filter(dedup, group_size > 1L)
message(
  "Traits in scope: ", nrow(dedup),
  " | duplicate groups: ", n_distinct(multi$dedup_group),
  " | traits in duplicate groups: ", nrow(multi),
  " | traits dropped: ", sum(!dedup$keep),
  " | traits kept: ", sum(dedup$keep)
)

pre_synonym_keys <- normalise_trait_names(
  traits$trait_name, rules[names(rules) != "synonym"]
)$key
unused_synonyms <- setdiff(names(synonym_map), pre_synonym_keys)
if (length(unused_synonyms) > 0) {
  message("\nSynonym entries matching no trait: ", paste(unused_synonyms, collapse = ", "))
}

# Same normalised name but different `category`: kept apart by the category
# guard, usually because of inconsistent trait metadata
category_split <- dedup |>
  mutate(name_key = sub(" \\| [^|]*$", "", dedup_key)) |>
  group_by(name_key) |>
  filter(n_distinct(category) > 1L) |>
  ungroup()
if (nrow(category_split) > 0) {
  message("\nSame name key split by category (not merged):")
  category_split |>
    arrange(name_key, category, group_rank) |>
    select(name_key, trait, trait_name, category, num_coloc_groups) |>
    as.data.frame() |>
    print(row.names = FALSE)
}

message("\nTraits merged per rule (traits in duplicate groups only):")
rule_counts <- table(unlist(strsplit(multi$match_rules[multi$match_rules != ""], ";")))
print(sort(rule_counts, decreasing = TRUE))

message("\nLargest duplicate groups:")
multi |>
  distinct(dedup_group, group_size, representative_trait, representative_trait_name) |>
  arrange(desc(group_size), dedup_group) |>
  head(20) |>
  as.data.frame() |>
  print(row.names = FALSE)

message("\nHighest-coverage singletons (check for missed synonyms):")
dedup |>
  filter(group_size == 1L) |>
  arrange(desc(num_coloc_groups)) |>
  select(trait, trait_name, category, sample_size, num_coloc_groups) |>
  head(n_review_singletons) |>
  as.data.frame() |>
  print(row.names = FALSE)

message("\nWrote ", output_path)
