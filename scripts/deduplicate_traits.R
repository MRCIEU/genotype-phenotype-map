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
# Proxy groups: beyond hard duplicates, many traits measure the same liability
# through another reporting route (self-reported illness, doctor diagnosed,
# ICD10 diagnosis, medication, family history, ...). For each duplicate group
# the representative's name has its reporting wrapper stripped
# (`proxy_wrappers`, `proxy_icd10_pattern`), is normalised with the same rules
# as above, and is mapped to a concept via `proxy_icd10_map` /
# `proxy_concept_map` (else the key itself), with `proxy_trait_overrides` for
# single traits. Concepts covering more than one duplicate group become a
# `proxy_group`; duplicates share a proxy group. Body-size measures are also
# grouped (`height`, `adiposity`) so a user can exclude several proxy groups at
# once, e.g. `body mass index` and `adiposity` for a BMI anchor.
#
# Output: one row per in-scope trait (singletons included). `duplicate_group`
# and `is_representative` give the duplicate groups (`is_representative ==
# TRUE` is the complete list of traits to use); `proxy_group` (and its display
# name `proxy_group_name`) is NA when a trait has no proxies. Downstream, an anchor trait keeps itself and drops
# every other trait sharing its duplicate_group or proxy_group. `match_rules`
# and `proxy_rules` record which rules produced each match, so merges made by a
# given rule (e.g. left_right, main_secondary) can be filtered out downstream.
#
# Usage (from scripts/): Rscript deduplicate_traits.R [--results_dir <dir>]
# After ingesting new traits: re-run, review the multi-member groups, the
# high-coverage singletons, the proxy groups and the missed/suspicious proxy
# lists printed at the end, then extend `spelling_map` / `synonym_map` for
# duplicates and `proxy_concept_map` / `proxy_icd10_map` /
# `proxy_trait_overrides` for proxies.

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
n_review_proxy_groups <- 40L      # largest proxy groups printed for review
n_review_proxy_pairs <- 100L      # missed proxy candidates printed for review
proxy_review_min_shared <- 10L    # min shared coloc groups for a missed proxy candidate
proxy_review_min_overlap <- 0.5   # min shared / smaller trait's coloc groups for a candidate
proxy_review_max_member_overlap <- 0.1  # proxy members below this overlap are flagged

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

# ---- Proxy settings ---------------------------------------------------------
# Reporting wrappers stripped from lowercased names to expose the underlying
# condition. Each pattern has one capture group; the first match wins.
proxy_wrappers <- c(
  "^(?:non-cancer illness code|cancer code), self-reported: (.+)$",
  "^.+ diagnosed by (?:doctor|a professional): (.+)$",
  "^illnesses of (?:adopted )?(?:mother|father|siblings): (.+)$",
  "^age (.+) diagnosed(?: by doctor)?$",
  "^(?:doctor diagnosed|self-reported) (.+)$",
  "^(.+) diagnosed by doctor$",
  "^treatment/medication code: (.+)$",
  "^medication use \\((.+)\\)$",
  "^medication for [^:]+: (.+)$",
  "^eye problems/disorders: (.+)$",
  "^operation code: (.+)$"
)
# ICD10 diagnosis names: group 1 is the 3 character code, group 2 the description
proxy_icd10_pattern <- "^(?:diagnoses - (?:main|secondary) icd10|type of cancer: icd10): ([a-z][0-9]{2})[0-9.]* (.+)$"

# Wrapped answers that say nothing about a condition never form a proxy group
proxy_non_answers <- c(
  "none of the above", "none of the above group 1", "none of the above group 2",
  "do not know", "prefer not to answer", "unclassifiable"
)

# Normalised keys (after wrapper stripping and the name rules above) mapped to
# a shared proxy concept. Each trait has one concept, so generic traits go to
# the dominant concept (e.g. unspecified diabetes -> type 2 diabetes).
proxy_concept_map <- c(
  "diabetes" = "type 2 diabetes",
  "drugs used in diabetes" = "type 2 diabetes",
  "metformin" = "type 2 diabetes",
  "gliclazide" = "type 2 diabetes",
  "diabetes related eye disease" = "type 2 diabetes",
  # Insulin use colocalises with type 1, not type 2, diabetes
  "insulin product" = "type 1 diabetes",
  "started insulin within one year diagnosis of diabetes" = "type 1 diabetes",
  "obesity" = "body mass index",
  # Body-size measures are not reporting routes, but most move with BMI (or
  # height) at almost every locus, so they are grouped to be excluded together
  "sitting height" = "height",
  "weight" = "adiposity",
  "weight manual entry" = "adiposity",
  "childhood body mass index" = "adiposity",
  "comparative body size at age 10" = "adiposity",
  "waist circumference" = "adiposity",
  "hip circumference" = "adiposity",
  "body fat percentage" = "adiposity",
  "whole body fat mass" = "adiposity",
  "trunk fat mass" = "adiposity",
  "arm fat mass" = "adiposity",
  "leg fat mass" = "adiposity",
  "trunk fat percentage" = "adiposity",
  "arm fat percentage" = "adiposity",
  "leg fat percentage" = "adiposity",
  "whole body fat free mass" = "adiposity",
  "trunk fat free mass" = "adiposity",
  "arm fat free mass" = "adiposity",
  "leg fat free mass" = "adiposity",
  "trunk predicted mass" = "adiposity",
  "arm predicted mass" = "adiposity",
  "leg predicted mass" = "adiposity",
  "whole body water mass" = "adiposity",
  "basal metabolic rate" = "adiposity",
  "impedance of whole body" = "adiposity",
  "impedance of arm" = "adiposity",
  "impedance of leg" = "adiposity",
  "ankle spacing width" = "adiposity",
  "high blood pressure" = "hypertension",
  "essential primary hypertension" = "hypertension",
  "blood pressure medication" = "hypertension",
  "antihypertensives" = "hypertension",
  "agents acting on the renin angiotensin system" = "hypertension",
  "calcium channel blockers" = "hypertension",
  "beta blocking agents" = "hypertension",
  "diuretics" = "hypertension",
  "bendroflumethiazide" = "hypertension",
  "amlodipine" = "hypertension",
  "atenolol" = "hypertension",
  "lisinopril" = "hypertension",
  "ramipril" = "hypertension",
  "enalapril" = "hypertension",
  "adalat 5mg capsule" = "hypertension",
  "nifedipine" = "hypertension",
  "lercanidipine" = "hypertension",
  "doxazosin" = "hypertension",
  "indapamide" = "hypertension",
  "bisoprolol" = "hypertension",
  "high cholesterol" = "hypercholesterolaemia",
  "pure hypercholesterolaemia" = "hypercholesterolaemia",
  "cholesterol lowering medication" = "hypercholesterolaemia",
  "cholesterol lowering medication use" = "hypercholesterolaemia",
  "hmg coa reductase inhibitors" = "hypercholesterolaemia",
  "simvastatin" = "hypercholesterolaemia",
  "atorvastatin" = "hypercholesterolaemia",
  "lipitor 10mg tablet" = "hypercholesterolaemia",
  "rosuvastatin" = "hypercholesterolaemia",
  "pravastatin" = "hypercholesterolaemia",
  "ezetimibe" = "hypercholesterolaemia",
  "familial combined hyperlipidemia defined by brunzell criteria" = "hypercholesterolaemia",
  "familial combined hyperlipidemia defined by consensus criteria" = "hypercholesterolaemia",
  "familial combined hyperlipidemia defined by dutch criteria" = "hypercholesterolaemia",
  "familial combined hyperlipidemia defined by goldstein criteria" = "hypercholesterolaemia",
  "familial combined hyperlipidemia defined by mexico criteria" = "hypercholesterolaemia",
  "heart attack" = "coronary artery disease",
  "heart attack myocardial infarction" = "coronary artery disease",
  "myocardial infarction" = "coronary artery disease",
  "coronary heart disease" = "coronary artery disease",
  "coronary atherosclerosis" = "coronary artery disease",
  "angina" = "coronary artery disease",
  "angina pectoris" = "coronary artery disease",
  "heart disease" = "coronary artery disease",
  "coronary artery bypass grafts" = "coronary artery disease",
  "coronary angioplasty ptca pos neg stent" = "coronary artery disease",
  "vasodilators used in cardiac diseases" = "coronary artery disease",
  "nicorandil" = "coronary artery disease",
  "gtn 400micrograms spray" = "coronary artery disease",
  "glyceryl trinitrate" = "coronary artery disease",
  "isosorbide mononitrate" = "coronary artery disease",
  "ischaemic heart disease wide definition" = "coronary artery disease",
  "digoxin" = "atrial fibrillation",
  "adrenergics inhalants" = "asthma",
  "salbutamol" = "asthma",
  "salmeterol" = "asthma",
  "beclometasone" = "asthma",
  "ventolin 100micrograms inhaler" = "asthma",
  "salamol 100micrograms inhaler" = "asthma",
  "bricanyl 250mcg inhaler" = "asthma",
  "serevent 25mcg inhaler" = "asthma",
  "seretide 50 evohaler" = "asthma",
  "becotide 50 inhaler" = "asthma",
  "beclazone 50 inhaler" = "asthma",
  "qvar 50 inhaler" = "asthma",
  "symbicort 100 6 turbohaler" = "asthma",
  "thyroid preparations" = "hypothyroidism",
  "levothyroxine sodium" = "hypothyroidism",
  "thyroxine product" = "hypothyroidism",
  "thyroxine sodium" = "hypothyroidism",
  "hypothyroidism myxoooedema" = "hypothyroidism",
  "hypothyroidism or myxoooedema" = "hypothyroidism",
  "copd chronic obstructive pulmonary disease" = "chronic obstructive pulmonary disease",
  "emphysema chronic bronchitis" = "chronic obstructive pulmonary disease",
  "chronic bronchitis emphysema" = "chronic obstructive pulmonary disease",
  "chronic bronchitis" = "chronic obstructive pulmonary disease",
  "tiotropium" = "chronic obstructive pulmonary disease",
  "alzheimer s disease dementia" = "alzheimer s disease",
  "any parental history of alzheimer s disease" = "alzheimer s disease",
  "alzheimer s disease or family history of alzheimer s disease" = "alzheimer s disease",
  "severe depression" = "depression",
  "depression broad" = "depression",
  "major depressive disorder" = "depression",
  "major depressive disorder icd 10 coded" = "depression",
  "major depressive disorder probable" = "depression",
  "paroxetine" = "depression",
  "mirtazapine" = "depression",
  "venlafaxine" = "depression",
  "dosulepin" = "depression",
  "hayfever allergic rhinitis or eczema" = "allergic disease",
  "hayfever or allergic rhinitis" = "allergic disease",
  "hayfever allergic rhinitis" = "allergic disease",
  "hay fever rhinitis or eczema" = "allergic disease",
  "allergic rhinitis" = "allergic disease",
  "allergic disease asthma hay fever or eczema" = "allergic disease",
  "eczema" = "allergic disease",
  "antihistamines for systemic use" = "allergic disease",
  "flixonase 50micrograms aqueous nasal spray" = "allergic disease",
  "omeprazole" = "gastroesophageal reflux disease",
  "omeprazole e g zanprol" = "gastroesophageal reflux disease",
  "esomeprazole" = "gastroesophageal reflux disease",
  "drugs for peptic ulcer and gastro oesophageal reflux disease" = "gastroesophageal reflux disease",
  "cholelithiasis" = "cholelithiasis gall stones",
  "cholecystitis" = "cholelithiasis gall stones",
  "cholecystectomy gall bladder removal" = "cholelithiasis gall stones",
  "allopurinol" = "gout",
  "sumatriptan" = "migraine",
  "antimigraine preparations" = "migraine",
  "xalatan 0 005 pct eye drops" = "glaucoma",
  "antiglaucoma preparations and miotics" = "glaucoma",
  "alendronate sodium" = "osteoporosis",
  "drugs affecting bone structure and mineralization" = "osteoporosis",
  "alfuzosin" = "hyperplasia of prostate",
  "tamsulosin" = "hyperplasia of prostate",
  "finasteride" = "hyperplasia of prostate",
  "malignant neoplasm of prostate" = "prostate cancer",
  "cancer of male genital" = "prostate cancer",
  "malignant neoplasm of male genital organs" = "prostate cancer",
  "malignant neoplasm of rectum" = "bowel cancer",
  "colon cancer sigmoid cancer" = "bowel cancer",
  "colourectal cancer" = "bowel cancer",
  "any parental history of lung cancer" = "lung cancer",
  "uterine endometrial cancer" = "endometrial cancer",
  "malignant neoplasm of skin" = "skin cancer",
  "basal cell carcinoma" = "skin cancer",
  "cancer of skin" = "skin cancer",
  "other malignant neoplasms of skin" = "skin cancer",
  "malignant melanoma" = "malignant melanoma of skin",
  "cerebral infarction" = "stroke",
  "retinal detachments and breaks" = "retinal detachment",
  "arthrosis" = "osteoarthritis",
  "polyarthrosis" = "osteoarthritis",
  "low back pain" = "back pain",
  "back problem" = "back pain"
)

# ICD10 code prefixes mapped to a proxy concept; these win over the description
proxy_icd10_map <- c(
  E10 = "type 1 diabetes",
  E11 = "type 2 diabetes",
  E14 = "type 2 diabetes",
  E66 = "body mass index",
  E78 = "hypercholesterolaemia",
  I10 = "hypertension",
  I20 = "coronary artery disease",
  I21 = "coronary artery disease",
  I25 = "coronary artery disease",
  I83 = "varicose veins",
  J44 = "chronic obstructive pulmonary disease",
  J45 = "asthma",
  F32 = "depression",
  F41 = "anxiety disorders",
  K21 = "gastroesophageal reflux disease",
  K80 = "cholelithiasis gall stones",
  K81 = "cholelithiasis gall stones",
  M10 = "gout",
  M81 = "osteoporosis",
  N40 = "hyperplasia of prostate",
  C18 = "bowel cancer",
  C20 = "bowel cancer",
  C34 = "lung cancer",
  C43 = "malignant melanoma of skin",
  C44 = "skin cancer",
  C50 = "breast cancer",
  C54 = "endometrial cancer",
  C61 = "prostate cancer",
  C67 = "bladder cancer"
)

# Per-trait overrides applied last, keyed by the normalised full trait name
# (dedup key without category); NA removes a trait from any proxy group
proxy_trait_overrides <- c(
  "medication for cholesterol blood pressure or diabetes insulin" = "type 1 diabetes",
  "medication for cholesterol blood pressure diabetes or take exogenous hormones insulin" = "type 1 diabetes",
  "age diabetes diagnosed" = "type 1 diabetes",
  # "None of the above" answers mirror the dominant answer of their field
  "vascular heart problems diagnosed by doctor none of the above" = "hypertension",
  "medication for cholesterol blood pressure or diabetes none of the above" = "hypertension",
  "medication for cholesterol blood pressure diabetes or take exogenous hormones none of the above" = "hypertension",
  "blood clot dvt bronchitis emphysema asthma rhinitis eczema allergy diagnosed by doctor none of the above" =
    "allergic disease",
  # Composite outcomes dominated by one condition
  "cardiovascular disease" = "hypertension",
  "respiratory or ear nose throat disease" = "asthma"
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
  transmute(dedup_key, duplicate_group = row_number())

# ---- Proxy groups -----------------------------------------------------------
# A proxy concept is assigned per duplicate group from its representative's
# name, so duplicates always share a proxy group
representatives <- filter(ranked, keep)
proxy_name <- rules$lowercase(representatives$trait_name)

is_icd10 <- grepl(proxy_icd10_pattern, proxy_name, perl = TRUE)
icd10_code <- rep(NA_character_, length(proxy_name))
icd10_code[is_icd10] <- toupper(sub(proxy_icd10_pattern, "\\1", proxy_name[is_icd10], perl = TRUE))
proxy_name[is_icd10] <- sub(proxy_icd10_pattern, "\\2", proxy_name[is_icd10], perl = TRUE)
# Generic descriptions ("E11.9 Without complications") and descriptions used by
# more than one code ("F10.1 / F17.1 Harmful use") fall back to the code so
# they do not merge across different codes
shared_icd10_text <- tibble::tibble(text = proxy_name[is_icd10], code = icd10_code[is_icd10]) |>
  distinct() |>
  count(text) |>
  filter(n > 1L)
generic_icd10 <- is_icd10 &
  (grepl("^(without|with|other|unspecified)\\b", proxy_name) | proxy_name %in% shared_icd10_text$text)
proxy_name[generic_icd10] <- paste("icd10", tolower(icd10_code[generic_icd10]))

is_wrapped <- rep(FALSE, length(proxy_name))
for (pattern in proxy_wrappers) {
  hit <- !is_icd10 & !is_wrapped & grepl(pattern, proxy_name, perl = TRUE)
  proxy_name[hit] <- sub(pattern, "\\1", proxy_name[hit], perl = TRUE)
  is_wrapped <- is_wrapped | hit
}
proxy_name <- sub(",? unspecified$", "", proxy_name)

proxy_key <- normalise_trait_names(proxy_name, rules)$key
proxy_concept <- proxy_key
in_concept_map <- proxy_key %in% names(proxy_concept_map)
proxy_concept[in_concept_map] <- unname(proxy_concept_map[proxy_key[in_concept_map]])
in_icd10_map <- !is.na(icd10_code) & icd10_code %in% names(proxy_icd10_map)
proxy_concept[in_icd10_map] <- unname(proxy_icd10_map[icd10_code[in_icd10_map]])
proxy_concept[proxy_key %in% proxy_non_answers] <- NA_character_
full_name_key <- sub(" \\| [^|]*$", "", representatives$dedup_key)
in_overrides <- full_name_key %in% names(proxy_trait_overrides)
proxy_concept[in_overrides] <- unname(proxy_trait_overrides[full_name_key[in_overrides]])

proxy_rules <- vapply(seq_along(proxy_name), function(i) {
  applied <- c(
    wrapper = is_wrapped[i], icd10 = is_icd10[i],
    concept_map = in_concept_map[i], icd10_map = in_icd10_map[i], override = in_overrides[i]
  )
  return(paste(names(applied)[applied], collapse = ";"))
}, character(1))

# Only concepts spanning more than one duplicate group are proxy groups
proxy_assignments <- tibble::tibble(
  dedup_key = representatives$dedup_key,
  proxy_group = proxy_concept,
  proxy_rules = proxy_rules
) |>
  add_count(proxy_group, name = "proxy_group_size") |>
  mutate(
    proxy_group = if_else(proxy_group_size > 1L, proxy_group, NA_character_),
    proxy_rules = if_else(is.na(proxy_group), "", proxy_rules),
    # Display name: the concept in sentence case, with possessives and ICD10
    # codes restored ("alzheimer s disease" -> "Alzheimer's disease")
    proxy_group_name = gsub("\\b([a-z]+) s\\b", "\\1's", proxy_group),
    proxy_group_name = sub("^icd10 ([a-z][0-9]{2})$", "ICD10 \\U\\1", proxy_group_name, perl = TRUE),
    proxy_group_name = paste0(toupper(substr(proxy_group_name, 1, 1)), substring(proxy_group_name, 2)),
    proxy_group_name = if_else(is.na(proxy_group), NA_character_, proxy_group_name)
  ) |>
  select(dedup_key, proxy_group, proxy_group_name, proxy_rules)

dedup <- ranked |>
  inner_join(group_order, by = "dedup_key") |>
  inner_join(proxy_assignments, by = "dedup_key") |>
  rename(is_representative = keep) |>
  arrange(duplicate_group, group_rank) |>
  select(
    trait, trait_name, trait_category, category, sample_size,
    num_study_extractions, num_coloc_groups, duplicate_group, dedup_key,
    group_size, group_rank, is_representative, representative_trait,
    representative_trait_name, match_rules, proxy_group, proxy_group_name, proxy_rules
  )

# ---- Checks -----------------------------------------------------------------
group_checks <- dedup |>
  group_by(duplicate_group) |>
  summarise(
    n_keep = sum(is_representative),
    rep_in_group = all(representative_trait %in% trait),
    n_category = n_distinct(category),
    n_proxy_group = n_distinct(proxy_group),
    .groups = "drop"
  )
stopifnot(
  !anyDuplicated(dedup$trait),
  nrow(dedup) == nrow(traits),
  all(group_checks$n_keep == 1L),
  all(group_checks$rep_in_group),
  all(group_checks$n_category == 1L),
  all(group_checks$n_proxy_group == 1L)
)

readr::write_tsv(dedup, output_path, na = "")

# ---- Summary ----------------------------------------------------------------
multi <- filter(dedup, group_size > 1L)
message(
  "Traits in scope: ", nrow(dedup),
  " | duplicate groups: ", n_distinct(multi$duplicate_group),
  " | traits in duplicate groups: ", nrow(multi),
  " | traits dropped: ", sum(!dedup$is_representative),
  " | traits kept: ", sum(dedup$is_representative)
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
  distinct(duplicate_group, group_size, representative_trait, representative_trait_name) |>
  arrange(desc(group_size), duplicate_group) |>
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

# ---- Proxy summary ----------------------------------------------------------
proxied <- filter(dedup, !is.na(proxy_group))
message(
  "\nProxy groups: ", n_distinct(proxied$proxy_group),
  " | duplicate groups in proxy groups: ", n_distinct(proxied$duplicate_group),
  " | traits in proxy groups: ", nrow(proxied)
)

unused_concepts <- setdiff(names(proxy_concept_map), proxy_key)
if (length(unused_concepts) > 0) {
  message("\nProxy concept entries matching no trait: ", paste(unused_concepts, collapse = ", "))
}
unused_icd10 <- setdiff(names(proxy_icd10_map), icd10_code)
if (length(unused_icd10) > 0) {
  message("\nProxy ICD10 entries matching no trait: ", paste(unused_icd10, collapse = ", "))
}
unused_overrides <- setdiff(names(proxy_trait_overrides), full_name_key)
if (length(unused_overrides) > 0) {
  message("\nProxy trait overrides matching no trait: ", paste(unused_overrides, collapse = ", "))
}

message("\nLargest proxy groups:")
proxied |>
  filter(is_representative) |>
  group_by(proxy_group) |>
  arrange(desc(num_coloc_groups), .by_group = TRUE) |>
  summarise(
    n_duplicate_groups = n(),
    members = paste(head(trait_name, 6), collapse = " | "),
    .groups = "drop"
  ) |>
  arrange(desc(n_duplicate_groups)) |>
  head(n_review_proxy_groups) |>
  as.data.frame() |>
  print(row.names = FALSE, right = FALSE)

# Shared coloc groups between representatives, used only to flag proxies for
# review: high overlap alone cannot separate a proxy from a real mechanism
trait_coloc <- study_extractions |>
  inner_join(select(studies, study_name, trait), by = c("study" = "study_name")) |>
  inner_join(coloc_groups, by = "unique_study_id") |>
  filter(trait %in% representatives$trait, !is.na(coloc_group_id)) |>
  distinct(trait, coloc_group_id)
coloc_group_index <- match(trait_coloc$coloc_group_id, unique(trait_coloc$coloc_group_id))
membership <- Matrix::sparseMatrix(
  i = match(trait_coloc$trait, representatives$trait),
  j = coloc_group_index,
  x = 1,
  dims = c(nrow(representatives), max(coloc_group_index))
)
shared <- methods::as(Matrix::tcrossprod(membership), "generalMatrix")
n_trait_coloc_groups <- Matrix::diag(shared)
rep_proxy_group <- proxy_assignments$proxy_group
shared_pairs <- methods::as(shared, "TsparseMatrix")

message("\nMissed proxy candidates (high coloc overlap, a wrapper/ICD10 trait, not in the same proxy group):")
data.frame(i = shared_pairs@i + 1L, j = shared_pairs@j + 1L, x = shared_pairs@x) |>
  filter(i < j, x >= proxy_review_min_shared) |>
  mutate(
    overlap = x / pmin(n_trait_coloc_groups[i], n_trait_coloc_groups[j]),
    is_wrapper_pair = is_wrapped[i] | is_icd10[i] | is_wrapped[j] | is_icd10[j],
    same_proxy_group = coalesce(rep_proxy_group[i] == rep_proxy_group[j], FALSE)
  ) |>
  filter(overlap >= proxy_review_min_overlap, is_wrapper_pair, !same_proxy_group) |>
  arrange(desc(x)) |>
  head(n_review_proxy_pairs) |>
  transmute(
    trait_name_a = representatives$trait_name[i],
    proxy_group_a = rep_proxy_group[i],
    trait_name_b = representatives$trait_name[j],
    proxy_group_b = rep_proxy_group[j],
    shared_coloc_groups = x,
    overlap = round(overlap, 2)
  ) |>
  print(row.names = FALSE, right = FALSE)

message("\nSuspicious proxy members (low coloc overlap with the top trait of their proxy group):")
tibble::tibble(
  index = seq_len(nrow(representatives)),
  trait_name = representatives$trait_name,
  proxy_group = rep_proxy_group,
  num_coloc_groups = n_trait_coloc_groups
) |>
  filter(!is.na(proxy_group)) |>
  group_by(proxy_group) |>
  arrange(desc(num_coloc_groups), .by_group = TRUE) |>
  mutate(top_index = first(index), top_trait_name = first(trait_name)) |>
  ungroup() |>
  filter(index != top_index, num_coloc_groups >= proxy_review_min_shared) |>
  mutate(
    overlap = shared[cbind(top_index, index)] / pmin(n_trait_coloc_groups[top_index], num_coloc_groups)
  ) |>
  filter(overlap < proxy_review_max_member_overlap) |>
  arrange(proxy_group, overlap) |>
  transmute(proxy_group, trait_name, top_trait_name, num_coloc_groups, overlap = round(overlap, 2)) |>
  as.data.frame() |>
  print(row.names = FALSE, right = FALSE)

message("\nWrote ", output_path)
