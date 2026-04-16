options(error = function() traceback(20))
Sys.setenv("VROOM_CONNECTION_SIZE" = 500000)

data_dir <- Sys.getenv("DATA_DIR")
gwas_upload_dir <- Sys.getenv("DATA_DIR")
results_dir <- Sys.getenv("RESULTS_DIR")
oracle_api_server <- Sys.getenv("ORACLE_SERVER")
oracle_upload_server <- Sys.getenv("ORACLE_UPLOAD_SERVER")
TEST_RUN <- Sys.getenv("TEST_RUN", NA)
is_test_run <- !is.na(TEST_RUN)

min_p_allowed_for_worker <- 1e-6
genome_wide_p_value_threshold <- 5e-8
lowest_p_value_threshold <- 1.5e-4

minimum_extraction_size_for_dense_coverage <- 150
minimum_extraction_size_for_sparse_coverage <- 4

posterior_prob_h4_threshold <- 0.8
posterior_prob_threshold_minimum <- 0.5


gpm_website_data <- list(
  url = "https://gpmap.opengwas.io",
  contact = "https://gpmap.opengwas.io/contact",
  name = "The Genotype-Phenotype Map Team"
)

latest_results_dir <- glue::glue("{results_dir}latest/")
current_results_dir <- glue::glue("{results_dir}current/")
results_analysis_dir <- glue::glue("{latest_results_dir}analysis/")
pipeline_metadata_dir <- glue::glue("{data_dir}pipeline_metadata/")
ld_block_data_dir <- glue::glue("{data_dir}ld_blocks/")
ld_reference_panel_dir <- glue::glue("{data_dir}ld_reference_panel_hg38/")
liftover_dir <- glue::glue("{data_dir}liftover/")
extracted_study_dir <- glue::glue("{data_dir}study/")
variant_annotation_dir <- glue::glue("{data_dir}variant_annotation/")
static_web_dir <- glue::glue("{current_results_dir}static_web/")
svg_dir <- glue::glue("{static_web_dir}svgs/")

oracle_bucket_name <- Sys.getenv("ORACLE_BUCKET_NAME")
server_sync_dir <- file.path(data_dir, "rsync_to_server")
oracle_data_dir <- "/oradiskvdb1/data/"

bespoke_parsing_options <- list(
  none = "none",
  gtex_sqtl = "gtex_sqtl",
  interval_sqtl = "interval_sqtl",
  godmc = "godmc"
)

data_types <- list(
  splice_variant = "splice_variant",
  transcript = "transcript",
  gene_expression = "gene_expression",
  protein = "protein",
  methylation = "methylation",
  metabolome = "metabolome",
  cell_trait = "cell_trait",
  plasma_protein = "plasma_protein",
  phenotype = "phenotype"
)

non_qtl_data_types <- c(
  data_types$metabolome,
  data_types$cell_trait,
  data_types$plasma_protein,
  data_types$phenotype
)
qtl_data_types <- c(
  data_types$splice_variant,
  data_types$transcript,
  data_types$gene_expression,
  data_types$protein,
  data_types$methylation
)

data_type_names <- list(
  splice_variant = "sQTL",
  transcript = "tQTL",
  gene_expression = "eQTL",
  protein = "pQTL",
  methylation = "methQTL",
  metabolome = "metaQTL",
  cell_trait = "Cell Trait",
  plasma_protein = "Targeted Protein Measure",
  phenotype = "Phenotype"
)
cell_types <- c(
  "B IN",
  "B Mem",
  "CD4 ET",
  "CD4 NC",
  "CD4 SOX4",
  "CD8 ET",
  "CD8 NC",
  "CD8 S100B",
  "DC",
  "Mono C",
  "Mono NC",
  "NK",
  "NKR",
  "Plasma"
)
study_categories <- list(continuous = "continuous", categorical = "categorical")
data_formats <- list(opengwas = "opengwas", besd = "besd", tsv = "tsv")
cis_trans <- list(cis_only = "cis", trans_only = "trans", cis_trans = "cis_trans")
variant_types <- list(common = "common", rare_exome = "rare_exome", rare_wgs = "rare_wgs")
ancestry_map <- list(EUR = "European", EAS = "East Asian", AFR = "African", SAS = "South Asian")
reverse_ancestry_map <- setNames(names(ancestry_map), ancestry_map)

reference_builds <- list(GRCh37 = "GRCh37", GRCh38 = "GRCh38")
available_liftover_conversions <- list(
  "GRCh36GRCh37" = glue::glue("{liftover_dir}hg18ToHg19.over.chain.gz"),
  "GRCh38GRCh37" = glue::glue("{liftover_dir}hg38ToHg19.over.chain.gz"),
  "GRCh37GRCh38" = glue::glue("{liftover_dir}/hg19ToHg38.over.chain.gz")
)
extraction_file_types <- list(vcf = "vcf", csv = "csv")
coverage_types <- list(dense = "dense", sparse = "sparse")

standardised_gwas_columns <- c("CHR", "BP", "EA", "OA", "EAF", "BETA", "SE", "P", "SNP", "Z", "GENE")
required_columns <- c("CHR", "BP", "EA", "OA", "EAF", "BETA", "SE", "P")
beta_columns <- c("BETA", "SE")
or_columns <- c("OR", "OR_LB", "OR_UB")

standardised_column_types <- vroom::cols(
  chr = vroom::col_character(),
  bp = vroom::col_number(),
  sample_size = vroom::col_number(),
  p_value_threshold = vroom::col_number(),
  snps_removed_by_reference_panel = vroom::col_number(),
  eaf_from_reference_panel = vroom::col_logical(),
  time_taken = vroom::col_character()
)

imputed_column_types <- vroom::cols(
  chr = vroom::col_character(),
  bp = vroom::col_number(),
  sample_size = vroom::col_number(),
  p_value_threshold = vroom::col_number(),
  time_taken = vroom::col_character()
)

finemapped_column_types <- vroom::cols(
  chr = vroom::col_character(),
  bp = vroom::col_number(),
  min_p = vroom::col_number(),
  sample_size = vroom::col_number(),
  p_value_threshold = vroom::col_number(),
  first_finemap_num_results = vroom::col_number(),
  second_finemap_num_results = vroom::col_number(),
  qc_step_run = vroom::col_logical(),
  snps_removed_by_qc = vroom::col_number(),
  time_taken = vroom::col_character(),
  cis_trans = vroom::col_character(),
  ignore = vroom::col_logical(),
  file = vroom::col_character(),
  svg_file = vroom::col_character(),
  file_with_lbfs = vroom::col_character(),
  coverage = vroom::col_character()
)

studies_processed_column_types <- vroom::cols(
  data_type = vroom::col_character(),
  data_format = vroom::col_character(),
  study_name = vroom::col_character(),
  trait = vroom::col_character(),
  ancestry = vroom::col_character(),
  sample_size = vroom::col_number(),
  category = vroom::col_character(),
  study_location = vroom::col_character(),
  extracted_location = vroom::col_character(),
  reference_build = vroom::col_character(),
  p_value_threshold = vroom::col_number(),
  gene = vroom::col_character(),
  probe = vroom::col_character(),
  tissue = vroom::col_character(),
  source = vroom::col_character(),
  variant_type = vroom::col_character(),
  trait_name = vroom::col_character(),
  ensg = vroom::col_character(),
  coverage = vroom::col_character(),
  heritability = vroom::col_double(),
  heritability_se = vroom::col_double()
)

coloc_clustered_results_column_types <- vroom::cols(
  unique_study_id = vroom::col_character(),
  component = vroom::col_number(),
  ld_block = vroom::col_character(),
  snp = vroom::col_character(),
  h4_connectedness = vroom::col_double(),
  h3_connectedness = vroom::col_double()
)

coloc_pairwise_results_column_types <- vroom::cols(
  unique_study_a = vroom::col_character(),
  study_a = vroom::col_character(),
  unique_study_b = vroom::col_character(),
  study_b = vroom::col_character(),
  bp_distance = vroom::col_double(),
  ignore = vroom::col_logical(),
  false_positive = vroom::col_logical(),
  false_negative = vroom::col_logical(),
  nsnps = vroom::col_number(),
  hit1 = vroom::col_character(),
  hit2 = vroom::col_character(),
  PP.H0.abf = vroom::col_double(),
  PP.H1.abf = vroom::col_double(),
  PP.H2.abf = vroom::col_double(),
  PP.H3.abf = vroom::col_double(),
  PP.H4.abf = vroom::col_double(),
  idx1 = vroom::col_number(),
  idx2 = vroom::col_number(),
  h4 = vroom::col_double(),
  ld_block = vroom::col_character()
)

file_prefix <- function(file_path) {
  file_name <- basename(file_path)
  file_prefix <- sub("\\..*", "", file_name)
  return(file_prefix)
}

ld_block_dirs <- function(block) {
  ld_info <- data.frame(
    block = block,
    ld_block_data = glue::glue("{ld_block_data_dir}{block}"),
    ld_reference_panel_prefix = glue::glue("{ld_reference_panel_dir}{block}")
  )
  # ld_info <- dplyr::bind_cols(ld_info, ld_block_components(block))
  return(ld_info)
}

construct_ld_block <- function(ancestry, chr, start, stop) {
  block <- ld_block_string(ancestry, chr, start, stop)
  ld_info <- ld_block_dirs(block)
  ld_info <- dplyr::bind_cols(ld_info, data.frame(ancestry = ancestry, chr = chr, start = start, stop = stop))
  return(ld_info)
}

ld_block_string <- function(ancestry, chr, start, stop) {
  return(glue::glue("{ancestry}/{chr}/{start}-{stop}"))
}

ld_block_components <- function(ld_block) {
  components <- strsplit(ld_block, "[/-]")[[1]]
  return(data.frame(
    ancestry = components[1],
    chr = components[2],
    start = as.numeric(components[3]),
    stop = as.numeric(components[4])
  ))
}

flattened_ld_block_name <- function(ld_block_string) {
  return(gsub("[/-]", "_", ld_block_string))
}

update_directories_for_worker <- function(worker_guid) {
  ld_block_data_dir <<- glue::glue("{gwas_upload_dir}ld_blocks/gwas_upload/{worker_guid}/")
  extracted_study_dir <<- glue::glue("{gwas_upload_dir}gwas_upload/{worker_guid}/")
  pipeline_metadata_dir <<- glue::glue("{gwas_upload_dir}pipeline_metadata/gwas_upload/{worker_guid}/")
  return()
}

diff_time_taken <- function(start_time) {
  return(hms::as_hms(difftime(Sys.time(), start_time)))
}

safe_lapply <- function(X, FUN, ...) {
  return(lapply(X, function(x) {
    return(tryCatch(
      FUN(x, ...),
      error = function(e) {
        stop("An error occurred: ", conditionMessage(e))
      }
    ))
  }))
}

replace_except_first_two_dashes <- function(x) {
  dash_positions <- gregexpr("-", x)[[1]]
  if (length(dash_positions) <= 2) {
    return(x)
  }
  result <- x
  for (i in 3:length(dash_positions)) {
    pos <- dash_positions[i]
    substr(result, pos, pos) <- "_"
  }
  return(result)
}

is_study_blocked <- function(block_list, study_name, cis_trans) {
  if (is.null(block_list) || nrow(block_list) == 0) {
    return(rep(FALSE, length(study_name)))
  }

  return(vapply(seq_along(study_name), function(i) {
    s <- study_name[i]
    c <- cis_trans[i]

    return(any(vapply(seq_len(nrow(block_list)), function(j) {
      pattern <- block_list[["id_pattern"]][j]
      pattern_regex <- gsub("\\.", "\\\\.", gsub("\\*", ".*", pattern))

      if (!grepl(pattern_regex, s)) {
        return(FALSE)
      }
      block_cis <- block_list[["cis_trans"]][j]
      if (is.na(block_cis) || identical(block_cis, "NA") || trimws(as.character(block_cis)) == "") {
        return(TRUE)
      }
      return(identical(as.character(c), as.character(block_cis)))
    }, logical(1))))
  }, logical(1)))
}

get_block_list_name <- function(block_list) {
  if (!is.null(block_list) && !is.na(block_list) && file.exists(block_list)) {
    return(sub("_block_list.csv", "", basename(block_list)))
  }
  return(NULL)
}