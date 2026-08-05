#!/usr/bin/env Rscript

# Rank LD blocks by truthset T-I pairs and analyse pairwise-coloc gaps vs clustering.
#
# Uses compiled results from results/1.0.0/ (study_extractions, coloc_pairwise_results,
# coloc_clustered_results, studies_processed) for fast I/O.
#
# Usage:
#   Rscript rank_truthset_ld_blocks.R
#   Rscript rank_truthset_ld_blocks.R --top_n 20
#   Rscript rank_truthset_ld_blocks.R --skip_gap_detail
#   Rscript rank_truthset_ld_blocks.R --block_list ../../../pipeline_steps/data/1_0_0_block_list.csv
#   Rscript rank_truthset_ld_blocks.R --results_version 1.0.0

cmd_args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", cmd_args[grep("--file=", cmd_args)])
if (length(script_path) == 1L && nzchar(script_path)) {
  setwd(dirname(normalizePath(script_path)))
}

if (file.exists("../../../.env")) {
  readRenviron("../../../.env")
}
source("../../../pipeline_steps/constants.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(vroom)
})

parser <- argparser::arg_parser(
  "Rank LD blocks by truthset T-I pairs and analyse pairwise-coloc gaps"
)
parser <- argparser::add_argument(
  parser,
  "--top_n",
  help = "Number of top LD blocks to report",
  type = "integer",
  default = 10L
)
parser <- argparser::add_argument(
  parser,
  "--output_dir",
  help = "Directory for ranking / gap summary outputs",
  type = "character",
  default = "../data"
)
parser <- argparser::add_argument(
  parser,
  "--results_version",
  help = "Version label under RESULTS_DIR, or absolute path to a results directory",
  type = "character",
  default = "1.0.0"
)
parser <- argparser::add_argument(
  parser,
  "--results_dir",
  help = "Deprecated alias for --results_version",
  type = "character",
  default = NA
)
parser <- argparser::add_argument(
  parser,
  "--block_list",
  help = "Optional CSV of studies to exclude",
  type = "character",
  default = NA
)
parser <- argparser::add_argument(
  parser,
  "--skip_gap_detail",
  help = "Skip pairwise gap detail extraction",
  type = "logical",
  flag = TRUE
)
args <- argparser::parse_args(parser)

top_n <- as.integer(args$top_n)
output_dir <- args$output_dir
skip_gap_detail <- isTRUE(args$skip_gap_detail)
block_list_path <- args$block_list
results_subdir <- if (!is.na(args$results_dir) && nzchar(args$results_dir)) {
  sub("/$", "", args$results_dir)
} else {
  sub("/$", "", args$results_version)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (grepl("^/", results_subdir)) {
  published_results_dir <- paste0(results_subdir, "/")
} else {
  published_results_dir <- glue::glue("{results_dir}{results_subdir}/")
}

truthset_file <- file.path(published_results_dir, "analysis", "clustering", "tipairs_launched_truthset.tsv")
if (!file.exists(truthset_file)) {
  fallback_truthset <- file.path("..", "data", "tipairs_launched_truthset.tsv")
  if (file.exists(fallback_truthset)) {
    message("Using fallback truthset: ", fallback_truthset)
    truthset_file <- fallback_truthset
  }
}
study_extractions_file <- file.path(published_results_dir, "study_extractions.tsv.gz")
studies_file <- file.path(published_results_dir, "studies_processed.tsv.gz")
coloc_global_file <- file.path(published_results_dir, "coloc_pairwise_results.tsv.gz")
clusters_global_file <- file.path(published_results_dir, "coloc_clustered_results.tsv.gz")

block_list <- NULL
block_list_name <- NULL
if (!is.na(block_list_path) && nzchar(block_list_path) && file.exists(block_list_path)) {
  block_list <- vroom::vroom(block_list_path, show_col_types = FALSE)
  block_list_name <- get_block_list_name(block_list_path)
  message("Using block list: ", block_list_path, " (", block_list_name, ")")
}

pair_id_undirected <- function(a, b) {
  return(ifelse(a < b, paste(a, b, sep = "||"), paste(b, a, sep = "||")))
}

read_coloc_global_filtered <- function(coloc_path, blocks_keep, uid_keep) {
  blocks_keep <- unique(blocks_keep[nzchar(blocks_keep)])
  uid_keep <- unique(uid_keep[!is.na(uid_keep) & nzchar(uid_keep)])
  if (length(blocks_keep) == 0 || length(uid_keep) == 0 || !file.exists(coloc_path)) {
    return(data.table::data.table())
  }

  blocks_file <- tempfile("ld_blocks_", tmpdir = tempdir())
  uids_file <- tempfile("uids_", tmpdir = tempdir())
  writeLines(blocks_keep, blocks_file)
  writeLines(uid_keep, uids_file)

  # Columns: unique_study_a, unique_study_b, PP.H0-H4, ld_block (8), ...
  awk_prog <- sprintf(
    paste(
      "BEGIN {",
      "  while ((getline < \"%s\") > 0) block[$1]=1; close(\"%s\");",
      "  while ((getline < \"%s\") > 0) uid[$1]=1; close(\"%s\")",
      "}",
      "FNR==1 { print; next }",
      "($8 in block && $1 in uid && $2 in uid)"
    ),
    blocks_file, blocks_file, uids_file, uids_file
  )

  cmd <- sprintf("zcat '%s' | awk -F'\\t' '%s'", coloc_path, awk_prog)
  df <- tryCatch(
    {
      out <- data.table::fread(cmd = cmd, showProgress = FALSE)
      unlink(c(blocks_file, uids_file))
      out
    },
    error = function(e) {
      unlink(c(blocks_file, uids_file))
      return(data.table::data.table())
    }
  )
  if (nrow(df) == 0) {
    return(df)
  }
  if (!"h4" %in% names(df) && "PP.H4.abf" %in% names(df)) {
    df[, h4 := PP.H4.abf]
  }
  if (!"ld_block" %in% names(df)) {
    names(df)[8] <- "ld_block"
  }
  return(df)
}

read_block_clusters_filtered <- function(ld_block, uid_keep) {
  uid_keep <- unique(uid_keep[!is.na(uid_keep)])
  if (length(uid_keep) == 0) {
    return(NULL)
  }

  block_dir <- file.path(ld_block_data_dir, ld_block)
  cluster_fnames <- if (!is.null(block_list_name)) {
    c(
      glue::glue("coloc_clustered_results_bfdr_{block_list_name}.tsv.gz"),
      glue::glue("coloc_clustered_results_{block_list_name}.tsv.gz")
    )
  } else {
    c("coloc_clustered_results_bfdr.tsv.gz", "coloc_clustered_results.tsv.gz")
  }

  for (fname in cluster_fnames) {
    path <- file.path(block_dir, fname)
    if (!file.exists(path)) {
      next
    }
    df <- data.table::fread(
      path,
      select = c("unique_study_id", "component", "h4_connectedness"),
      showProgress = FALSE
    )
    df <- df[unique_study_id %in% uid_keep]
    if (nrow(df) == 0) {
      return(NULL)
    }
    df[, ld_block := ld_block]
    return(df)
  }
  return(NULL)
}

read_block_bfdr_info <- function(ld_block) {
  block_dir <- file.path(ld_block_data_dir, ld_block)
  rds_fnames <- if (!is.null(block_list_name)) {
    c(
      glue::glue("igraph_clustered_results_{block_list_name}.rds"),
      "igraph_clustered_results.rds"
    )
  } else {
    "igraph_clustered_results.rds"
  }

  for (fname in rds_fnames) {
    path <- file.path(block_dir, fname)
    if (!file.exists(path)) {
      next
    }
    obj <- tryCatch(readRDS(path), error = function(e) NULL)
    if (is.null(obj) || is.null(obj$bfdr_info)) {
      return(NULL)
    }
    return(obj$bfdr_info)
  }
  return(NULL)
}

message("Published results dir: ", published_results_dir)
message("Loading truthset and study metadata...")

for (required_file in c(truthset_file, study_extractions_file, studies_file)) {
  if (!file.exists(required_file)) {
    stop("Required file not found: ", required_file)
  }
}

truthset <- vroom::vroom(truthset_file, show_col_types = FALSE)
studies_processed <- vroom::vroom(studies_file, show_col_types = FALSE)

gene_study_lookup <- studies_processed |>
  dplyr::filter(data_type != "phenotype", !is.na(gene)) |>
  dplyr::select(study_name, gene, data_type, tissue)

study_metadata <- studies_processed |>
  dplyr::select(
    study_name, data_type, gene, tissue, trait_name, variant_type, category
  )

message("Loading study_extractions from ", study_extractions_file)
study_extractions <- data.table::fread(study_extractions_file, showProgress = FALSE)

finemap_index <- tibble::as_tibble(study_extractions) |>
  dplyr::left_join(study_metadata, by = c("study" = "study_name")) |>
  dplyr::filter(variant_type == "common")

if (!is.null(block_list) && nrow(block_list) > 0) {
  blocked <- is_study_blocked(
    block_list,
    finemap_index$study,
    finemap_index$cis_trans
  )
  n_blocked <- sum(blocked)
  if (n_blocked > 0) {
    message("Excluding ", n_blocked, " finemapped rows matching block list")
    finemap_index <- finemap_index |> dplyr::filter(!blocked)
  }
}

message("Finemap index rows (common, from study_extractions): ", nrow(finemap_index))

study_blocks <- finemap_index |>
  dplyr::distinct(study, ld_block)

# --- 1. Support tier summary ---
truthset_pairwise <- truthset |>
  dplyr::filter(pairwisecoloc_support == TRUE) |>
  dplyr::mutate(
    support_tier = dplyr::case_when(
      gpmap_support & preinfomap_support ~ "pairwise_preinfomap_and_gpmap",
      !gpmap_support & preinfomap_support ~ "pairwise_and_preinfomap_only",
      gpmap_support & !preinfomap_support ~ "pairwise_and_gpmap_only",
      TRUE ~ "pairwise_only"
    ),
    missing_gpmap_cluster = !gpmap_support,
    missing_preinfomap_cluster = !preinfomap_support
  )

support_summary <- truthset_pairwise |>
  dplyr::count(support_tier, name = "n_rows") |>
  dplyr::arrange(dplyr::desc(n_rows))

support_summary_ti <- truthset_pairwise |>
  dplyr::count(support_tier, name = "n_rows") |>
  dplyr::mutate(by = "truthset_row") |>
  dplyr::bind_rows(
    truthset_pairwise |>
      dplyr::distinct(ti_uid, support_tier) |>
      dplyr::count(support_tier, name = "n_rows") |>
      dplyr::mutate(by = "ti_uid")
  )

message("\nTruthset rows with pairwise coloc support:")
print(support_summary)

truthset_gap <- truthset_pairwise |>
  dplyr::filter(missing_gpmap_cluster | missing_preinfomap_cluster)

message("\nGap rows (pairwise but missing gpmap and/or preinfomap): ", nrow(truthset_gap))
message("Unique ti_uids in gap set: ", dplyr::n_distinct(truthset_gap$ti_uid))

# --- 2. LD block rankings ---
truthset_gap_genes <- truthset_gap |>
  dplyr::inner_join(
    gene_study_lookup,
    by = c("target" = "gene"),
    relationship = "many-to-many"
  )

mol_in_block <- truthset_gap_genes |>
  dplyr::inner_join(
    study_blocks |> dplyr::rename(study_name = study),
    by = "study_name",
    relationship = "many-to-many"
  )

truthset_in_block <- mol_in_block |>
  dplyr::inner_join(
    study_blocks |> dplyr::rename(trait_gpmap = study),
    by = c("trait_gpmap", "ld_block")
  )

rank_by_ti_uid <- truthset_in_block |>
  dplyr::distinct(ld_block, ti_uid) |>
  dplyr::count(ld_block, name = "n_truthset_ti_uids") |>
  dplyr::arrange(dplyr::desc(n_truthset_ti_uids))

rank_by_row <- truthset_in_block |>
  dplyr::distinct(ld_block, ti_uid, trait_gpmap) |>
  dplyr::count(ld_block, name = "n_truthset_rows") |>
  dplyr::arrange(dplyr::desc(n_truthset_rows))

rankings <- rank_by_ti_uid |>
  dplyr::left_join(rank_by_row, by = "ld_block") |>
  dplyr::arrange(dplyr::desc(n_truthset_ti_uids), dplyr::desc(n_truthset_rows))

message("\nTop ", top_n, " LD blocks by truthset ti_uids:")
print(dplyr::slice_head(rankings, n = top_n))

rankings_file <- file.path(output_dir, "truthset_ld_block_rankings.tsv")
support_file <- file.path(output_dir, "truthset_support_summary.tsv")
vroom::vroom_write(rankings, rankings_file)
vroom::vroom_write(support_summary_ti, support_file)

if (skip_gap_detail) {
  message("\nSkipping gap detail (--skip_gap_detail). Wrote rankings only.")
  message("  ", rankings_file)
  message("  ", support_file)
  message("\nDone.")
  quit(save = "no", status = 0)
}

# --- 3. Gap detail ---
message("\nBuilding pairwise gap detail with finemap metadata...")

gap_blocks <- unique(truthset_in_block$ld_block)
message("Gap rows map to ", length(gap_blocks), " LD blocks")

gap_detail <- truthset_in_block |>
  dplyr::distinct(
    ld_block, ti_uid, target, trait_gpmap, trait_name_gpmap,
    indication_mesh_term, combined_max_phase,
    gpmap_support, preinfomap_support, pairwisecoloc_support,
    study_name, data_type, tissue
  ) |>
  dplyr::rename(mol_study_name = study_name, mol_data_type = data_type, mol_tissue = tissue)

needed_mol <- gap_detail |> dplyr::distinct(ld_block, mol_study_name)
needed_pheno <- gap_detail |> dplyr::distinct(ld_block, trait_gpmap)

mol_finemap <- finemap_index |>
  dplyr::inner_join(needed_mol, by = c("ld_block", "study" = "mol_study_name")) |>
  dplyr::rename_with(~ paste0("mol_", .), -c(study, ld_block)) |>
  dplyr::rename(mol_study_name = study)

pheno_finemap <- finemap_index |>
  dplyr::inner_join(needed_pheno, by = c("ld_block", "study" = "trait_gpmap")) |>
  dplyr::rename_with(~ paste0("pheno_", .), -c(study, ld_block)) |>
  dplyr::rename(trait_gpmap = study)

# Multiple credible sets per study per block -> many-to-many is expected
gap_detail <- gap_detail |>
  dplyr::left_join(
    mol_finemap,
    by = c("ld_block", "mol_study_name"),
    relationship = "many-to-many"
  ) |>
  dplyr::left_join(
    pheno_finemap,
    by = c("ld_block", "trait_gpmap"),
    relationship = "many-to-many"
  )

bfdr_info_by_block <- stats::setNames(
  lapply(gap_blocks, read_block_bfdr_info),
  gap_blocks
)

gap_detail <- gap_detail |>
  dplyr::mutate(
    block_bfdr_threshold = vapply(ld_block, function(b) {
      info <- bfdr_info_by_block[[b]]
      if (is.null(info)) {
        return(NA_real_)
      }
      return(info$threshold)
    }, numeric(1)),
    h4_threshold_used = dplyr::coalesce(block_bfdr_threshold, posterior_prob_h4_threshold)
  )

uid_keep_by_block <- dplyr::bind_rows(
  gap_detail |> dplyr::select(ld_block, uid = mol_unique_study_id),
  gap_detail |> dplyr::select(ld_block, uid = pheno_unique_study_id)
) |>
  dplyr::filter(!is.na(uid)) |>
  dplyr::distinct()

all_uids <- unique(uid_keep_by_block$uid)

message("Loading coloc pairs from global file (single pass): ", coloc_global_file)
if (!file.exists(coloc_global_file)) {
  warning("Global coloc file not found: ", coloc_global_file)
  coloc_all <- data.table::data.table()
} else {
  coloc_all <- read_coloc_global_filtered(coloc_global_file, gap_blocks, all_uids)
  message("Coloc pairs loaded: ", nrow(coloc_all))
}

if (!is.null(block_list_name)) {
  message("Loading per-block cluster assignments (block_list: ", block_list_name, ")")
  cluster_all <- lapply(gap_blocks, function(b) {
    uids <- uid_keep_by_block$uid[uid_keep_by_block$ld_block == b]
    if (length(uids) == 0) {
      return(NULL)
    }
    return(read_block_clusters_filtered(b, unique(uids)))
  })
  cluster_all <- data.table::rbindlist(cluster_all, use.names = TRUE, fill = TRUE)
} else {
  message("Loading cluster assignments from global file: ", clusters_global_file)
  if (!file.exists(clusters_global_file)) {
    cluster_all <- data.table::data.table()
  } else {
    cluster_all <- data.table::fread(clusters_global_file, showProgress = FALSE)
    cluster_cols <- intersect(
      names(cluster_all),
      c("unique_study_id", "component", "ld_block", "h4_connectedness")
    )
    cluster_all <- cluster_all[, ..cluster_cols]
    cluster_all <- cluster_all[
      ld_block %in% gap_blocks & unique_study_id %in% all_uids
    ]
  }
  message("Cluster rows loaded: ", nrow(cluster_all))
}

if (nrow(coloc_all) > 0) {
  coloc_pairs <- tibble::as_tibble(coloc_all) |>
    dplyr::mutate(
      pair_id = pair_id_undirected(unique_study_a, unique_study_b),
      h4_val = dplyr::coalesce(.data$h4, .data$PP.H4.abf)
    )

  for (col in c("ignore", "false_positive", "false_negative", "bp_distance")) {
    if (!col %in% names(coloc_pairs)) {
      if (col == "bp_distance") {
        coloc_pairs[[col]] <- NA_real_
      } else {
        coloc_pairs[[col]] <- NA
      }
    }
  }

  coloc_pairs <- coloc_pairs |>
    dplyr::select(
      ld_block, pair_id, h4_val, ignore, false_positive, false_negative, bp_distance
    )

  gap_detail <- gap_detail |>
    dplyr::mutate(
      pair_id = pair_id_undirected(mol_unique_study_id, pheno_unique_study_id)
    ) |>
    dplyr::left_join(coloc_pairs, by = c("ld_block", "pair_id"))
} else {
  gap_detail$h4_val <- NA_real_
  gap_detail$ignore <- NA
  gap_detail$false_positive <- NA
  gap_detail$false_negative <- NA
  gap_detail$bp_distance <- NA_real_
}

if (nrow(cluster_all) > 0) {
  cluster_lookup <- tibble::as_tibble(cluster_all) |>
    dplyr::select(dplyr::any_of(c("ld_block", "unique_study_id", "component", "h4_connectedness")))

  gap_detail <- gap_detail |>
    dplyr::left_join(
      cluster_lookup |> dplyr::rename(mol_unique_study_id = unique_study_id, mol_component = component),
      by = c("ld_block", "mol_unique_study_id"),
      relationship = "many-to-many"
    ) |>
    dplyr::left_join(
      cluster_lookup |> dplyr::rename(pheno_unique_study_id = unique_study_id, pheno_component = component),
      by = c("ld_block", "pheno_unique_study_id"),
      relationship = "many-to-many"
    ) |>
    dplyr::mutate(
      same_cluster_component = !is.na(mol_component) &
        !is.na(pheno_component) &
        mol_component == pheno_component,
      missing_from_clusters = is.na(mol_component) | is.na(pheno_component)
    )
} else {
  gap_detail$same_cluster_component <- NA
  gap_detail$missing_from_clusters <- NA
}

gap_detail <- gap_detail |>
  dplyr::mutate(
    bp_distance_finemap = abs(mol_bp - pheno_bp),
    weaker_min_p = pmax(mol_min_p, pheno_min_p),
    either_ignored = dplyr::coalesce(mol_ignore, FALSE) | dplyr::coalesce(pheno_ignore, FALSE),
    gap_reason_hint = dplyr::case_when(
      is.na(h4_val) ~ "no_coloc_pair_in_block",
      dplyr::coalesce(ignore, FALSE) ~ "coloc_pair_ignored",
      h4_val < h4_threshold_used ~ "coloc_h4_below_block_threshold",
      !dplyr::coalesce(same_cluster_component, FALSE) & !dplyr::coalesce(missing_from_clusters, TRUE) ~
        "coloc_ok_but_different_clusters",
      dplyr::coalesce(missing_from_clusters, FALSE) ~ "study_missing_from_cluster_output",
      either_ignored ~ "study_flagged_ignore_in_finemap",
      weaker_min_p > lowest_p_value_threshold ~ "weak_finemap_min_p",
      TRUE ~ "other"
    )
  ) |>
  dplyr::relocate(
    ti_uid, target, trait_gpmap, ld_block, gap_reason_hint,
    gpmap_support, preinfomap_support, pairwisecoloc_support,
    mol_study_name, mol_unique_study_id, mol_min_p, mol_bp, mol_cis_trans, mol_ignore,
    pheno_unique_study_id, pheno_min_p, pheno_bp, pheno_ignore,
    h4_val, ignore, same_cluster_component, missing_from_clusters,
    h4_threshold_used, block_bfdr_threshold,
    .before = dplyr::everything()
  )

gap_row_summary <- gap_detail |>
  dplyr::group_by(ti_uid, trait_gpmap, target) |>
  dplyr::summarise(
    n_blocks = dplyr::n_distinct(ld_block),
    n_mol_studies_in_blocks = dplyr::n_distinct(mol_study_name),
    best_h4 = if (all(is.na(h4_val))) NA_real_ else max(h4_val, na.rm = TRUE),
    any_same_cluster = any(same_cluster_component %in% TRUE, na.rm = TRUE),
    any_coloc_pair = any(!is.na(h4_val)),
    median_mol_min_p = median(mol_min_p, na.rm = TRUE),
    median_pheno_min_p = median(pheno_min_p, na.rm = TRUE),
    top_gap_reason = names(sort(table(gap_reason_hint), decreasing = TRUE))[1],
    .groups = "drop"
  ) |>
  dplyr::left_join(
    truthset_gap |>
      dplyr::select(
        ti_uid, trait_gpmap, gpmap_support, preinfomap_support,
        indication_mesh_term, combined_max_phase
      ),
    by = c("ti_uid", "trait_gpmap")
  )

gap_detail_file <- file.path(output_dir, "truthset_pairwise_gap_detail.tsv")
gap_summary_file <- file.path(output_dir, "truthset_pairwise_gap_summary.tsv")
vroom::vroom_write(gap_detail, gap_detail_file)
vroom::vroom_write(gap_row_summary, gap_summary_file)

message("\nWrote:")
message("  ", rankings_file)
message("  ", support_file)
message("  ", gap_detail_file)
message("  ", gap_summary_file)

message("\nGap reason counts (detail rows):")
print(gap_detail |> dplyr::count(gap_reason_hint, sort = TRUE))

message("\nDone.")
