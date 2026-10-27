#!/usr/bin/env Rscript
# Build pre-infomap coloc clusters (connected components of each LD block's H4 graph, before
# infomap splitting and pruning) as coloc_clustered_results_no_pruning.tsv.gz, read by
# 03_ti_pair_validation_of_coloc.Rmd.
#
# Usage:
#   Rscript 00_build_preinfomap_clusters.R --results_version current
#   Rscript 00_build_preinfomap_clusters.R --results_version 1.0.0 --out /tmp/no_pruning.tsv.gz

cmd_args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", cmd_args[grep("--file=", cmd_args)])
if (length(script_path) == 1L && nzchar(script_path)) {
  setwd(dirname(normalizePath(script_path)))
}

source("../../../pipeline_steps/constants.R")

parser <- argparser::arg_parser("Build pre-infomap coloc clusters from per-block unpruned H4 graphs")
parser <- argparser::add_argument(
  parser,
  "--results_version",
  help = "Version label under RESULTS_DIR, or absolute path to a results directory",
  type = "character",
  default = "current"
)
parser <- argparser::add_argument(
  parser,
  "--rds_suffix",
  help = "Suffix of the per-block igraph_clustered_results{suffix}.rds files to read",
  type = "character",
  default = ""
)
parser <- argparser::add_argument(
  parser,
  "--out",
  help = "Output file (default: <results>/coloc_clustered_results_no_pruning.tsv.gz)",
  type = "character",
  default = NA
)
args <- argparser::parse_args(parser)

versioned_results_dir <- if (grepl("^/", args$results_version)) {
  args$results_version
} else {
  file.path(sub("/$", "", results_dir), args$results_version)
}
out_file <- if (is.na(args$out)) {
  file.path(versioned_results_dir, "coloc_clustered_results_no_pruning.tsv.gz")
} else {
  args$out
}
if (file.exists(out_file)) {
  stop("Refusing to overwrite existing ", out_file)
}

rds_name <- glue::glue("igraph_clustered_results{args$rds_suffix}.rds")
rds_files <- Sys.glob(glue::glue("{ld_block_data_dir}*/*/*/{rds_name}"))
message("Reading ", length(rds_files), " ", rds_name, " files")

pre_infomap_clusters <- do.call("rbind", lapply(rds_files, function(rds_file) {
  graph <- readRDS(rds_file)$unpruned_graph
  if (is.null(graph) || igraph::vcount(graph) == 0) {
    return(NULL)
  }
  graph_components <- igraph::components(graph)
  # Singletons are not clusters
  keep <- graph_components$csize[graph_components$membership] > 1
  if (!any(keep)) {
    return(NULL)
  }
  return(data.frame(
    unique_study_id = names(graph_components$membership)[keep],
    component = unname(graph_components$membership[keep]),
    ld_block = sub(paste0("^", ld_block_data_dir), "", dirname(rds_file))
  ))
}))

study_metadata <- data.table::fread(file.path(versioned_results_dir, "studies_processed.tsv.gz"))

pre_infomap_clusters <- pre_infomap_clusters |>
  dplyr::mutate(study = sub("_EUR.*", "", unique_study_id), .before = unique_study_id) |>
  dplyr::group_by(ld_block, component) |>
  dplyr::mutate(coloc_group_id = dplyr::cur_group_id()) |>
  dplyr::ungroup() |>
  dplyr::arrange(coloc_group_id) |>
  dplyr::left_join(study_metadata[, c("study_name", "data_type", "gene")], by = c("study" = "study_name"))

message(
  "Pre-infomap clusters: ", nrow(pre_infomap_clusters), " rows, ",
  length(unique(pre_infomap_clusters$coloc_group_id)), " groups, ",
  length(unique(pre_infomap_clusters$ld_block)), " LD blocks"
)
vroom::vroom_write(pre_infomap_clusters, out_file)
message("Wrote: ", out_file)
