#!/usr/bin/env Rscript
# Generate a list of launched T-I pairs that have colocalisation evidence as a
# "truth set" to assess clustering methodologies.
#
# Usage:
#   Rscript 04_pull_truth_links.R
#   Rscript 04_pull_truth_links.R --results_version 1.0.0

cmd_args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", cmd_args[grep("--file=", cmd_args)])
if (length(script_path) == 1L && nzchar(script_path)) {
  setwd(dirname(normalizePath(script_path)))
}

source("../../../pipeline_steps/constants.R")

parser <- argparser::arg_parser(
  "Build launched T-I pair truth set from scored GPMAP and pre-infomap support"
)
parser <- argparser::add_argument(
  parser,
  "--results_version",
  help = "Version label under RESULTS_DIR, or absolute path to a results directory",
  type = "character",
  default = "1.0.0"
)
args <- argparser::parse_args(parser)

versioned_results_dir <- if (grepl("^/", args$results_version)) {
  args$results_version
} else {
  file.path(sub("/$", "", results_dir), args$results_version)
}
analysis_dir <- file.path(versioned_results_dir, "analysis", "clustering")
dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)

## Extended table of T-I pairs with GPMAP evidence (from 02_gpmap_support_for_ti_pairs.Rmd)
ti_pairs <- data.table::fread(file.path(analysis_dir, "GPMAP_T-Ipairs_allmatchedstudies.tsv"))

## Colocalisation group evidence (studies linked by H4>0.8) pre infomap refinement of clusters
preinfomap_file <- file.path(ti_pairs_data_dir, "tipairs_preinfomap.rda")
if (!file.exists(preinfomap_file)) {
  stop("Missing ", preinfomap_file)
}
load(preinfomap_file) # ti_preinfomap
message("Loaded pre-infomap support from: ", preinfomap_file)

gc_support_preinfomap <- do.call("rbind", lapply(ti_preinfomap, function(x) {
  return(x[[1]])
})) |>
  dplyr::mutate(ti_uid = ti_pairs$ti_uid, .before = "target") |>
  dplyr::rename("GC_support_preinfomap" = preinfomap_GC_support)

truth_set <- cbind(
  ti_pairs |> dplyr::select(
    "ti_uid",
    "target",
    "indication_mesh_id",
    "indication_mesh_term",
    "combined_max_phase",
    "trait_gpmap",
    "trait_name_gpmap",
    gpmap_support = "GC_support",
    "coloc_groups",
    pairwisecoloc_support = "GC_support_colocpairs"
  ),
  gc_support_preinfomap |> dplyr::select(
    preinfomap_support = "GC_support_preinfomap"
  )
) |>
  dplyr::filter(combined_max_phase == "Launched")

out_file <- file.path(analysis_dir, "tipairs_launched_truthset.tsv")
data.table::fwrite(truth_set, out_file, sep = "\t")
message("Wrote: ", out_file)
