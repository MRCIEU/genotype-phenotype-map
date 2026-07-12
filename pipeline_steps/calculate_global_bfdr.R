source("constants.R")

# Calculate a single global Bayesian FDR (BFDR) threshold across the pairwise coloc
# results of every LD block, and store it in the pipeline_metadata directory so the
# downstream clustering step can apply one consistent H4 threshold everywhere.

parser <- argparser::arg_parser("Calculate the global Bayesian FDR threshold across all LD blocks")
parser <- argparser::add_argument(
  parser,
  "--output_file",
  help = "File to write the global BFDR threshold and metadata to",
  type = "character",
  default = NA
)
parser <- argparser::add_argument(
  parser,
  "--target_bfdr",
  help = "Target Bayesian FDR",
  type = "numeric",
  default = target_global_bfdr
)
parser <- argparser::add_argument(
  parser,
  "--block_list",
  help = "CSV of studies excluded from colocalisation (only used for output naming parity)",
  type = "character",
  default = NA
)
parser <- argparser::add_argument(
  parser,
  "--completed_output_file",
  help = "Sentinel file to write once the global BFDR threshold has been calculated",
  type = "character",
  default = NA
)

args <- argparser::parse_args(parser)

main <- function() {
  start_time <- Sys.time()

  output_file <- args$output_file
  if (is.null(output_file) || is.na(output_file)) {
    output_file <- global_bfdr_file_path(args$block_list)
  }

  pairwise_files <- Sys.glob(glue::glue("{ld_block_data_dir}*/*/*/coloc_pairwise_results.tsv.gz"))
  message(glue::glue("Found {length(pairwise_files)} LD blocks with pairwise coloc results"))

  h4_values <- unlist(lapply(pairwise_files, read_block_h4), use.names = FALSE)
  message(glue::glue(
    "Collected {length(h4_values)} pairwise H4 values in {diff_time_taken(start_time)}"
  ))

  bfdr <- calculate_global_bfdr_threshold(
    posterior_probs = h4_values,
    target_bfdr = args$target_bfdr,
    minimum_threshold = posterior_prob_threshold_minimum
  )

  global_bfdr <- data.frame(
    threshold = bfdr$threshold,
    estimated_bfdr = bfdr$estimated_bfdr,
    n_discoveries = bfdr$n_discoveries,
    target_bfdr = args$target_bfdr,
    minimum_threshold = posterior_prob_threshold_minimum,
    n_pairs = length(h4_values),
    n_blocks = length(pairwise_files),
    default_threshold = posterior_prob_h4_threshold
  )

  message(glue::glue(
    "Global BFDR threshold: {signif(bfdr$threshold, 3)} ",
    "(default {signif(posterior_prob_h4_threshold, 3)}, target BFDR={args$target_bfdr}, ",
    "discoveries={bfdr$n_discoveries}, estimated BFDR={signif(bfdr$estimated_bfdr, 3)})"
  ))

  dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
  vroom::vroom_write(global_bfdr, output_file)
  message(glue::glue("Wrote global BFDR threshold to {output_file}"))

  completed_output_file <- args$completed_output_file
  if (!is.null(completed_output_file) && !is.na(completed_output_file)) {
    dir.create(dirname(completed_output_file), showWarnings = FALSE, recursive = TRUE)
    vroom::vroom_write(data.frame(), completed_output_file)
    message(glue::glue("Wrote completion sentinel to {completed_output_file}"))
  }

  return()
}

calculate_global_bfdr_threshold <- function(posterior_probs, target_bfdr = 0.05, minimum_threshold = 0.5) {
  valid_probs <- posterior_probs[is.finite(posterior_probs)]
  valid_probs <- valid_probs[!is.na(valid_probs)]
  if (length(valid_probs) == 0) {
    return(list(
      threshold = 1,
      estimated_bfdr = NA_real_,
      n_discoveries = 0
    ))
  }

  ordered_probs <- sort(valid_probs, decreasing = TRUE)
  cumulative_bfdr <- cumsum(1 - ordered_probs) / seq_along(ordered_probs)
  passing_indices <- which(cumulative_bfdr <= target_bfdr)

  if (length(passing_indices) == 0) {
    return(list(
      threshold = 1,
      estimated_bfdr = NA_real_,
      n_discoveries = 0
    ))
  }

  threshold <- ordered_probs[max(passing_indices)]
  threshold <- max(threshold, minimum_threshold)
  discoveries <- valid_probs[valid_probs >= threshold]
  estimated_bfdr <- if (length(discoveries) == 0) NA_real_ else sum(1 - discoveries) / length(discoveries)

  return(list(
    threshold = threshold,
    estimated_bfdr = estimated_bfdr,
    n_discoveries = length(discoveries)
  ))
}

#' Read the eligible H4 values from a block's pairwise coloc results.
#' Self-pairs (study_a == study_b) and ignored pairs are excluded from the global pool,
#' since self-colocs are ~1 by construction and would bias the FDR threshold downward.
read_block_h4 <- function(pairwise_file) {
  if (!file.exists(pairwise_file) || file.size(pairwise_file) == 0) {
    return(numeric(0))
  }

  header <- tryCatch(
    names(data.table::fread(pairwise_file, nrows = 0, showProgress = FALSE)),
    error = function(e) character(0)
  )
  h4_column <- if ("h4" %in% header) "h4" else if ("PP.H4.abf" %in% header) "PP.H4.abf" else NA_character_
  if (is.na(h4_column)) {
    return(numeric(0))
  }

  select_columns <- intersect(c(h4_column, "study_a", "study_b", "ignore"), header)
  pairs <- tryCatch(
    data.table::fread(
      pairwise_file,
      select = select_columns,
      showProgress = FALSE,
      verbose = FALSE
    ),
    error = function(e) {
      message(glue::glue("Could not read {pairwise_file}: {e$message}"))
      return(NULL)
    }
  )
  if (is.null(pairs) || nrow(pairs) == 0) {
    return(numeric(0))
  }

  keep <- rep(TRUE, nrow(pairs))
  if (all(c("study_a", "study_b") %in% names(pairs))) {
    keep <- keep & !(!is.na(pairs$study_a) & !is.na(pairs$study_b) & pairs$study_a == pairs$study_b)
  }
  if ("ignore" %in% names(pairs)) {
    keep <- keep & !(pairs$ignore %in% TRUE)
  }

  return(as.numeric(pairs[[h4_column]][keep]))
}

invisible(main())
