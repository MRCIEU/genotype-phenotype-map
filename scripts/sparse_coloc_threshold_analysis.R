source("../pipeline_steps/constants.R")


parser <- argparser::arg_parser(
  paste0(
    "Evaluate sparse finemapped coloc reliability: mask SNPs by reporting p-value by setting LBF=0 ",
    "(mimics released sumstats cutoffs; uses finemapped regional files with LBF only)."
  )
)
parser <- argparser::add_argument(parser, "--ld_block", help = "LD block identifier", type = "character")
parser <- argparser::add_argument(parser, "--sparse_study", help = "unique_study_id", type = "character")
parser <- argparser::add_argument(parser, "--output_file", type = "character", default = NA)
args <- argparser::parse_args(parser)
# args$sparse_study <- "godmc-methylation-cg00768179_EUR/1/7451118-9307591_1"

ensure_p_column <- function(gwas) {
  if ("P" %in% names(gwas)) {
    return(gwas)
  }
  if (all(c("BETA", "SE") %in% names(gwas))) {
    b <- as.numeric(gwas$BETA)
    s <- as.numeric(gwas$SE)
    z <- ifelse(!is.na(s) & s > 0, b / s, NA_real_)
    gwas$P <- 2 * stats::pnorm(-abs(z))
    return(gwas)
  }
  stop(
    "Finemapped GWAS must contain column P (or BETA and SE to derive P) for p-value masking. ",
    "Found columns: ", paste(names(gwas), collapse = ", ")
  )
}

#' Mimic \"only variants with P < Pthr released\": among originally real rows, set LBF = 0 where P >= p_reporting
mask_lbf_by_p_threshold <- function(gwas, real_indices, p_reporting) {
  out <- gwas
  pad_idx <- real_indices[out$P[real_indices] >= p_reporting]
  if (length(pad_idx) > 0) {
    out$LBF[pad_idx] <- 0
  }
  return(out)
}

lbf_vector_from_gwas <- function(gwas) {
  v <- as.numeric(gwas$LBF)
  names(v) <- gwas$SNP
  v[is.na(v)] <- 0
  return(v)
}

#' Read coloc_pairwise_results for one sparse study: partner unique_study_ids and file H4 per partner.
#' Returns list(h4_by_control, control_study_ids). Stops if file missing or no rows for sparse_uid.
pairwise_file_h4_lookup <- function(pairwise_path, ld_block, sparse_uid) {
  if (!file.exists(pairwise_path)) {
    stop(
      "Missing ", pairwise_path, ". ",
      "coloc_pairwise_results.tsv.gz is required to choose control studies.",
      call. = FALSE
    )
  }
  cw <- vroom::vroom(
    pairwise_path,
    show_col_types = FALSE,
    col_types = coloc_pairwise_results_column_types
  )
  if ("ld_block" %in% names(cw)) {
    cw <- dplyr::filter(cw, is.na(.data$ld_block) | .data$ld_block == ld_block)
  }
  h4_by_control <- cw |>
    dplyr::filter(.data$unique_study_a == sparse_uid | .data$unique_study_b == sparse_uid) |>
    dplyr::mutate(
      control_study = dplyr::if_else(
        .data$unique_study_a == sparse_uid,
        .data$unique_study_b,
        .data$unique_study_a
      ),
      pairwise_file_h4 = dplyr::coalesce(.data$PP.H4.abf, .data$h4)
    ) |>
    dplyr::filter(.data$control_study != sparse_uid) |>
    dplyr::group_by(.data$control_study) |>
    dplyr::arrange(dplyr::coalesce(.data$ignore, FALSE)) |>
    dplyr::slice_head(n = 1) |>
    dplyr::ungroup() |>
    dplyr::select(.data$control_study, .data$pairwise_file_h4)

  if (nrow(h4_by_control) == 0) {
    stop(
      "No coloc_pairwise rows for sparse study '", sparse_uid, "' in ", pairwise_path, ".",
      call. = FALSE
    )
  }

  return(list(
    h4_by_control = h4_by_control,
    control_study_ids = unique(h4_by_control$control_study)
  ))
}

main <- function() {
  ld_info <- ld_block_dirs(args$ld_block)

  finemapped_studies <- vroom::vroom(
    glue::glue("{ld_info$ld_block_data}/finemapped_studies.tsv"),
    show_col_types = FALSE,
    col_types = finemapped_column_types
  )

  sparse_finemapped <- dplyr::filter(finemapped_studies, .data$unique_study_id == args$sparse_study & !ignore)
  if (nrow(sparse_finemapped) == 0) {
    stop(glue::glue(
      "Sparse study '{args$sparse_study}' not found in finemapped_studies.tsv ",
      "or all rows marked ignore."
    ))
  }
  if (nrow(sparse_finemapped) > 1) {
    message(glue::glue(
      "Note: {nrow(sparse_finemapped)} finemapped rows for this study; using the first (",
      "{sparse_finemapped$unique_study_id[1]}). Pass a narrower query if needed."
    ))
  }
  sparse_finemapped <- sparse_finemapped[1, ]

  pairwise_path <- file.path(ld_info$ld_block_data, "coloc_pairwise_results.tsv.gz")
  pairwise_coloc <- pairwise_file_h4_lookup(
    pairwise_path, args$ld_block, sparse_finemapped$unique_study_id[1]
  )
  pairwise_h4_by_control <- pairwise_coloc$h4_by_control
  control_names <- pairwise_coloc$control_study_ids
  message(glue::glue(
    "Using {length(control_names)} control study id(s) from coloc_pairwise (partners of sparse study)."
  ))

  sparse_gwas <- vroom::vroom(sparse_finemapped$file, show_col_types = FALSE)
  if (!"LBF" %in% names(sparse_gwas) || !"SNP" %in% names(sparse_gwas)) {
    stop("Finemapped file must contain SNP and LBF columns: ", sparse_finemapped$file)
  }
  sparse_gwas <- ensure_p_column(sparse_gwas)

  real_indices <- which(!(sparse_gwas$BETA == 0 & sparse_gwas$SE == 1))
  if (length(real_indices) == 0) {
    real_indices <- seq_len(nrow(sparse_gwas))
    message("No BETA/SE padding pattern; treating all rows as eligible for p-value masking.")
  }
  padded_indices <- setdiff(seq_len(nrow(sparse_gwas)), real_indices)
  n_real <- length(real_indices)
  n_lbf0_baseline <- sum(sparse_gwas$LBF[real_indices] != 0, na.rm = TRUE)
  message(glue::glue(
    "{n_real} non-padded rows (by BETA/SE); {length(padded_indices)} padded; ",
    "{n_lbf0_baseline} of those with LBF != 0 at baseline"
  ))

  if (n_real == 0) stop("No rows to test")

  p_grid <- c(0.05, 1e-3, 1e-5, 5e-6, 5e-8)

  control_finemapped <- dplyr::filter(
    finemapped_studies,
    .data$unique_study_id %in% control_names & !ignore
  )
  control_finemapped <- dplyr::filter(control_finemapped, min_p <= .data$p_value_threshold)

  skipped_partners <- setdiff(control_names, control_finemapped$unique_study_id)
  if (length(skipped_partners) > 0) {
    message(glue::glue(
      "{length(skipped_partners)} coloc_pairwise partner(s) skipped."
    ))
  }

  if (nrow(control_finemapped) == 0) {
    stop(
      "No control finemapped studies match coloc_pairwise partners after filters. ",
      "Check finemapped_studies.tsv and p-value gates for this block.",
      call. = FALSE
    )
  }
  message(glue::glue(
    "Running coloc vs {nrow(control_finemapped)} finemapped control(s); ",
    "pipeline H4 for {nrow(pairwise_h4_by_control)} pair(s) in coloc_pairwise_results."
  ))

  control_lbfs <- lapply(seq_len(nrow(control_finemapped)), function(i) {
    gwas <- vroom::vroom(control_finemapped$file[i], show_col_types = FALSE)
    if (!"LBF" %in% colnames(gwas)) {
      return(NULL)
    }
    return(lbf_vector_from_gwas(gwas))
  })
  names(control_lbfs) <- control_finemapped$unique_study_id
  control_lbfs <- control_lbfs[!sapply(control_lbfs, is.null)]

  if (length(control_lbfs) == 0) stop("No control studies have LBF data")
  message(glue::glue("Loaded LBF vectors for {length(control_lbfs)} control studies"))

  count_signal_snps <- function(gwas, idx = real_indices) {
    return(sum(gwas$LBF[idx] != 0, na.rm = TRUE))
  }

  run_one_mask <- function(test_gwas, p_label, p_reporting) {
    n_lbf_nonzero <- count_signal_snps(test_gwas)

    test_lbf <- lbf_vector_from_gwas(test_gwas)

    return(
      lapply(names(control_lbfs), function(ctrl_id) {
        ctrl_lbf <- control_lbfs[[ctrl_id]]
        shared_snps <- sort(intersect(names(test_lbf), names(ctrl_lbf)))
        if (length(shared_snps) < 50) {
          return(data.frame(
            control_study = ctrl_id,
            sparse_unique_study_id = sparse_finemapped$unique_study_id,
            p_reporting_threshold = p_reporting,
            p_reporting_label = p_label,
            n_real_data_points = n_lbf_nonzero,
            n_shared_snps = length(shared_snps),
            PP.H0.abf = NA, PP.H1.abf = NA, PP.H2.abf = NA, PP.H3.abf = NA, PP.H4.abf = NA
          ))
        }

        result <- coloc::coloc.bf_bf(
          bf1 = test_lbf[shared_snps],
          bf2 = ctrl_lbf[shared_snps]
        )

        return(data.frame(
          control_study = ctrl_id,
          sparse_unique_study_id = sparse_finemapped$unique_study_id,
          p_reporting_threshold = p_reporting,
          p_reporting_label = p_label,
          n_real_data_points = n_lbf_nonzero,
          n_shared_snps = length(shared_snps),
          PP.H0.abf = result$summary$PP.H0.abf,
          PP.H1.abf = result$summary$PP.H1.abf,
          PP.H2.abf = result$summary$PP.H2.abf,
          PP.H3.abf = result$summary$PP.H3.abf,
          PP.H4.abf = result$summary$PP.H4.abf
        ))
      }) |> dplyr::bind_rows()
    )
  }

  message("  Baseline (finemapped LBF as-is)...")
  results <- run_one_mask(sparse_gwas, "baseline", NA_real_)

  for (p_thr in p_grid) {
    message(glue::glue("  Zero LBF where P >= {signif(p_thr, 4)} (non-padded rows only)..."))
    test_gwas <- mask_lbf_by_p_threshold(sparse_gwas, real_indices, p_thr)
    results <- dplyr::bind_rows(results, run_one_mask(test_gwas, as.character(p_thr), p_thr))
  }

  results$sparse_study <- args$sparse_study
  results$ld_block <- args$ld_block
  results$total_real_rows_nonpadded <- n_real
  results$n_lbf_nonzero_baseline <- n_lbf0_baseline
  results$passes_h4_threshold <- results$PP.H4.abf >= posterior_prob_h4_threshold

  results <- dplyr::left_join(results, pairwise_h4_by_control, by = "control_study")
  results$passes_h4_threshold_both <- !is.na(results$pairwise_file_h4) &
    results$PP.H4.abf >= posterior_prob_h4_threshold &
    results$pairwise_file_h4 >= posterior_prob_h4_threshold

  thr <- posterior_prob_h4_threshold
  summary_tbl <- results |>
    dplyr::filter(!is.na(PP.H4.abf)) |>
    dplyr::group_by(p_reporting_label, p_reporting_threshold) |>
    dplyr::summarise(
      n_real_data_points = dplyr::first(n_real_data_points),
      n_controls = dplyr::n(),
      mean_h4 = mean(PP.H4.abf),
      median_h4 = median(PP.H4.abf),
      max_h4 = max(PP.H4.abf),
      fraction_passing_h4 = mean(PP.H4.abf >= thr),
      n_passing_h4 = sum(PP.H4.abf >= thr),
      mean_h4_pairwise_file = {
        v <- pairwise_file_h4[!is.na(pairwise_file_h4)]
        if (length(v) == 0) NA_real_ else mean(v)
      },
      n_pairwise_file_matched = sum(!is.na(pairwise_file_h4)),
      fraction_pairwise_file_pass_h4 = {
        m <- sum(!is.na(pairwise_file_h4))
        if (m == 0) NA_real_ else sum(pairwise_file_h4 >= thr, na.rm = TRUE) / m
      },
      fraction_both_pass_h4 = {
        m <- sum(!is.na(pairwise_file_h4))
        if (m == 0) {
          NA_real_
        } else {
          sum(
            !is.na(pairwise_file_h4) & PP.H4.abf >= thr & pairwise_file_h4 >= thr,
            na.rm = TRUE
          ) / m
        }
      },
      .groups = "drop"
    )

  summary_tbl <- summary_tbl |>
    dplyr::mutate(sort_key = dplyr::coalesce(p_reporting_threshold, Inf)) |>
    dplyr::arrange(dplyr::desc(is.na(p_reporting_threshold)), dplyr::desc(sort_key)) |>
    dplyr::select(-sort_key)

  message("\n=== Summary: H4 by reporting p-threshold ===")
  print(as.data.frame(summary_tbl), row.names = FALSE)

  baseline_h4 <- summary_tbl$mean_h4[summary_tbl$p_reporting_label == "baseline"]
  if (length(baseline_h4) > 0 && !is.na(baseline_h4) && baseline_h4 > 0) {
    dropped <- summary_tbl |>
      dplyr::filter(!is.na(p_reporting_threshold), mean_h4 < baseline_h4 * 0.9)
    if (nrow(dropped) > 0) {
      worst <- dropped |> dplyr::slice_max(order_by = p_reporting_threshold, n = 1, with_ties = FALSE)
      message(glue::glue(
        "\nInflection (approx): mean H4 drops >10% from baseline ({round(baseline_h4, 3)}) ",
        "once you zero LBF for P >= {signif(worst$p_reporting_threshold, 4)} ",
        "(~{worst$n_real_data_points} SNPs with non-zero LBF; {worst$n_passing_h4}/{worst$n_controls} pairs pass H4)"
      ))
    } else {
      message("\nNo inflection detected — H4 stays within 10% of baseline across all p cutoffs")
    }
  }

  padding_only <- summary_tbl |>
    dplyr::filter(!is.na(p_reporting_threshold)) |>
    dplyr::slice_min(order_by = n_real_data_points, n = 1, with_ties = FALSE)
  if (nrow(padding_only) > 0) {
    r <- padding_only[1, ]
    message(glue::glue(
      "\nStrictest cutoff in grid: P >= {signif(r$p_reporting_threshold, 4)} → ",
      "{r$n_real_data_points} SNPs with LBF != 0; mean H4 = {round(r$mean_h4, 4)}; ",
      "{r$n_passing_h4}/{r$n_controls} pairs pass H4 > {posterior_prob_h4_threshold}"
    ))
  }

  output_file <- args$output_file
  if (is.na(output_file)) {
    safe_ld <- gsub("/", "_", args$ld_block)
    safe_sparse <- gsub("/", "_", args$sparse_study)
    output_file <- file.path(
      getwd(),
      glue::glue("sparse_coloc_analysis_{safe_ld}_{safe_sparse}.tsv")
    )
  }
  vroom::vroom_write(results, output_file)
  message(glue::glue("\nFull results written to {output_file}"))
  return(invisible(NULL))
}

main()
