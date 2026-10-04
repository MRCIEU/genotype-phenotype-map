standardise_step <- "standardise_studies_in_ld_block.R"

run_standardise <- function(paths, ...) {
  return(run_step(
    standardise_step,
    "--ld_block", fixture_block,
    "--completed_output_file", paths$standardisation_complete,
    ...
  ))
}

panel_snps <- read_tsv(ld_block_file_paths(fixture_block)$ld_matrix_tsv)$SNP
truth <- read_fixture_gwas("study_a_imputed")

test_that("standardise: alleles are flipped back to the reference panel orientation", {
  paths <- local_block()
  file <- place_study("study_a", "extracted", read_fixture_gwas("study_a_extracted"))
  write_index(study_index_row("study_a", file), "extracted")

  result <- run_standardise(paths)
  expect_step_succeeded(result, paths$standardisation_complete)

  standardised_studies <- read_tsv(paths$standardised_studies)
  expect_equal(nrow(standardised_studies), 1)
  expect_equal(standardised_studies$file, sub("extracted", "standardised", file))
  expect_equal(standardised_studies$snps_removed_by_reference_panel, 5)
  expect_false(standardised_studies$eaf_from_reference_panel)

  standardised <- read_tsv(standardised_studies$file)
  expect_true(all(standardised$SNP %in% panel_snps))
  expect_equal(nrow(standardised), length(panel_snps))
  expect_equal(standardised$SNP, panel_snps)

  compared <- dplyr::inner_join(standardised, truth, by = "SNP", suffix = c("", "_truth"))
  expect_equal(compared$BETA, compared$BETA_truth, tolerance = 1e-6)
  expect_equal(compared$EAF, compared$EAF_truth, tolerance = 1e-6)
  expect_equal(compared$Z, compared$BETA / compared$SE)
})

test_that("standardise: messy input (lowercase alleles, duplicates, LP, zero SE) is cleaned", {
  paths <- local_block()
  extracted <- read_fixture_gwas("study_a_extracted")
  extracted$EA[1:20] <- tolower(extracted$EA[1:20])
  extracted <- dplyr::bind_rows(extracted, extracted[100:109, ])
  zero_se_bp <- extracted$BP[200]
  zero_se_and_beta_bp <- extracted$BP[201]
  expected_se <- extracted$SE[200]
  extracted$SE[200:201] <- 0
  extracted$BETA[201] <- 0
  extracted <- dplyr::mutate(extracted, LP = -log10(P)) |> dplyr::select(-P)

  file <- place_study("study_messy", "extracted", extracted)
  write_index(study_index_row("study_messy", file), "extracted")

  result <- run_standardise(paths)
  expect_step_succeeded(result, paths$standardisation_complete)

  standardised <- read_tsv(sub("extracted", "standardised", file))
  expect_equal(nrow(standardised), length(panel_snps))
  expect_equal(anyDuplicated(standardised$SNP), 0)
  expect_true(all(standardised$EA == toupper(standardised$EA)))
  expect_true("P" %in% colnames(standardised))
  expect_false("LP" %in% colnames(standardised))
  expect_true(all(standardised$P >= 0 & standardised$P <= 1))

  # zero SE is back-derived from BETA and P, unless BETA is 0, where a small sentinel SE is used
  expect_equal(standardised$SE[standardised$BP == zero_se_bp], expected_se, tolerance = 1e-4)
  expect_equal(standardised$SE[standardised$BP == zero_se_and_beta_bp], 0.00001)
  expect_true(all(is.finite(standardised$Z)))
})

test_that("standardise: EAF is taken from the reference panel when it is entirely missing", {
  paths <- local_block()
  extracted <- read_fixture_gwas("study_a_extracted") |> dplyr::mutate(EAF = NA)
  file <- place_study("study_no_eaf", "extracted", extracted)
  write_index(study_index_row("study_no_eaf", file), "extracted")

  result <- run_standardise(paths)
  expect_step_succeeded(result, paths$standardisation_complete)

  expect_true(read_tsv(paths$standardised_studies)$eaf_from_reference_panel)
  standardised <- read_tsv(sub("extracted", "standardised", file))
  expect_false(any(is.na(standardised$EAF)))
  compared <- dplyr::inner_join(standardised, truth, by = "SNP", suffix = c("", "_truth"))
  expect_equal(compared$EAF, compared$EAF_truth, tolerance = 1e-6)
})

test_that("standardise: variants with extreme EAF are removed", {
  paths <- local_block()
  extracted <- read_fixture_gwas("study_a_extracted")
  extreme_bps <- extracted$BP[c(300, 301, 302)]
  extracted$EAF[300] <- 0.001
  extracted$EAF[301] <- 0.999
  extracted$EAF[302] <- 0.005
  file <- place_study("study_extreme_eaf", "extracted", extracted)
  write_index(study_index_row("study_extreme_eaf", file), "extracted")

  result <- run_standardise(paths)
  expect_step_succeeded(result, paths$standardisation_complete)

  standardised <- read_tsv(sub("extracted", "standardised", file))
  expect_false(any(standardised$BP %in% extreme_bps))
  expect_true(all(standardised$EAF > 0.005 & standardised$EAF < 0.995))
})

test_that("standardise: invalid P or SE values fail the step without writing the sentinel", {
  for (bad_column in c("P", "SE")) {
    paths <- local_block()
    extracted <- read_fixture_gwas("study_a_extracted")
    if (bad_column == "P") extracted$P[10] <- 1.5 else extracted$SE[10] <- -0.1
    file <- place_study("study_invalid", "extracted", extracted)
    write_index(study_index_row("study_invalid", file), "extracted")

    result <- run_standardise(paths)
    expect_step_failed(result, paths$standardisation_complete, glue::glue("GWAS has some {bad_column} values outside"))
  }
})

test_that("standardise: rare variant studies bypass the reference panel and EAF filters", {
  paths <- local_block()
  extracted <- read_fixture_gwas("study_a_extracted")
  extracted$EAF[300] <- 0.0001
  file <- place_study("study_rare", "extracted", extracted)
  write_index(study_index_row("study_rare", file, variant_type = variant_types$rare_exome), "extracted")

  result <- run_standardise(paths)
  expect_step_succeeded(result, paths$standardisation_complete)

  standardised <- read_tsv(sub("extracted", "standardised", file))
  expect_equal(nrow(standardised), nrow(extracted))
  expect_true(any(!standardised$SNP %in% panel_snps))
  expect_equal(read_tsv(paths$standardised_studies)$snps_removed_by_reference_panel, 0)
})

test_that("standardise: small extractions are skipped by coverage rules", {
  paths <- local_block()
  extracted <- read_fixture_gwas("study_a_extracted")
  in_panel <- dplyr::filter(extracted, !startsWith(RSID, "rs_not_in_panel"))

  dense_small <- place_study("study_dense_small", "extracted", in_panel[1:100, ])
  sparse_tiny <- place_study("study_sparse_tiny", "extracted", in_panel[1:3, ])
  sparse_ok <- place_study("study_sparse_ok", "extracted", in_panel[1:10, ])
  write_index(dplyr::bind_rows(
    study_index_row("study_dense_small", dense_small),
    study_index_row("study_sparse_tiny", sparse_tiny, coverage = coverage_types$sparse),
    study_index_row("study_sparse_ok", sparse_ok, coverage = coverage_types$sparse)
  ), "extracted")

  result <- run_standardise(paths)
  expect_step_succeeded(result, paths$standardisation_complete)

  standardised_studies <- read_tsv(paths$standardised_studies)
  expect_equal(standardised_studies$study, "study_sparse_ok")

  skipped <- read_tsv(paths$standardised_skipped) |> dplyr::arrange(study)
  expect_equal(skipped$study, c("study_dense_small", "study_sparse_tiny"))
  expect_equal(skipped$reason, c("below_min_dense", "below_min_sparse"))
  expect_equal(skipped$n_variants, c(100, 3))
})

test_that("standardise: rerunning is a no-op, but a stale skip rule_version is reprocessed", {
  paths <- local_block()
  file <- place_study("study_a", "extracted", read_fixture_gwas("study_a_extracted"))
  write_index(study_index_row("study_a", file), "extracted")
  vroom::vroom_write(data.frame(
    study = "study_a", input_file = file, reason = "below_min_dense", n_variants = 10L,
    coverage = coverage_types$dense, variant_type = variant_types$common, ld_block = fixture_block,
    rule_version = "min_dense=1000,min_sparse=4"
  ), paths$standardised_skipped)

  first_run <- run_standardise(paths)
  expect_step_succeeded(first_run, paths$standardisation_complete)
  expect_equal(read_tsv(paths$standardised_studies)$study, "study_a")

  unlink(paths$standardisation_complete)
  second_run <- run_standardise(paths)
  expect_step_succeeded(second_run, paths$standardisation_complete)
  expect_true(any(grepl("already standardised or skipped", second_run$output)))
  expect_equal(nrow(read_tsv(paths$standardised_studies)), 1)
})

test_that("standardise: --block_list excludes matching studies", {
  paths <- local_block()
  file_a <- place_study("study_a", "extracted", read_fixture_gwas("study_a_extracted"))
  file_blocked <- place_study("blocked_study", "extracted", read_fixture_gwas("study_a_extracted"))
  write_index(dplyr::bind_rows(
    study_index_row("study_a", file_a),
    study_index_row("blocked_study", file_blocked)
  ), "extracted")
  block_list_file <- file.path(pipeline_metadata_dir, "unit_block_list.csv")
  vroom::vroom_write(data.frame(id_pattern = "blocked*", cis_trans = NA), block_list_file, delim = ",")

  result <- run_standardise(paths, "--block_list", block_list_file)
  expect_step_succeeded(result, paths$standardisation_complete)
  expect_equal(read_tsv(paths$standardised_studies)$study, "study_a")
})

test_that("standardise: a block without extracted studies writes empty outputs", {
  paths <- local_block()

  result <- run_standardise(paths)
  expect_step_succeeded(result, paths$standardisation_complete)
  expect_equal(nrow(read_tsv(paths$standardised_studies)), 0)
  expect_equal(nrow(read_tsv(paths$standardised_skipped)), 0)
})
