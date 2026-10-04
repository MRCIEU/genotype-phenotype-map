impute_step <- "impute_studies_in_ld_block.R"

run_impute <- function(paths, ...) {
  return(run_step(
    impute_step,
    "--ld_block", fixture_block,
    "--completed_output_file", paths$imputation_complete,
    ...
  ))
}

block_paths <- ld_block_file_paths(fixture_block)
panel <- read_tsv(block_paths$ld_matrix_tsv)
truth <- read_fixture_gwas("study_a_imputed")
standardised_a <- read_fixture_gwas("study_a_standardised")
dropped_snps <- read_tsv(fixture_path("gwas", "study_a_standardised_dropped_snps.tsv"), delim = "\t")$SNP

place_standardised <- function(study, gwas, ...) {
  file <- place_study(study, "standardised", gwas)
  return(study_index_row(study, file, ...))
}

test_that("impute: missing SNPs are imputed close to their true values", {
  paths <- local_block()
  write_index(place_standardised("study_a", standardised_a), "standardised")

  result <- run_impute(paths)
  expect_step_succeeded(result, paths$imputation_complete)

  imputed_studies <- read_tsv(paths$imputed_studies)
  expect_equal(nrow(imputed_studies), 1)
  expect_equal(imputed_studies$rows_imputed, length(dropped_snps))
  expect_gte(imputed_studies$b_cor, 0.7)

  imputed <- read_tsv(imputed_studies$file)
  expect_false(any(is.na(imputed$BETA) | is.na(imputed$SE) | is.na(imputed$P)))
  expect_true(all(imputed$SE > 0 & imputed$P >= 0 & imputed$P <= 1))
  # every dropped SNP is imputed; a few observed SNPs can also be re-imputed as SE outliers
  expect_true(all(imputed$IMPUTED[imputed$SNP %in% dropped_snps]))
  expect_lt(sum(imputed$IMPUTED & !imputed$SNP %in% dropped_snps), 10)
  # results are trimmed to inside the BP range of the original (standardised) GWAS
  expect_true(all(imputed$BP >= min(standardised_a$BP) & imputed$BP <= max(standardised_a$BP)))

  compared <- dplyr::inner_join(
    dplyr::filter(imputed, IMPUTED),
    truth,
    by = "SNP",
    suffix = c("", "_truth")
  )
  expect_gt(nrow(compared), 0.9 * length(dropped_snps))
  expect_gt(cor(compared$Z, compared$Z_truth), 0.7)
})

test_that("impute: no imputed SNP is much more significant than the observed SNPs it is in LD with", {
  paths <- local_block()
  write_index(place_standardised("study_a", standardised_a), "standardised")
  result <- run_impute(paths)
  expect_step_succeeded(result, paths$imputation_complete)

  imputed <- read_tsv(read_tsv(paths$imputed_studies)$file)
  ld_matrix <- as.matrix(read_tsv(block_paths$ld_matrix_vcor, col_names = FALSE, altrep = FALSE))
  ld_index <- match(imputed$SNP, panel$SNP)
  min_observed_p <- min(imputed$P[!imputed$IMPUTED])
  suspicious <- which(imputed$IMPUTED & imputed$P < min(min_observed_p, lowest_p_value_threshold))

  for (row in suspicious) {
    in_ld <- ld_index[!imputed$IMPUTED] %in% which(ld_matrix[ld_index[row], ] > 0.6)
    correlated_p <- imputed$P[!imputed$IMPUTED][in_ld]
    expect_true(length(correlated_p) > 0 && min(correlated_p) * 0.1 <= imputed$P[row], info = imputed$SNP[row])
  }
})

test_that("impute: a GWAS with nothing missing falls back to the standardised data", {
  paths <- local_block()
  complete_gwas <- dplyr::select(truth, -IMPUTED)
  write_index(place_standardised("study_complete", complete_gwas), "standardised")

  result <- run_impute(paths)
  expect_step_succeeded(result, paths$imputation_complete)

  imputed_studies <- read_tsv(paths$imputed_studies)
  expect_equal(imputed_studies$rows_imputed, 0)
  expect_true(is.na(imputed_studies$b_cor))
  imputed <- read_tsv(imputed_studies$file)
  expect_equal(nrow(imputed), nrow(complete_gwas))
  expect_false(any(imputed$IMPUTED))
  expect_equal(imputed$BETA, complete_gwas$BETA, tolerance = 1e-8)
})

test_that("impute: a constant Z score cannot be imputed, so missing values are padded", {
  paths <- local_block()
  constant_z <- dplyr::mutate(standardised_a, BETA = 2 * SE, Z = 2, P = 2 * pnorm(-2))
  write_index(place_standardised("study_constant_z", constant_z), "standardised")

  result <- run_impute(paths)
  expect_step_succeeded(result, paths$imputation_complete)

  imputed_studies <- read_tsv(paths$imputed_studies)
  expect_equal(imputed_studies$rows_imputed, 0)
  imputed <- read_tsv(imputed_studies$file)
  padded <- dplyr::filter(imputed, SNP %in% dropped_snps)
  expect_equal(nrow(padded), length(dropped_snps))
  expect_true(all(padded$IMPUTED & padded$BETA == 0 & padded$SE == 1 & padded$P == 1))
})

test_that("impute: sparse coverage studies are padded instead of imputed", {
  paths <- local_block()
  write_index(place_standardised("study_sparse", standardised_a, coverage = coverage_types$sparse), "standardised")

  result <- run_impute(paths)
  expect_step_succeeded(result, paths$imputation_complete)

  imputed_studies <- read_tsv(paths$imputed_studies)
  expect_equal(imputed_studies$rows_imputed, 0)
  expect_true(is.na(imputed_studies$significant_rows_imputed))
  imputed <- read_tsv(imputed_studies$file)
  expect_equal(nrow(imputed), nrow(panel))
  padded <- dplyr::filter(imputed, SNP %in% dropped_snps)
  expect_true(all(padded$IMPUTED & padded$BETA == 0 & padded$SE == 1 & padded$P == 1 & padded$Z == 0))
  expect_false(any(imputed$IMPUTED[!imputed$SNP %in% dropped_snps]))
})

test_that("impute: an EAF of 0 fails the step", {
  paths <- local_block()
  bad_eaf <- standardised_a
  bad_eaf$EAF[5] <- 0
  write_index(place_standardised("study_bad_eaf", bad_eaf), "standardised")

  result <- run_impute(paths)
  expect_step_failed(result, paths$imputation_complete, "funky EAF")
})

test_that("impute: only common variant studies are imputed, and reruns are skipped", {
  paths <- local_block()
  write_index(dplyr::bind_rows(
    place_standardised("study_a", standardised_a),
    place_standardised("study_rare", standardised_a, variant_type = variant_types$rare_exome)
  ), "standardised")

  first_run <- run_impute(paths)
  expect_step_succeeded(first_run, paths$imputation_complete)
  expect_equal(read_tsv(paths$imputed_studies)$study, "study_a")

  unlink(paths$imputation_complete)
  second_run <- run_impute(paths)
  expect_step_succeeded(second_run, paths$imputation_complete)
  expect_true(any(grepl("already imputed", second_run$output)))
  expect_equal(nrow(read_tsv(paths$imputed_studies)), 1)
})

test_that("impute: --block_list excludes matching studies", {
  paths <- local_block()
  write_index(dplyr::bind_rows(
    place_standardised("study_a", standardised_a),
    place_standardised("blocked_study", standardised_a)
  ), "standardised")
  block_list_file <- file.path(pipeline_metadata_dir, "unit_block_list.csv")
  vroom::vroom_write(data.frame(id_pattern = "blocked*", cis_trans = NA), block_list_file, delim = ",")

  result <- run_impute(paths, "--block_list", block_list_file)
  expect_step_succeeded(result, paths$imputation_complete)
  expect_equal(read_tsv(paths$imputed_studies)$study, "study_a")
})

test_that("impute: a block without standardised studies writes an empty index", {
  paths <- local_block()

  result <- run_impute(paths)
  expect_step_succeeded(result, paths$imputation_complete)
  expect_equal(nrow(read_tsv(paths$imputed_studies)), 0)
})
