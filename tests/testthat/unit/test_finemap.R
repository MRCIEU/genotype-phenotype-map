finemap_step <- "finemap_studies_in_ld_block.R"

run_finemap <- function(paths, ...) {
  return(run_step(
    finemap_step,
    "--ld_block", fixture_block,
    "--completed_output_file", paths$finemapping_complete,
    ...
  ))
}

place_imputed <- function(study, gwas, ...) {
  file <- place_study(study, "imputed", gwas)
  return(study_index_row(study, file, ...))
}

read_finemapped <- function(paths) {
  return(read_tsv(paths$finemapped_studies, col_types = finemapped_column_types) |> dplyr::arrange(unique_study_id))
}

causal_snp <- function(name) {
  return(causal_snps$snp[causal_snps$name == name])
}

causal_bp <- function(name) {
  return(causal_snps$bp[causal_snps$name == name])
}

fixture_leads <- read_tsv(fixture_path("finemapped", "finemapped_studies_template.tsv"))
study_a <- read_fixture_gwas("study_a_imputed")

test_that("finemap: a GWAS with 2 signals is split into 2 credible sets", {
  paths <- local_block()
  write_index(place_imputed("study_a", study_a), "imputed")

  result <- run_finemap(paths)
  expect_step_succeeded(result, paths$finemapping_complete)

  finemapped <- read_finemapped(paths)
  expected_leads <- dplyr::filter(fixture_leads, study == "study_a")
  expect_equal(nrow(finemapped), 2)
  expect_equal(finemapped$unique_study_id, paste0("study_a_", fixture_block, "_", 1:2))
  expect_true(all(finemapped$finemap_message == "success"))
  expect_true(all(finemapped$first_finemap_num_results == 2))
  expect_false(any(finemapped$qc_step_run))
  expect_equal(sort(finemapped$snp), sort(expected_leads$snp))
  expect_equal(sort(finemapped$snp), sort(c(causal_snp("causal_1"), causal_snp("causal_2"))))
  expect_equal(
    finemapped$min_p[order(finemapped$snp)],
    expected_leads$min_p[order(expected_leads$snp)],
    tolerance = 1e-6
  )

  # SNPs whose LBF implies a much lower P value than their actual P value are logged and removed per credible set.
  # This happens even for this cleanly simulated GWAS
  problematic_file <- file.path(paths$ld_block_data, "problematic_finemapped_snps.tsv")
  problematic <- data.frame(unique_study_id = character())
  if (file.exists(problematic_file)) problematic <- read_tsv(problematic_file)
  for (i in seq_len(nrow(finemapped))) {
    credible_set <- read_tsv(finemapped$file[i])
    expect_equal(colnames(credible_set)[1:8], c("SNP", "CHR", "BP", "BETA", "SE", "EAF", "IMPUTED", "LBF"))
    removed <- sum(problematic$unique_study_id == finemapped$unique_study_id[i])
    expect_equal(nrow(credible_set), nrow(study_a) - removed)
    expect_true(file.exists(finemapped$svg_file[i]))
  }
  with_lbfs <- read_tsv(unique(finemapped$file_with_lbfs))
  expect_true(all(c("LBF_1", "LBF_2", "EA", "OA", "Z", "P") %in% colnames(with_lbfs)))
  expect_true(file.exists(sub("_1.tsv.gz$", "_results.rds", finemapped$file[1])))
})

test_that("finemap: a single signal is kept unsplit with its credible set's lead SNP", {
  paths <- local_block()
  write_index(place_imputed("study_b", read_fixture_gwas("study_b_imputed")), "imputed")

  result <- run_finemap(paths)
  expect_step_succeeded(result, paths$finemapping_complete)

  finemapped <- read_finemapped(paths)
  expect_equal(nrow(finemapped), 1)
  expect_equal(finemapped$finemap_message, "less_than_2_cs")
  expect_equal(finemapped$snp, causal_snp("causal_1"))
  expect_equal(finemapped$bp, causal_bp("causal_1"))
  expect_equal(finemapped$first_finemap_num_results, 0)
  expect_true("LBF_1" %in% colnames(read_tsv(finemapped$file_with_lbfs)))
})

test_that("finemap: a GWAS with no signal uses its lowest P value SNP", {
  paths <- local_block()
  null_gwas <- read_fixture_gwas("study_null_imputed")
  write_index(place_imputed("study_null", null_gwas), "imputed")

  result <- run_finemap(paths)
  expect_step_succeeded(result, paths$finemapping_complete)

  finemapped <- read_finemapped(paths)
  min_p_row <- which.min(null_gwas$P)
  expect_equal(nrow(finemapped), 1)
  # SuSiE converges without a credible set, which is reported the same as finding 1
  expect_equal(finemapped$finemap_message, "less_than_2_cs")
  expect_equal(finemapped$bp, null_gwas$BP[min_p_row])
  expect_equal(finemapped$min_p, null_gwas$P[min_p_row])
})

test_that("finemap: small, tiny and sparse GWASes are not finemapped", {
  paths <- local_block()
  small <- study_a[seq(1, 800, by = 2), ]
  annotated_bp <- small$BP[10]

  small_unannotated <- small
  underflow_rows <- c(20, 30)
  small_unannotated$P[underflow_rows] <- 0
  small_unannotated$Z[underflow_rows] <- c(40, -45)

  write_index(dplyr::bind_rows(
    place_imputed("study_small", small, bp = annotated_bp),
    place_imputed("study_small_unannotated", small_unannotated, bp = 5),
    place_imputed("study_tiny", study_a[1:100, ]),
    place_imputed("study_sparse", study_a, coverage = coverage_types$sparse)
  ), "imputed")

  result <- run_finemap(paths)
  expect_step_succeeded(result, paths$finemapping_complete)

  finemapped <- read_finemapped(paths) |> dplyr::arrange(study)
  expect_equal(finemapped$study, c("study_small", "study_small_unannotated", "study_sparse"))
  expect_equal(finemapped$finemap_message, c("too_small_to_finemap", "too_small_to_finemap", "sparse_population"))
  expect_true(is.na(finemapped$first_finemap_num_results[finemapped$study == "study_sparse"]))

  # the SNP comes from the variant annotations at the study's bp, if there is one
  expect_equal(finemapped$snp[1], small$SNP[10])
  # otherwise it is the min P SNP, tie-broken on |Z| when P values underflow to 0
  expect_equal(finemapped$snp[2], small_unannotated$SNP[underflow_rows[2]])
})

test_that("finemap: 4 or more credible sets run the DENTIST QC step", {
  paths <- local_block()
  write_index(place_imputed("study_d", read_fixture_gwas("study_d_imputed")), "imputed")

  result <- run_finemap(paths)
  expect_step_succeeded(result, paths$finemapping_complete)
  expect_true(any(grepl("performing qc", result$output)))

  finemapped <- read_finemapped(paths)
  expect_gte(nrow(finemapped), 1)
  expect_true(all(finemapped$qc_step_run))
  expect_true(all(finemapped$first_finemap_num_results >= 4))
  expect_false(any(is.na(finemapped$snps_removed_by_qc)))
})

test_that("finemap: SNPs whose LBF disagrees with their P value are not chosen as lead SNPs", {
  paths <- local_block()
  inconsistent <- study_a
  inconsistent$P[inconsistent$SNP == causal_snp("causal_1")] <- 0.5
  write_index(place_imputed("study_a", inconsistent), "imputed")

  result <- run_finemap(paths)
  expect_step_succeeded(result, paths$finemapping_complete)

  finemapped <- read_finemapped(paths)
  expect_equal(nrow(finemapped), 2)
  expect_false(causal_snp("causal_1") %in% finemapped$snp)

  problematic_file <- file.path(paths$ld_block_data, "problematic_finemapped_snps.tsv")
  expect_true(file.exists(problematic_file))
  problematic <- read_tsv(problematic_file)
  expect_true(causal_snp("causal_1") %in% problematic$SNP)
  expect_true(all(problematic$unique_study_id %in% finemapped$unique_study_id))
})

test_that("finemap: cis studies are relabelled trans when the credible set lead is over 1Mb away", {
  paths <- local_block()
  write_index(dplyr::bind_rows(
    place_imputed("study_near", study_a, bp = causal_bp("causal_1"), cis_trans = cis_trans$cis_only),
    place_imputed("study_far", study_a, bp = 5000000, cis_trans = cis_trans$cis_only)
  ), "imputed")

  result <- run_finemap(paths)
  expect_step_succeeded(result, paths$finemapping_complete)

  finemapped <- read_finemapped(paths)
  expect_true(all(finemapped$cis_trans[finemapped$study == "study_near"] == cis_trans$cis_only))
  expect_true(all(finemapped$cis_trans[finemapped$study == "study_far"] == cis_trans$trans_only))
})

test_that("finemap: each credible set gets its own cis/trans label", {
  paths <- local_block()
  # causal_1 is just over 1Mb from this bp, causal_2 is well within 1Mb of it
  mixed_bp <- causal_bp("causal_1") + 1000001
  write_index(place_imputed("study_mixed", study_a, bp = mixed_bp, cis_trans = cis_trans$cis_only), "imputed")

  result <- run_finemap(paths)
  expect_step_succeeded(result, paths$finemapping_complete)

  finemapped <- read_finemapped(paths)
  expect_equal(finemapped$cis_trans[finemapped$snp == causal_snp("causal_1")], cis_trans$trans_only)
  expect_equal(finemapped$cis_trans[finemapped$snp == causal_snp("causal_2")], cis_trans$cis_only)
})

test_that("finemap: a GWAS that does not match the LD matrix fails the step", {
  paths <- local_block()
  duplicated_rows <- dplyr::bind_rows(study_a, study_a[1:5, ]) |> dplyr::arrange(BP)
  write_index(place_imputed("study_duplicated", duplicated_rows), "imputed")

  result <- run_finemap(paths)
  expect_step_failed(result, paths$finemapping_complete, "should match size")
})

test_that("finemap: a missing BETA fails the step", {
  paths <- local_block()
  missing_beta <- study_a
  missing_beta$BETA[100] <- NA
  write_index(place_imputed("study_missing_beta", missing_beta), "imputed")

  result <- run_finemap(paths)
  expect_step_failed(result, paths$finemapping_complete, "Missing SNP, BETA, or SE")
})

test_that("finemap: rare studies, blocked studies and already finemapped studies are skipped", {
  paths <- local_block()
  write_index(dplyr::bind_rows(
    place_imputed("study_b", read_fixture_gwas("study_b_imputed")),
    place_imputed("study_rare", study_a, variant_type = variant_types$rare_exome),
    place_imputed("blocked_study", study_a)
  ), "imputed")
  block_list_file <- file.path(pipeline_metadata_dir, "unit_block_list.csv")
  vroom::vroom_write(data.frame(id_pattern = "blocked*", cis_trans = NA), block_list_file, delim = ",")

  first_run <- run_finemap(paths, "--block_list", block_list_file)
  expect_step_succeeded(first_run, paths$finemapping_complete)
  expect_equal(read_finemapped(paths)$study, "study_b")

  unlink(paths$finemapping_complete)
  second_run <- run_finemap(paths, "--block_list", block_list_file)
  expect_step_succeeded(second_run, paths$finemapping_complete)
  expect_true(any(grepl("already finemapped", second_run$output)))
  expect_equal(read_finemapped(paths)$study, "study_b")
})

test_that("finemap: a block without imputed studies only writes the sentinel", {
  paths <- local_block()

  result <- run_finemap(paths)
  expect_step_succeeded(result, paths$finemapping_complete)
  expect_false(file.exists(paths$finemapped_studies))
})
