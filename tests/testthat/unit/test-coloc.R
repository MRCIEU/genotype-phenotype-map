coloc_step <- "coloc_studies_in_ld_block.R"

run_coloc <- function(paths, ...) {
  return(run_step(
    coloc_step,
    "--ld_block", fixture_block,
    "--completed_output_file", paths$coloc_complete,
    ...
  ))
}

read_coloc <- function(paths) {
  return(read_tsv(paths$coloc_pairwise) |> dplyr::arrange(unique_study_a, unique_study_b))
}

unique_id <- function(study, cs = 1) {
  return(paste0(study, "_", fixture_block, "_", cs))
}

find_pair <- function(coloc_results, id_a, id_b) {
  return(dplyr::filter(
    coloc_results,
    (unique_study_a == id_a & unique_study_b == id_b) | (unique_study_a == id_b & unique_study_b == id_a)
  ))
}

test_that("coloc: pairs within 50kb (or from the same study) are colocalised", {
  paths <- local_block()
  write_index(finemapped_fixture_rows(), "finemapped")

  result <- run_coloc(paths)
  expect_step_succeeded(result, paths$coloc_complete)

  coloc_results <- read_coloc(paths)
  # study_a's second signal is >50kb from study_b and study_c, so is only compared with study_a's first signal
  expect_equal(nrow(coloc_results), 4)
  expect_equal(nrow(find_pair(coloc_results, unique_id("study_a", 2), unique_id("study_b"))), 0)
  expect_equal(nrow(find_pair(coloc_results, unique_id("study_a", 2), unique_id("study_c"))), 0)
  expect_equal(nrow(find_pair(coloc_results, unique_id("study_a", 1), unique_id("study_a", 2))), 1)

  shared <- find_pair(coloc_results, unique_id("study_a", 1), unique_id("study_b"))
  distinct <- find_pair(coloc_results, unique_id("study_a", 1), unique_id("study_c"))
  expect_gt(shared$h4, 0.8)
  expect_equal(shared$h4, shared$PP.H4.abf)
  expect_lt(distinct$h4, 0.2)
  expect_lt(find_pair(coloc_results, unique_id("study_a", 1), unique_id("study_a", 2))$h4, 0.2)
  expect_true(all(coloc_results$nsnps == 1254))
  expect_true(all(coloc_results$ld_block == fixture_block))
  expect_false(any(coloc_results$ignore))
})

test_that("coloc: pairs with too few shared SNPs are kept but ignored; duplicate and NA LBFs are dropped", {
  paths <- local_block()
  finemapped <- finemapped_fixture_rows(c("study_a", "study_b"))

  b_lbfs <- read_tsv(fixture_path("finemapped", "study_b_1.tsv.gz"))
  few_snps <- place_study("study_few_snps", "finemapped", b_lbfs[1:40, ], "study_few_snps_1.tsv.gz")
  messy <- dplyr::bind_rows(b_lbfs, b_lbfs[1:10, ])
  messy$LBF[500:504] <- NA
  messy_file <- place_study("study_messy", "finemapped", messy, "study_messy_1.tsv.gz")

  b_row <- dplyr::filter(finemapped, study == "study_b")
  extra_rows <- dplyr::bind_rows(
    dplyr::mutate(b_row, study = "study_few_snps", unique_study_id = unique_id("study_few_snps"), file = few_snps),
    dplyr::mutate(b_row, study = "study_messy", unique_study_id = unique_id("study_messy"), file = messy_file)
  )
  first_signal_of_a <- dplyr::filter(finemapped, study == "study_a" & grepl("_1$", unique_study_id))
  write_index(dplyr::bind_rows(first_signal_of_a, extra_rows), "finemapped")

  result <- run_coloc(paths)
  expect_step_succeeded(result, paths$coloc_complete)

  coloc_results <- read_coloc(paths)
  too_few <- find_pair(coloc_results, unique_id("study_a"), unique_id("study_few_snps"))
  expect_equal(nrow(too_few), 1)
  expect_true(too_few$ignore)
  expect_true(is.na(too_few$PP.H4.abf) && is.na(too_few$h4))

  messy_pair <- find_pair(coloc_results, unique_id("study_a"), unique_id("study_messy"))
  expect_equal(messy_pair$nsnps, 1254 - 5)
  expect_gt(messy_pair$h4, 0.8)
})

test_that("coloc: existing pairs are kept and only new pairs are colocalised", {
  paths <- local_block()
  write_index(finemapped_fixture_rows(), "finemapped")
  existing <- data.frame(
    unique_study_a = unique_id("study_a"), study_a = "study_a",
    unique_study_b = unique_id("study_b"), study_b = "study_b",
    bp_distance = 0, ignore = FALSE, false_positive = FALSE, false_negative = FALSE,
    nsnps = 1254, h4 = 0.123, PP.H4.abf = 0.123, ld_block = fixture_block
  )
  vroom::vroom_write(existing, paths$coloc_pairwise)

  result <- run_coloc(paths)
  expect_step_succeeded(result, paths$coloc_complete)

  coloc_results <- read_coloc(paths)
  expect_equal(nrow(coloc_results), 4)
  expect_equal(find_pair(coloc_results, unique_id("study_a"), unique_id("study_b"))$h4, 0.123)
  expect_true(any(grepl("Found 3 study pairs to coloc", result$output)))
})

test_that("coloc: only significant studies that are not ignored are compared", {
  for (left_out in c("not_significant", "ignored")) {
    paths <- local_block()
    finemapped <- finemapped_fixture_rows()
    if (left_out == "not_significant") {
      finemapped$min_p[finemapped$study == "study_c"] <- 0.01
    } else {
      finemapped$ignore[finemapped$study == "study_c"] <- TRUE
    }
    write_index(finemapped, "finemapped")

    result <- run_coloc(paths)
    expect_step_succeeded(result, paths$coloc_complete)
    coloc_results <- read_coloc(paths)
    expect_equal(nrow(coloc_results), 2, info = left_out)
    expect_false(any(coloc_results$study_a == "study_c" | coloc_results$study_b == "study_c"), info = left_out)
  }
})

test_that("coloc: --block_list excludes matching studies", {
  paths <- local_block()
  write_index(finemapped_fixture_rows(), "finemapped")
  block_list_file <- file.path(pipeline_metadata_dir, "unit_block_list.csv")
  vroom::vroom_write(data.frame(id_pattern = "study_c*", cis_trans = NA), block_list_file, delim = ",")

  result <- run_coloc(paths, "--block_list", block_list_file)
  expect_step_succeeded(result, paths$coloc_complete)
  coloc_results <- read_coloc(paths)
  expect_equal(nrow(coloc_results), 2)
  expect_false(any(coloc_results$study_a == "study_c" | coloc_results$study_b == "study_c"))
})

test_that("coloc: blocks that were not updated, or have nothing finemapped, are skipped", {
  paths <- local_block()
  write_index(finemapped_fixture_rows(), "finemapped")
  vroom::vroom_write(
    read_tsv(pipeline_metadata_file_paths()$updated_ld_blocks)[0, ],
    pipeline_metadata_file_paths()$updated_ld_blocks
  )
  result <- run_coloc(paths)
  expect_step_succeeded(result, paths$coloc_complete)
  expect_true(any(grepl("Nothing to coloc", result$output)))
  expect_false(file.exists(paths$coloc_pairwise))

  paths <- local_block()
  result <- run_coloc(paths)
  expect_step_succeeded(result, paths$coloc_complete)
  expect_false(file.exists(paths$coloc_pairwise))
})
