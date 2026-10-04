cluster_step <- "cluster_studies_in_ld_block.R"

run_cluster <- function(paths, ...) {
  return(run_step(
    cluster_step,
    "--ld_block", fixture_block,
    "--completed_output_file", paths$clustering_complete,
    ...
  ))
}

cluster_snps <- c("1:1000_A_G", "1:2000_A_G", "1:3000_A_G")
clustered_columns <- c("unique_study_id", "component", "ld_block", "snp", "h4_connectedness", "h3_connectedness")

#' Writes a synthetic finemapped_studies.tsv (with small LBF files) and coloc_pairwise_results.tsv.gz.
#' studies: unique_study_id, study, min_p, bp, and an lbf list column over cluster_snps
#' pairs: a, b, h4 (and optionally ignore)
write_cluster_inputs <- function(studies, pairs) {
  files <- vapply(seq_len(nrow(studies)), function(i) {
    lbf_gwas <- data.frame(SNP = cluster_snps, LBF = studies$lbf[[i]])
    return(place_study(studies$study[i], "finemapped", lbf_gwas, paste0(studies$unique_study_id[i], ".tsv.gz")))
  }, character(1))

  finemapped <- data.frame(
    study = studies$study,
    unique_study_id = studies$unique_study_id,
    ld_block = fixture_block,
    variant_type = variant_types$common,
    file = files,
    ancestry = "EUR",
    chr = "1",
    bp = studies$bp,
    snp = cluster_snps[1],
    p_value_threshold = lowest_p_value_threshold,
    min_p = studies$min_p,
    category = study_categories$continuous,
    sample_size = fixture_sample_size,
    cis_trans = NA,
    finemap_message = "success",
    first_finemap_num_results = NA,
    second_finemap_num_results = NA,
    qc_step_run = FALSE,
    snps_removed_by_qc = NA,
    time_taken = "00:00:01",
    svg_file = NA,
    file_with_lbfs = NA,
    ignore = FALSE,
    coverage = coverage_types$dense
  )
  write_index(finemapped, "finemapped")

  if (!"ignore" %in% colnames(pairs)) pairs$ignore <- FALSE
  study_of <- setNames(studies$study, studies$unique_study_id)
  coloc_results <- data.frame(
    unique_study_a = pairs$a,
    study_a = unname(study_of[pairs$a]),
    unique_study_b = pairs$b,
    study_b = unname(study_of[pairs$b]),
    bp_distance = 0,
    ignore = pairs$ignore,
    false_positive = FALSE,
    false_negative = FALSE,
    nsnps = 1000,
    PP.H0.abf = 0,
    PP.H1.abf = 0,
    PP.H2.abf = 0,
    PP.H3.abf = ifelse(is.na(pairs$h4), NA, 1 - pairs$h4),
    PP.H4.abf = pairs$h4,
    h4 = pairs$h4,
    ld_block = fixture_block
  )
  vroom::vroom_write(coloc_results, ld_block_file_paths(fixture_block)$coloc_pairwise)
  return(invisible(NULL))
}

make_studies <- function(ids, studies = ids, min_p = 1e-10, bp = 1000, lbf = NULL) {
  if (is.null(lbf)) lbf <- rep(list(c(5, 1, 0)), length(ids))
  result <- data.frame(unique_study_id = ids, study = studies, min_p = min_p, bp = bp)
  result$lbf <- lbf
  return(result)
}

clique_pairs <- function(ids, h4 = 0.95) {
  combinations <- utils::combn(ids, 2)
  return(data.frame(a = combinations[1, ], b = combinations[2, ], h4 = h4))
}

read_clusters <- function(paths) {
  return(read_tsv(paths$clustered) |> dplyr::arrange(unique_study_id))
}

component_of <- function(clusters, id) {
  return(clusters$component[clusters$unique_study_id == id])
}

test_that("cluster: two separate groups of colocalising studies become two components", {
  paths <- local_block()
  ids <- paste0("s", 1:7)
  studies <- make_studies(ids, lbf = list(
    c(5, 4, 0), c(0, 4, 1), c(1, 3, 0),
    c(0, 1, 9), c(0, 2, 8), c(1, 0, 7),
    c(9, 0, 0)
  ))
  pairs <- dplyr::bind_rows(
    clique_pairs(c("s1", "s2", "s3")),
    clique_pairs(c("s4", "s5", "s6")),
    data.frame(a = c("s1", "s3", "s7"), b = c("s4", "s6", "s1"), h4 = c(0.1, 0.2, 0.1))
  )
  pairs$h4[pairs$a == "s1" & pairs$b == "s3"] <- 0.5
  write_cluster_inputs(studies, pairs)

  result <- run_cluster(paths)
  expect_step_succeeded(result, paths$clustering_complete)

  clusters <- read_clusters(paths)
  expect_equal(colnames(clusters), clustered_columns)
  # s7 only has a weak coloc with s1, so is dropped as a singleton
  expect_equal(clusters$unique_study_id, paste0("s", 1:6))
  expect_equal(length(unique(clusters$component)), 2)
  expect_equal(length(unique(clusters$component[clusters$unique_study_id %in% c("s1", "s2", "s3")])), 1)
  expect_false(component_of(clusters, "s1") == component_of(clusters, "s4"))

  # the cluster SNP is the one with the highest summed LBF, not each study's own top SNP
  expect_equal(unique(clusters$snp[clusters$unique_study_id %in% c("s1", "s2", "s3")]), cluster_snps[2])
  expect_equal(unique(clusters$snp[clusters$unique_study_id %in% c("s4", "s5", "s6")]), cluster_snps[3])
  expect_equal(unique(clusters$h4_connectedness[clusters$unique_study_id == "s1"]), 2 / 3)
  expect_equal(unique(clusters$h4_connectedness[clusters$unique_study_id == "s4"]), 1)

  # s1 and s3 end up in the same cluster without colocalising with each other
  coloc_results <- read_tsv(paths$coloc_pairwise)
  expect_true(coloc_results$false_negative[coloc_results$unique_study_a == "s1" & coloc_results$unique_study_b == "s3"])
  expect_equal(sum(coloc_results$false_negative), 1)
  expect_true(file.exists(paths$igraph_clustered))
})

test_that("cluster: a single group of colocalising studies becomes one component", {
  paths <- local_block()
  write_cluster_inputs(make_studies(c("s1", "s2", "s3")), clique_pairs(c("s1", "s2", "s3")))

  result <- run_cluster(paths)
  expect_step_succeeded(result, paths$clustering_complete)

  clusters <- read_clusters(paths)
  expect_equal(clusters$unique_study_id, c("s1", "s2", "s3"))
  expect_true(all(clusters$component == 1))
})

test_that("cluster: an edge bridging two communities is pruned and marked as a false positive", {
  paths <- local_block()
  ids <- paste0("s", 1:8)
  pairs <- dplyr::bind_rows(
    clique_pairs(c("s1", "s2", "s3")),
    clique_pairs(c("s4", "s5", "s6")),
    clique_pairs(c("s7", "s8")),
    data.frame(a = "s3", b = "s4", h4 = 0.95)
  )
  write_cluster_inputs(make_studies(ids), pairs)

  result <- run_cluster(paths)
  expect_step_succeeded(result, paths$clustering_complete)

  clusters <- read_clusters(paths)
  expect_false(component_of(clusters, "s3") == component_of(clusters, "s4"))
  coloc_results <- read_tsv(paths$coloc_pairwise)
  expect_true(coloc_results$false_positive[coloc_results$unique_study_a == "s3" & coloc_results$unique_study_b == "s4"])
  expect_equal(sum(coloc_results$false_positive), 1)
})

test_that("cluster: of two credible sets from one study in a cluster, the most significant is kept", {
  paths <- local_block()
  ids <- c("study_x_1", "study_x_2", "s2", "s3", "s4", "s5")
  studies <- make_studies(
    ids,
    studies = c("study_x", "study_x", "s2", "s3", "s4", "s5"),
    min_p = c(1e-10, 1e-20, 1e-10, 1e-10, 1e-10, 1e-10)
  )
  pairs <- dplyr::bind_rows(
    clique_pairs(c("study_x_1", "s2", "s3")),
    clique_pairs(c("study_x_2", "s2", "s3")),
    clique_pairs(c("s4", "s5")),
    data.frame(a = "study_x_1", b = "study_x_2", h4 = 0.1)
  ) |> dplyr::distinct(a, b, .keep_all = TRUE)
  write_cluster_inputs(studies, pairs)

  result <- run_cluster(paths)
  expect_step_succeeded(result, paths$clustering_complete)
  clusters <- read_clusters(paths)
  expect_true("study_x_2" %in% clusters$unique_study_id)
  expect_false("study_x_1" %in% clusters$unique_study_id)
})

test_that("cluster: when credible sets from one study tie on min_p (both 0), the lowest bp is kept", {
  paths <- local_block()
  ids <- c("study_x_1", "study_x_2", "s2", "s3", "s4", "s5")
  studies <- make_studies(
    ids,
    studies = c("study_x", "study_x", "s2", "s3", "s4", "s5"),
    min_p = c(0, 0, 1e-10, 1e-10, 1e-10, 1e-10),
    bp = c(2000, 1000, 1000, 1000, 1000, 1000)
  )
  pairs <- dplyr::bind_rows(
    clique_pairs(c("study_x_1", "s2", "s3")),
    clique_pairs(c("study_x_2", "s2", "s3")),
    clique_pairs(c("s4", "s5")),
    data.frame(a = "study_x_1", b = "study_x_2", h4 = 0.1)
  ) |> dplyr::distinct(a, b, .keep_all = TRUE)
  write_cluster_inputs(studies, pairs)

  result <- run_cluster(paths)
  expect_step_succeeded(result, paths$clustering_complete)
  clusters <- read_clusters(paths)
  expect_true("study_x_2" %in% clusters$unique_study_id)
  expect_false("study_x_1" %in% clusters$unique_study_id)
})

test_that("cluster: duplicate credible sets from one study are pruned when there is only one component", {
  paths <- local_block()
  ids <- c("study_x_1", "study_x_2", "s2", "s3")
  studies <- make_studies(ids, studies = c("study_x", "study_x", "s2", "s3"), min_p = c(1e-10, 1e-20, 1e-10, 1e-10))
  pairs <- dplyr::bind_rows(
    clique_pairs(c("study_x_1", "s2", "s3")),
    clique_pairs(c("study_x_2", "s2", "s3")),
    data.frame(a = "study_x_1", b = "study_x_2", h4 = 0.1)
  ) |> dplyr::distinct(a, b, .keep_all = TRUE)
  write_cluster_inputs(studies, pairs)

  result <- run_cluster(paths)
  expect_step_succeeded(result, paths$clustering_complete)
  clusters <- read_clusters(paths)
  expect_equal(clusters$unique_study_id, c("s2", "s3", "study_x_2"))
  expect_true(all(clusters$component == 1))

  # the removed credible set's colocalising edges are marked as false positives
  coloc_results <- read_tsv(paths$coloc_pairwise)
  removed_pairs <- coloc_results$unique_study_a == "study_x_1" | coloc_results$unique_study_b == "study_x_1"
  expect_true(all(coloc_results$false_positive[removed_pairs & coloc_results$h4 > 0.8]))
})

test_that("cluster: credible sets of one study that colocalise with each other drop the weaker one", {
  paths <- local_block()
  ids <- c("study_x_1", "study_x_2", "s2", "s3")
  studies <- make_studies(ids, studies = c("study_x", "study_x", "s2", "s3"), min_p = c(1e-20, 1e-8, 1e-10, 1e-10))
  write_cluster_inputs(studies, clique_pairs(ids))

  result <- run_cluster(paths)
  expect_step_succeeded(result, paths$clustering_complete)

  clusters <- read_clusters(paths)
  expect_equal(clusters$unique_study_id, c("s2", "s3", "study_x_1"))
  coloc_results <- read_tsv(paths$coloc_pairwise)
  same_study_pair <- coloc_results$unique_study_a == "study_x_1" & coloc_results$unique_study_b == "study_x_2"
  expect_true(coloc_results$ignore[same_study_pair])
})

test_that("cluster: weak and missing H4 values do not join studies together", {
  paths <- local_block()
  ids <- paste0("s", 1:4)
  pairs <- data.frame(
    a = c("s1", "s1", "s1"),
    b = c("s2", "s3", "s4"),
    h4 = c(0.95, 0.79, NA),
    ignore = c(FALSE, FALSE, TRUE)
  )
  write_cluster_inputs(make_studies(ids), pairs)

  result <- run_cluster(paths)
  expect_step_succeeded(result, paths$clustering_complete)
  expect_equal(read_clusters(paths)$unique_study_id, c("s1", "s2"))
})

test_that("cluster: --block_list removes pairs with blocked studies", {
  paths <- local_block()
  ids <- c("s1", "s2", "blocked_study")
  write_cluster_inputs(make_studies(ids), clique_pairs(ids))
  block_list_file <- file.path(pipeline_metadata_dir, "unit_block_list.csv")
  vroom::vroom_write(data.frame(id_pattern = "blocked*", cis_trans = NA), block_list_file, delim = ",")

  result <- run_cluster(paths, "--block_list", block_list_file)
  expect_step_succeeded(result, paths$clustering_complete)

  block_list_paths <- ld_block_file_paths(fixture_block, block_list = block_list_file)
  expect_equal(basename(block_list_paths$clustered), "clustered_results_unit.tsv.gz")
  expect_equal(read_clusters(block_list_paths)$unique_study_id, c("s1", "s2"))
})

test_that("cluster: a block with no coloc results writes an empty clustered file", {
  paths <- local_block()
  write_index(finemapped_fixture_rows(), "finemapped")

  result <- run_cluster(paths)
  expect_step_succeeded(result, paths$clustering_complete)
  clusters <- read_tsv(paths$clustered)
  expect_equal(nrow(clusters), 0)
  expect_equal(colnames(clusters), clustered_columns)
})
