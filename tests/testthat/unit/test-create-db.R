withr::with_dir(pipeline_steps_dir, source("database_definitions.R"))

db_names <- c(
  "studies", "associations_full", "associations_specific", "coloc_pairs_full",
  "coloc_pairs_significant", "ld", "gwas_upload"
)
db_file_names <- c(
  studies = "studies.db", associations_full = "associations_full.db", associations_specific = "associations.db",
  coloc_pairs_full = "coloc_pairs_full.db", coloc_pairs_significant = "coloc_pairs.db", ld = "ld.db",
  gwas_upload = "gwas_upload.db"
)
result_names <- c(
  "studies_processed", "traits_processed", "study_extractions",
  "coloc_clustered_results", "coloc_pairwise_results", "rare_results"
)
latest_studies_db <- file.path(latest_results_dir, "studies.db")

# The compiled results of a real e2e run, and the LD block / study files they point to
local_block()
for (sub_dir in c("ld_blocks", "study")) {
  file.copy(
    list.files(fixture_path("create_db/data", sub_dir), full.names = TRUE),
    file.path(data_dir, sub_dir),
    recursive = TRUE
  )
}

read_fixture_results <- function(name) {
  return(read_tsv(
    fixture_path("create_db/results", paste0(name, ".tsv.gz")),
    col_types = vroom::cols(.default = "c")
  ))
}

#' Copies the compiled results into a fresh results dir, applying any edits (a named list of functions)
prepare_results <- function(scenario, edits = list()) {
  scenario_dir <- file.path(results_dir, paste0("create_db_", scenario))
  unlink(scenario_dir, recursive = TRUE)
  dir.create(scenario_dir, recursive = TRUE)
  for (name in result_names) {
    results <- read_fixture_results(name)
    if (!is.null(edits[[name]])) results <- edits[[name]](results)
    vroom::vroom_write(results, file.path(scenario_dir, paste0(name, ".tsv.gz")), na = "NA")
  }
  return(scenario_dir)
}

db_path <- function(scenario_dir, name) {
  return(file.path(scenario_dir, db_file_names[[name]]))
}

run_create_db <- function(scenario_dir) {
  db_args <- unlist(lapply(db_names, function(name) {
    return(c(paste0("--", name, "_db_file"), db_path(scenario_dir, name)))
  }))
  return(run_step(
    "create_db_from_results.R",
    "--results_dir", scenario_dir,
    db_args,
    "--completed_output_file", file.path(scenario_dir, "create_dbs_done")
  ))
}

db_query <- function(db_file, sql) {
  conn <- duckdb::dbConnect(duckdb::duckdb(), db_file, read_only = TRUE)
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  return(DBI::dbGetQuery(conn, sql))
}

unlink(latest_studies_db)
baseline_dir <- prepare_results("baseline")
baseline_result <- run_create_db(baseline_dir)
baseline_studies_db <- db_path(baseline_dir, "studies")

test_that("create_db: all DBs are created from the compiled results", {
  expect_step_succeeded(baseline_result, file.path(baseline_dir, "create_dbs_done"))
  for (name in db_names) {
    expect_true(file.exists(db_path(baseline_dir, name)), info = name)
  }

  tables <- db_query(baseline_studies_db, "SELECT table_name FROM information_schema.tables")$table_name
  expected_tables <- c(
    vapply(studies_db, function(table) table$name, character(1)),
    vapply(additional_studies_tables, function(table) table$name, character(1))
  )
  expect_true(all(expected_tables %in% tables))

  for (table in gwas_upload_db) {
    count <- db_query(db_path(baseline_dir, "gwas_upload"), glue::glue("SELECT count(*) AS n FROM {table$name}"))$n
    expect_equal(count, 0, info = table$name)
  }
})

test_that("create_db: studies without any extractions are left out", {
  studies <- db_query(baseline_studies_db, "SELECT * FROM studies")
  expect_equal(nrow(studies), 5)
  expect_false("test-height" %in% studies$study_name)
  traits <- db_query(baseline_studies_db, "SELECT * FROM traits")
  expect_equal(sort(traits$id), sort(studies$trait_id))
})

test_that("create_db: every foreign key in the studies DB resolves", {
  ids <- function(table) db_query(baseline_studies_db, glue::glue("SELECT id FROM {table}"))$id
  study_ids <- ids("studies")
  variant_ids <- ids("variant_annotations")
  ld_block_ids <- ids("ld_blocks")
  extraction_ids <- ids("study_extractions")

  study_extractions <- db_query(baseline_studies_db, "SELECT * FROM study_extractions")
  expect_equal(nrow(study_extractions), 6)
  expect_true(all(study_extractions$study_id %in% study_ids))
  expect_true(all(study_extractions$variant_id %in% variant_ids))
  expect_true(all(study_extractions$ld_block_id %in% ld_block_ids))
  # svg and LBF file paths are stored relative to DATA_DIR
  expect_false(any(startsWith(na.omit(study_extractions$svg_file), "/")))
  expect_false(any(startsWith(na.omit(study_extractions$file_with_lbfs), "/")))

  coloc_groups <- db_query(baseline_studies_db, "SELECT * FROM coloc_groups")
  expect_equal(nrow(coloc_groups), 3)
  expect_equal(length(unique(coloc_groups$coloc_group_id)), 1)
  expect_true(all(coloc_groups$study_extraction_id %in% extraction_ids))
  expect_true(all(coloc_groups$variant_id %in% variant_ids))

  rare_results <- db_query(baseline_studies_db, "SELECT * FROM rare_results")
  expect_equal(nrow(rare_results), 2)
  expect_true(all(rare_results$study_extraction_id %in% extraction_ids))
})

test_that("create_db: coloc pair DBs only hold usable pairs", {
  coloc_pairs_full <- db_query(db_path(baseline_dir, "coloc_pairs_full"), "SELECT * FROM coloc_pairs")
  expect_equal(nrow(coloc_pairs_full), 3)
  expect_false(any(is.na(coloc_pairs_full$h4)))
  coloc_pairs <- db_query(db_path(baseline_dir, "coloc_pairs_significant"), "SELECT * FROM coloc_pairs")
  expect_true(all(coloc_pairs$h4 > posterior_prob_threshold_minimum))
})

test_that("create_db: associations have no missing values", {
  full_tables <- db_query(db_path(baseline_dir, "associations_full"), "SELECT * FROM associations_metadata")
  expect_gte(nrow(full_tables), 1)
  general <- db_query(db_path(baseline_dir, "associations_full"), glue::glue(
    "SELECT * FROM {full_tables$associations_table_name[1]}"
  ))
  specific <- db_query(db_path(baseline_dir, "associations_specific"), "SELECT * FROM associations")
  expect_gt(nrow(general), 0)
  expect_gt(nrow(specific), 0)
  for (associations in list(general, specific)) {
    expect_false(any(is.na(associations[, c("beta", "se", "p", "eaf")])))
  }
})

test_that("create_db: LD correlations are sign-flipped to match flipped variants", {
  ld <- db_query(db_path(baseline_dir, "ld"), "SELECT * FROM ld")
  expect_gt(nrow(ld), 0)
  expect_true(all(abs(ld$r) <= 1))

  variants <- db_query(baseline_studies_db, "SELECT id, snp, flipped FROM variant_annotations")
  ld_blocks <- db_query(baseline_studies_db, "SELECT id, ld_block FROM ld_blocks")
  block <- "EUR/1/1170341-1730405"
  ld <- dplyr::filter(ld, ld_block_id == ld_blocks$id[ld_blocks$ld_block == block])
  ld$lead <- variants$snp[match(ld$lead_variant_id, variants$id)]
  ld$proxy <- variants$snp[match(ld$proxy_variant_id, variants$id)]
  lead_flipped <- variants$flipped[match(ld$lead_variant_id, variants$id)]
  ld$flip <- lead_flipped != variants$flipped[match(ld$proxy_variant_id, variants$id)]

  block_paths <- ld_block_file_paths(block)
  vars <- read_tsv(block_paths$ld_matrix_vars, col_names = "snp", delim = " ")$snp
  ld_matrix <- as.matrix(read_tsv(block_paths$ld_matrix_vcor, col_names = FALSE, altrep = FALSE))
  raw_r <- ld_matrix[cbind(match(ld$lead, vars), match(ld$proxy, vars))]
  expect_equal(ld$r, ifelse(ld$flip, -raw_r, raw_r), tolerance = 1e-6)
})

test_that("create_db: opengwas source urls point at the dataset page", {
  source_urls <- db_query(baseline_studies_db, paste(
    "SELECT DISTINCT studies.study_name, coloc_groups_wide.source_url FROM coloc_groups_wide",
    "JOIN studies ON coloc_groups_wide.study_id = studies.id"
  ))
  expect_equal(
    source_urls$source_url[source_urls$study_name == "ebi-a-GCST90028992"],
    "https://opengwas.io/datasets/ebi-a-GCST90028992"
  )
})

test_that("create_db: without KEGG or STRING files the pathway tables are empty", {
  expect_equal(db_query(baseline_studies_db, "SELECT count(*) AS n FROM pathway_mappings")$n, 0)
  expect_equal(db_query(baseline_studies_db, "SELECT count(*) AS n FROM pathway_sizes")$n, 0)
})

test_that("create_db: ids are kept from the latest studies DB, and new rows get the next id", {
  expect_step_succeeded(baseline_result)
  # a latest studies DB with only the persisted id columns, where studies and study_extractions have different ids
  # to the ones the baseline run gave them, and one study is missing
  dir.create(latest_results_dir, recursive = TRUE, showWarnings = FALSE)
  unlink(latest_studies_db)
  withr::defer(unlink(latest_studies_db))
  conn <- duckdb::dbConnect(duckdb::duckdb(), latest_studies_db)
  for (table in studies_db) {
    if (is.na(table$persist_id_from)) next
    ids <- db_query(baseline_studies_db, glue::glue("SELECT id, {table$persist_id_from} FROM {table$name}"))
    if (table$name %in% c("studies", "study_extractions")) ids$id <- ids$id + 100
    if (table$name == "studies") ids <- dplyr::filter(ids, study_name != "ukb-wes-bm-00000297")
    DBI::dbWriteTable(conn, table$name, ids)
  }
  DBI::dbDisconnect(conn, shutdown = TRUE)

  latest_studies <- db_query(latest_studies_db, "SELECT id, study_name FROM studies")
  latest_extractions <- db_query(latest_studies_db, "SELECT id, unique_study_id FROM study_extractions")

  persisted_dir <- prepare_results("persisted_ids")
  result <- run_create_db(persisted_dir)
  expect_step_succeeded(result, file.path(persisted_dir, "create_dbs_done"))

  studies <- db_query(db_path(persisted_dir, "studies"), "SELECT id, study_name FROM studies")
  kept <- dplyr::inner_join(studies, latest_studies, by = "study_name", suffix = c("", "_latest"))
  expect_equal(nrow(kept), 4)
  expect_equal(kept$id, kept$id_latest)
  expect_equal(studies$id[studies$study_name == "ukb-wes-bm-00000297"], max(latest_studies$id) + 1)

  extractions <- db_query(db_path(persisted_dir, "studies"), "SELECT id, unique_study_id FROM study_extractions")
  kept_extractions <- dplyr::inner_join(
    extractions, latest_extractions,
    by = "unique_study_id", suffix = c("", "_latest")
  )
  expect_equal(nrow(kept_extractions), 6)
  expect_equal(kept_extractions$id, kept_extractions$id_latest)
})

test_that("create_db: rows that cannot be linked are dropped and logged", {
  unknown_extraction <- "unknown-study_EUR/1/1170341-1730405_1"
  edits <- list(
    study_extractions = function(results) {
      missing_snp <- results$study == "ebi-a-GCST90028992" & results$ld_block == "EUR/1/1730405-3355587"
      results$snp[missing_snp] <- "1:2539237_A_T"
      ignored <- results[missing_snp, ] |>
        dplyr::mutate(unique_study_id = sub("_1$", "_2", unique_study_id), snp = "1:2539236_C_T", ignore = "TRUE")
      return(dplyr::bind_rows(results, ignored))
    },
    coloc_clustered_results = function(results) {
      return(dplyr::bind_rows(results, dplyr::mutate(results[1, ], unique_study_id = unknown_extraction)))
    },
    coloc_pairwise_results = function(results) {
      return(dplyr::bind_rows(results, dplyr::mutate(results[1, ], unique_study_b = unknown_extraction)))
    },
    rare_results = function(results) {
      return(dplyr::bind_rows(results, dplyr::mutate(results[1, ], traits = "unknown-rare_EUR_1_2591624-C-T")))
    }
  )
  scenario_dir <- prepare_results("unlinked_rows", edits)
  result <- run_create_db(scenario_dir)
  expect_step_succeeded(result, file.path(scenario_dir, "create_dbs_done"))
  studies_db_file <- db_path(scenario_dir, "studies")

  study_extractions <- db_query(studies_db_file, "SELECT unique_study_id FROM study_extractions")$unique_study_id
  expect_equal(length(study_extractions), 5)
  expect_false("ebi-a-GCST90028992_EUR/1/1730405-3355587_1" %in% study_extractions)
  expect_false("ebi-a-GCST90028992_EUR/1/1730405-3355587_2" %in% study_extractions)
  expect_true(file.exists(file.path(scenario_dir, "missing_snps_in_study_extractions.tsv")))

  expect_equal(db_query(studies_db_file, "SELECT count(*) AS n FROM coloc_groups")$n, 3)
  expect_true(file.exists(file.path(scenario_dir, "missing_stuff_in_clustered_colocs.tsv")))

  coloc_pairs_full <- db_query(db_path(scenario_dir, "coloc_pairs_full"), "SELECT * FROM coloc_pairs")
  expect_equal(nrow(coloc_pairs_full), 3)
  expect_true(file.exists(file.path(scenario_dir, "missing_study_extractions_in_pairwise_colocs.tsv")))

  expect_equal(db_query(studies_db_file, "SELECT count(*) AS n FROM rare_results")$n, 2)
  expect_true(file.exists(file.path(scenario_dir, "missing_study_extractions_in_rare_results.tsv")))
})

test_that("create_db: genes are resolved by ensg first, then by name or unambiguous alias", {
  tnfrsf4 <- "ENSG00000186827"
  tnfrsf18 <- "ENSG00000186891"
  edits <- list(studies_processed = function(results) {
    eqtl <- results$study_name == "qtl-GTEx-eQTL-v10-Brain-Cortex-ENSG00000186827-11"
    sqtl <- grepl("GTEx-sQTL", results$study_name)
    phenotype <- results$study_name == "ebi-a-GCST90028992"
    methylation <- results$study_name == "ukb-wes-bm-00000040"
    # ensg wins over a gene name pointing at a different gene
    results$gene[eqtl] <- "TNFRSF18"
    # no ensg: resolved through the HGNC alias GITR
    results$ensg[sqtl] <- NA
    results$gene[sqtl] <- "GITR"
    # FLIP is an alias of 2 genes, so it is left unresolved
    results$gene[phenotype] <- "FLIP"
    # methylation studies never get a gene
    results$data_type[methylation] <- data_types$methylation
    results$gene[methylation] <- "TNFRSF4"
    results$ensg[methylation] <- tnfrsf4
    return(results)
  })
  scenario_dir <- prepare_results("genes", edits)
  result <- run_create_db(scenario_dir)
  expect_step_succeeded(result, file.path(scenario_dir, "create_dbs_done"))

  studies <- db_query(db_path(scenario_dir, "studies"), paste(
    "SELECT studies.study_name, studies.gene, gene_annotations.ensembl_id FROM studies",
    "LEFT JOIN gene_annotations ON studies.gene_id = gene_annotations.id"
  ))
  gene_of <- function(pattern) studies$ensembl_id[grepl(pattern, studies$study_name)]
  expect_equal(gene_of("GTEx-eQTL"), tnfrsf4)
  expect_equal(gene_of("GTEx-sQTL"), tnfrsf18)
  expect_true(is.na(gene_of("ebi-a-GCST90028992")))
  expect_true(is.na(gene_of("ukb-wes-bm-00000040")))
  expect_true(is.na(studies$gene[studies$study_name == "ukb-wes-bm-00000040"]))
})

test_that("create_db: duplicate studies or study extractions fail the step", {
  duplicate_extraction_dir <- prepare_results("duplicate_extractions", list(
    study_extractions = function(results) dplyr::bind_rows(results, results[1, ])
  ))
  result <- run_create_db(duplicate_extraction_dir)
  expect_step_failed(result, file.path(duplicate_extraction_dir, "create_dbs_done"), "duplicate unique_study_id")
  expect_true(file.exists(file.path(duplicate_extraction_dir, "duplicate_study_extractions.tsv")))

  duplicate_study_dir <- prepare_results("duplicate_studies", list(
    studies_processed = function(results) dplyr::bind_rows(results, results[1, ])
  ))
  result <- run_create_db(duplicate_study_dir)
  expect_step_failed(result, file.path(duplicate_study_dir, "create_dbs_done"), "duplicate study_name")
  expect_true(file.exists(file.path(duplicate_study_dir, "duplicate_studies.tsv")))
})
