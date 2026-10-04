# Shared setup for the per-step tests. Every pipeline step is run as `Rscript <step>.R ...` from pipeline_steps/,
# the same way the Snakefile runs it, against the permanent fixtures built by tests/generate_unit_fixtures.R.
library(testthat)

repo_root <- normalizePath("../../..")
pipeline_steps_dir <- file.path(repo_root, "pipeline_steps")
unit_root_dir <- "/local-scratch/projects/genotype-phenotype-map/test/unit"
unit_fixtures_dir <- file.path(unit_root_dir, "fixtures")

Sys.setenv("TEST_RUN" = "test")
# vroom defaults to one thread per core, which on a many-core box makes every small read slow (~15s per step)
Sys.setenv("VROOM_THREADS" = 4)
Sys.setenv("DATA_DIR" = paste0(unit_root_dir, "/data/"))
Sys.setenv("RESULTS_DIR" = paste0(unit_root_dir, "/results/"))
withr::with_dir(pipeline_steps_dir, source("constants.R"))

fixture_block <- "EUR/1/16103-1170341"
coloc_block <- "EUR/1/1170341-1730405"
simulation_info <- vroom::vroom(file.path(unit_fixtures_dir, "gwas/simulation_info.tsv"), show_col_types = FALSE)
fixture_sample_size <- simulation_info$sample_size
causal_snps <- vroom::vroom(file.path(unit_fixtures_dir, "gwas/causal_snps.tsv"), show_col_types = FALSE)

fixture_path <- function(...) {
  return(file.path(unit_fixtures_dir, ...))
}

read_fixture_gwas <- function(name) {
  return(vroom::vroom(fixture_path("gwas", paste0(name, ".tsv.gz")), show_col_types = FALSE))
}

read_tsv <- function(file, ...) {
  return(vroom::vroom(file, show_col_types = FALSE, ...))
}

#' Runs a pipeline step script with Rscript, from pipeline_steps/, as the Snakefile does
run_step <- function(script, ...) {
  step_args <- as.character(unlist(list(...)))
  output <- withr::with_dir(
    pipeline_steps_dir,
    suppressWarnings(system2("Rscript", shQuote(c(script, step_args)), stdout = TRUE, stderr = TRUE))
  )
  status <- attr(output, "status")
  if (is.null(status)) status <- 0
  return(list(status = status, output = output))
}

#' The Snakefile ignores exit codes, so a step has only succeeded if it also wrote its sentinel file
expect_step_succeeded <- function(result, completed_output_file = NULL) {
  step_output <- paste(utils::tail(result$output, 40), collapse = "\n")
  expect_equal(result$status, 0, info = step_output)
  if (!is.null(completed_output_file)) {
    expect_true(file.exists(completed_output_file), info = step_output)
  }
  return(invisible(result))
}

#' constants.R sets options(error = traceback), so a failing step still exits with status 0: the error message
#' and the missing sentinel file are the only signs of failure
expect_step_failed <- function(result, completed_output_file = NULL, error_pattern = "Error") {
  expect_true(
    any(grepl(error_pattern, result$output)),
    info = paste(utils::tail(result$output, 40), collapse = "\n")
  )
  if (!is.null(completed_output_file)) {
    expect_false(file.exists(completed_output_file))
  }
  return(invisible(result))
}

#' Wipes the unit data dir at the start of each test (so the last test's state is left for debugging),
#' then creates a fresh LD block dir and marks it as updated, like organise_extracted_regions_into_ld_blocks.R
local_block <- function(ld_block = fixture_block) {
  for (dir in c(ld_block_data_dir, extracted_study_dir, pipeline_metadata_dir)) {
    unlink(list.files(dir, full.names = TRUE), recursive = TRUE)
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  }
  paths <- ld_block_file_paths(ld_block)
  dir.create(paths$ld_block_data, recursive = TRUE, showWarnings = FALSE)

  components <- ld_block_components(ld_block)
  updated_ld_blocks <- dplyr::mutate(
    components,
    ld_block = ld_block,
    data_dir = as.character(paths$ld_block_data)
  )
  vroom::vroom_write(updated_ld_blocks, pipeline_metadata_file_paths()$updated_ld_blocks)
  return(paths)
}

#' Writes a GWAS into data/study/<study>/<form>/ (creating the dirs the steps expect) and returns its path
place_study <- function(study, form, gwas, file_name = "EUR_1_fixture.tsv.gz") {
  study_dir <- file.path(extracted_study_dir, study)
  for (sub_dir in c("extracted", "standardised", "imputed", "finemapped", "svgs/extractions")) {
    dir.create(file.path(study_dir, sub_dir), recursive = TRUE, showWarnings = FALSE)
  }
  file <- file.path(study_dir, form, file_name)
  vroom::vroom_write(gwas, file)
  return(file)
}

#' One row of an extracted / standardised / imputed studies index
study_index_row <- function(study,
                            file,
                            ld_block = fixture_block,
                            bp = 1000000,
                            coverage = coverage_types$dense,
                            variant_type = variant_types$common,
                            cis_trans = NA,
                            sample_size = fixture_sample_size,
                            category = study_categories$continuous) {
  return(data.frame(
    study = study,
    file = file,
    ancestry = "EUR",
    chr = "1",
    bp = bp,
    p_value_threshold = lowest_p_value_threshold,
    category = category,
    sample_size = sample_size,
    cis_trans = cis_trans,
    reference_build = reference_builds$GRCh38,
    ld_block = ld_block,
    variant_type = variant_type,
    coverage = coverage
  ))
}

write_index <- function(rows, form, ld_block = fixture_block) {
  index_file <- ld_block_file_paths(ld_block)[[paste0(form, "_studies")]]
  vroom::vroom_write(rows, index_file)
  return(index_file)
}

#' Copies the simulated finemapped (LBF) fixtures into data/study and returns a full finemapped_studies index
finemapped_fixture_rows <- function(studies = c("study_a", "study_b", "study_c"), ld_block = fixture_block) {
  template <- read_tsv(fixture_path("finemapped", "finemapped_studies_template.tsv")) |>
    dplyr::filter(study %in% studies)
  files <- vapply(seq_len(nrow(template)), function(i) {
    gwas <- read_tsv(fixture_path("finemapped", template$file_name[i]))
    return(place_study(template$study[i], "finemapped", gwas, template$file_name[i]))
  }, character(1))

  return(data.frame(
    study = template$study,
    unique_study_id = paste0(template$study, "_", ld_block, "_", template$cs),
    ld_block = ld_block,
    variant_type = variant_types$common,
    file = files,
    ancestry = "EUR",
    chr = "1",
    bp = template$bp,
    snp = template$snp,
    p_value_threshold = lowest_p_value_threshold,
    min_p = template$min_p,
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
  ))
}

#' One row of pipeline_metadata/studies_to_process.tsv, as written by identify_studies_to_process.R
studies_to_process_row <- function(study_name,
                                   data_format,
                                   study_location,
                                   reference_build = reference_builds$GRCh38,
                                   p_value_threshold = lowest_p_value_threshold,
                                   variant_type = variant_types$common,
                                   data_type = data_types$phenotype,
                                   probe = NA,
                                   file_type = NA,
                                   column_names = NA) {
  return(data.frame(
    data_type = data_type,
    data_format = data_format,
    source = "unit_test",
    study_name = study_name,
    trait = study_name,
    trait_name = study_name,
    ancestry = "EUR",
    sample_size = fixture_sample_size,
    category = study_categories$continuous,
    study_location = study_location,
    extracted_location = paste0(extracted_study_dir, study_name, "/"),
    reference_build = reference_build,
    p_value_threshold = p_value_threshold,
    variant_type = variant_type,
    gene = NA,
    probe = probe,
    tissue = NA,
    cell_type = NA,
    coverage = coverage_types$dense,
    heritability = NA,
    heritability_se = NA,
    ensg = NA,
    file_type = file_type,
    column_names = column_names
  ))
}

write_studies_to_process <- function(rows) {
  vroom::vroom_write(rows, pipeline_metadata_file_paths()$studies_to_process)
  return(invisible(rows))
}
