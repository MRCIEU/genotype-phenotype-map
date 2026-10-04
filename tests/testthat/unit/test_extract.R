all_ld_blocks <- withr::with_dir(pipeline_steps_dir, read_tsv("data/ld_blocks.tsv")) |>
  dplyr::filter(ancestry == "EUR") |>
  dplyr::mutate(ld_block = ld_block_string(ancestry, chr, start, stop))
block_bounds <- function(ld_block) {
  return(all_ld_blocks[all_ld_blocks$ld_block == ld_block, ])
}
second_block <- "EUR/1/1170341-1730405"
third_block <- "EUR/1/1730405-3355587"

run_extraction <- function(script, study_row, ...) {
  return(run_step(
    script,
    "--extracted_study_location", study_row$extracted_location,
    "--extracted_output_file", extracted_output_file(study_row),
    ...
  ))
}

extracted_output_file <- function(study_row) {
  return(file.path(study_row$extracted_location, "extracted_snps.tsv"))
}

read_extracted_snps <- function(study_row) {
  file <- extracted_output_file(study_row)
  if (file.size(file) == 0) {
    return(data.frame())
  }
  return(read_tsv(file, col_types = vroom::cols(.default = "?", chr = "c", cis_trans = "c")))
}

#' Every extracted region file only holds SNPs from inside its LD block
expect_extractions_inside_blocks <- function(extracted_snps) {
  for (i in seq_len(nrow(extracted_snps))) {
    bounds <- block_bounds(extracted_snps$ld_block[i])
    extraction <- read_tsv(extracted_snps$file[i])
    expect_true(nrow(extraction) > 0)
    expect_true(all(extraction$BP >= bounds$start & extraction$BP <= bounds$stop), info = extracted_snps$file[i])
  }
  return(invisible(NULL))
}

input_dir <- function() {
  dir <- file.path(extracted_study_dir, "unit_inputs")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  return(dir)
}

#' A GWAS over 3 LD blocks: block 1 has a 1e-27 hit, block 2 a 1e-5 hit, block 3 nothing below 1e-3
summary_stats_gwas <- function() {
  block_1 <- read_fixture_gwas("study_a_imputed") |> dplyr::select(CHR, BP, EA, OA, EAF, BETA, SE, P)
  panel_2 <- read_tsv(ld_block_file_paths(second_block)$ld_matrix_tsv)[1:300, ]
  # the panels of neighbouring blocks share their boundary SNP, so only keep it once
  block_2 <- dplyr::transmute(panel_2, CHR, BP, EA, OA, EAF, BETA = 0.001, SE = 0.01, P = 0.5) |>
    dplyr::filter(!BP %in% block_1$BP)
  block_2$P[150] <- 1e-5
  panel_3 <- read_tsv(ld_block_file_paths(third_block)$ld_matrix_tsv)[1:200, ]
  block_3 <- dplyr::transmute(panel_3, CHR, BP, EA, OA, EAF, BETA = 0.001, SE = 0.01, P = 0.5)
  block_3$P[100] <- 1e-3
  return(dplyr::bind_rows(block_1, block_2, block_3))
}

standard_column_names <- '{"CHR":"CHR","BP":"BP","EA":"EA","OA":"OA","P":"P","BETA":"BETA","SE":"SE","EAF":"EAF"}'

summary_stats_study <- function(gwas,
                                study_name = "unit-summary-stats",
                                column_names = standard_column_names,
                                ...) {
  study_location <- file.path(input_dir(), paste0(study_name, ".tsv.gz"))
  vroom::vroom_write(gwas, study_location)
  study_row <- studies_to_process_row(
    study_name, data_formats$summary_stats, study_location,
    column_names = column_names, ...
  )
  return(study_row)
}

run_summary_stats <- function(study_row, ...) {
  write_studies_to_process(study_row)
  return(run_extraction("extract_regions_from_summary_stats.R", study_row, ...))
}

# summary_stats ---------------------------------------------------------------------------------------------------

test_that("extract summary_stats: one region per LD block with a hit, sorted by significance", {
  local_block()
  study_row <- summary_stats_study(summary_stats_gwas())

  result <- run_summary_stats(study_row)
  expect_step_succeeded(result, extracted_output_file(study_row))

  extracted_snps <- read_extracted_snps(study_row)
  expect_equal(extracted_snps$ld_block, c(fixture_block, second_block))
  expect_equal(extracted_snps$bp[1], causal_snps$bp[causal_snps$name == "causal_1"])
  expect_true(all(diff(extracted_snps$log_p) <= 0))
  expect_true(all(is.na(extracted_snps$cis_trans)))
  expect_extractions_inside_blocks(extracted_snps)
  expect_true(file.exists(file.path(study_row$extracted_location, "svgs/full.zip")))
  expect_true(file.exists(file.path(study_row$extracted_location, "svgs/full.json")))
  for (sub_dir in c("extracted", "standardised", "imputed", "finemapped")) {
    expect_true(dir.exists(file.path(study_row$extracted_location, sub_dir)))
  }
})

test_that("extract summary_stats: the p value threshold decides which blocks are extracted", {
  local_block()
  thresholds <- list(genome_wide = 5e-8, default = NA, nothing = 1e-300)
  expected_blocks <- list(genome_wide = fixture_block, default = c(fixture_block, second_block), nothing = character())

  for (name in names(thresholds)) {
    study_row <- summary_stats_study(
      summary_stats_gwas(),
      study_name = paste0("unit-threshold-", name),
      p_value_threshold = thresholds[[name]]
    )
    result <- run_summary_stats(study_row)
    expect_step_succeeded(result, extracted_output_file(study_row))
    if (name == "nothing") {
      expect_lt(file.size(extracted_output_file(study_row)), 70)
    } else {
      expect_equal(read_extracted_snps(study_row)$ld_block, expected_blocks[[name]], info = name)
    }
  }
})

test_that("extract summary_stats: columns are renamed from the metadata column_names", {
  local_block()
  gwas <- dplyr::rename(
    summary_stats_gwas(),
    chromosome = CHR, position = BP, effect_allele = EA, other_allele = OA,
    pval = P, beta = BETA, stderr = SE, freq = EAF
  )
  column_names <- paste0(
    '{"CHR":"chromosome","BP":"position","EA":"effect_allele","OA":"other_allele","P":"pval",',
    '"BETA":"beta","SE":"stderr","EAF":"freq","OR":null,"RSID":""}'
  )
  study_row <- summary_stats_study(gwas, column_names = column_names, file_type = "csv")

  result <- run_summary_stats(study_row)
  expect_step_succeeded(result, extracted_output_file(study_row))

  extracted_snps <- read_extracted_snps(study_row)
  expect_equal(extracted_snps$ld_block, c(fixture_block, second_block))
  extraction <- read_tsv(extracted_snps$file[1])
  expect_true(all(c("CHR", "BP", "EA", "OA", "P", "BETA", "SE", "EAF") %in% colnames(extraction)))
  expect_false(any(c("chromosome", "position", "pval") %in% colnames(extraction)))
})

test_that("extract summary_stats: a numeric N in column_names sets the sample size column", {
  local_block()
  column_names <- '{"CHR":"CHR","BP":"BP","EA":"EA","OA":"OA","P":"P","BETA":"BETA","SE":"SE","EAF":"EAF","N":5000}'
  study_row <- summary_stats_study(summary_stats_gwas(), column_names = column_names)

  result <- run_summary_stats(study_row)
  expect_step_succeeded(result, extracted_output_file(study_row))
  extraction <- read_tsv(read_extracted_snps(study_row)$file[1])
  expect_true("N" %in% colnames(extraction))
  expect_true(all(extraction$N == 5000))
})

test_that("extract summary_stats: chromosome names, P of 0 and odds ratios are standardised", {
  local_block()
  gwas <- summary_stats_gwas()
  gwas$CHR <- paste0("chr", gwas$CHR)
  gwas$P[gwas$P == 1e-5] <- 0
  gwas <- dplyr::mutate(gwas, OR = exp(BETA), OR_LB = exp(BETA - 1.959964 * SE), OR_UB = exp(BETA + 1.959964 * SE)) |>
    dplyr::select(-BETA, -SE)
  other_chromosomes <- data.frame(
    CHR = c("X", "Y", "MT", "chrX", "foo"), BP = 1500000, EA = "A", OA = "G", EAF = 0.3, P = 1e-50,
    OR = 2, OR_LB = 1.5, OR_UB = 2.5
  )
  gwas <- dplyr::bind_rows(gwas, other_chromosomes)
  column_names <- '{"CHR":"CHR","BP":"BP","EA":"EA","OA":"OA","P":"P","OR":"OR","EAF":"EAF"}'
  study_row <- summary_stats_study(gwas, column_names = column_names)

  result <- run_summary_stats(study_row)
  expect_step_succeeded(result, extracted_output_file(study_row))

  extracted_snps <- read_extracted_snps(study_row)
  # the P = 0 hit is treated as the smallest double, so it becomes the most significant region
  expect_equal(extracted_snps$ld_block, c(second_block, fixture_block))
  expect_equal(extracted_snps$log_p[1], -log10(.Machine$double.xmin))
  for (file in extracted_snps$file) {
    extraction <- read_tsv(file)
    expect_true(all(extraction$CHR == 1))
    expect_true(all(c("BETA", "SE") %in% colnames(extraction)))
  }
  # the fixture's last SNP sits on the boundary with the next block, so it is extracted into that block instead
  truth <- read_fixture_gwas("study_a_imputed") |> dplyr::filter(BP < block_bounds(fixture_block)$stop)
  block_1_extraction <- read_tsv(extracted_snps$file[2]) |>
    dplyr::inner_join(dplyr::select(truth, BP, EA, OA, BETA_truth = BETA, SE_truth = SE), by = c("BP", "EA", "OA"))
  expect_equal(nrow(block_1_extraction), nrow(truth))
  expect_equal(block_1_extraction$BETA, block_1_extraction$BETA_truth, tolerance = 1e-6)
  expect_equal(block_1_extraction$SE, block_1_extraction$SE_truth, tolerance = 1e-4)
})

test_that("extract summary_stats: GRCh37 GWASes are lifted over to GRCh38", {
  local_block()
  gwas_38 <- summary_stats_gwas()

  # build a GRCh37 version of the GWAS with the liftOver binary, the reverse of what the extractor does
  bed_38 <- file.path(input_dir(), "hg38.bed")
  bed_37 <- file.path(input_dir(), "hg19.bed")
  vroom::vroom_write(
    data.frame(chr = "chr1", start = gwas_38$BP, end = gwas_38$BP + 1, name = gwas_38$BP),
    bed_38, col_names = FALSE
  )
  system2(
    file.path(liftover_dir, "liftOver"),
    c(bed_38, available_liftover_conversions$GRCh38GRCh37, bed_37, file.path(input_dir(), "unmapped.bed")),
    stdout = FALSE, stderr = FALSE
  )
  lifted <- read_tsv(bed_37, col_names = c("chr", "start", "end", "bp_38")) |> dplyr::distinct(bp_38, .keep_all = TRUE)
  gwas_37 <- dplyr::inner_join(gwas_38, dplyr::select(lifted, bp_38, BP_37 = start), by = c("BP" = "bp_38")) |>
    dplyr::mutate(BP = BP_37) |>
    dplyr::select(-BP_37)
  unmappable <- gwas_37[1, ] |> dplyr::mutate(BP = 249000000)
  gwas_37 <- dplyr::bind_rows(gwas_37, unmappable)

  study_row <- summary_stats_study(gwas_37, reference_build = reference_builds$GRCh37)
  result <- run_summary_stats(study_row)
  expect_step_succeeded(result, extracted_output_file(study_row))
  expect_true(any(grepl("the GWAS lost 1 rows", result$output)))

  extracted_snps <- read_extracted_snps(study_row)
  expect_equal(extracted_snps$ld_block, c(fixture_block, second_block))
  expect_equal(extracted_snps$bp[1], causal_snps$bp[causal_snps$name == "causal_1"])
  block_1_extraction <- read_tsv(extracted_snps$file[1])
  expect_true(all(block_1_extraction$BP %in% gwas_38$BP))
})

test_that("extract summary_stats: VCF files and missing studies fail", {
  local_block()
  study_row <- summary_stats_study(summary_stats_gwas(), file_type = "vcf")
  result <- run_summary_stats(study_row)
  expect_step_failed(result, extracted_output_file(study_row), "VCF support not yet implemented")

  local_block()
  study_row <- summary_stats_study(summary_stats_gwas())
  write_studies_to_process(dplyr::mutate(study_row, extracted_location = "/somewhere/else/"))
  result <- run_extraction("extract_regions_from_summary_stats.R", study_row)
  expect_step_failed(result, extracted_output_file(study_row), "cant find study to process")

  result <- run_step("extract_regions_from_summary_stats.R")
  expect_step_failed(result, NULL, "--worker_guid or both --extracted_study_location")
})

test_that("extract summary_stats: worker uploads are read from study_metadata.json", {
  local_block()
  guid <- "unit-worker-guid"
  worker_dir <- file.path(gwas_upload_dir, "gwas_upload", guid)
  unlink(worker_dir, recursive = TRUE)
  dir.create(worker_dir, recursive = TRUE)
  withr::defer(unlink(worker_dir, recursive = TRUE))

  gwas_file <- file.path(worker_dir, "upload.tsv.gz")
  vroom::vroom_write(summary_stats_gwas(), gwas_file)
  jsonlite::write_json(list(
    guid = guid, file_location = gwas_file, ancestry = "EUR", sample_size = fixture_sample_size,
    category = "continuous", file_type = "csv", reference_build = "GRCh38", p_value_threshold = 1.5e-4,
    column_names = jsonlite::fromJSON(standard_column_names)
  ), file.path(worker_dir, "study_metadata.json"), auto_unbox = TRUE)

  result <- run_step("extract_regions_from_summary_stats.R", "--worker_guid", guid)
  expect_step_succeeded(result, file.path(worker_dir, "extracted_snps.tsv"))
  expect_equal(read_tsv(file.path(worker_dir, "extracted_snps.tsv"))$ld_block, c(fixture_block, second_block))
  expect_true(any(grepl(paste("from", guid), result$output)))
  expect_false(dir.exists(file.path(worker_dir, "svgs")))
})

test_that("extract summary_stats: a SNP on the boundary of two LD blocks is extracted into the block it starts", {
  local_block()
  gwas <- summary_stats_gwas()
  boundary <- block_bounds(second_block)$start
  boundary_hit <- data.frame(CHR = 1, BP = boundary, EA = "A", OA = "G", EAF = 0.3, BETA = 1, SE = 0.1, P = 1e-30)
  gwas <- dplyr::bind_rows(gwas, boundary_hit)
  study_row <- summary_stats_study(gwas)

  result <- run_summary_stats(study_row)
  expect_step_succeeded(result, extracted_output_file(study_row))
  extracted_snps <- read_extracted_snps(study_row)
  expect_equal(sum(extracted_snps$bp == boundary), 1)
  expect_equal(extracted_snps$ld_block[extracted_snps$bp == boundary], second_block)
  expect_equal(anyDuplicated(extracted_snps$file), 0)
  expect_extractions_inside_blocks(extracted_snps)
  block_1_extraction <- read_tsv(extracted_snps$file[extracted_snps$ld_block == fixture_block])
  expect_false(boundary %in% block_1_extraction$BP)
})

test_that("extract summary_stats: a GWAS without P or BETA (e.g. only LOG_P or Z) fails", {
  local_block()
  log_p_only <- dplyr::mutate(summary_stats_gwas(), LOG_P = -log10(P)) |> dplyr::select(-P)
  column_names <- '{"CHR":"CHR","BP":"BP","EA":"EA","OA":"OA","BETA":"BETA","SE":"SE","EAF":"EAF"}'
  study_row <- summary_stats_study(log_p_only, study_name = "unit-log-p-only", column_names = column_names)
  result <- run_summary_stats(study_row)
  expect_step_failed(result, extracted_output_file(study_row))

  local_block()
  z_only <- dplyr::mutate(summary_stats_gwas(), Z = BETA / SE) |> dplyr::select(-BETA)
  column_names <- '{"CHR":"CHR","BP":"BP","EA":"EA","OA":"OA","P":"P","SE":"SE","EAF":"EAF"}'
  study_row <- summary_stats_study(z_only, study_name = "unit-z-only", column_names = column_names)
  result <- run_summary_stats(study_row)
  expect_step_failed(result, extracted_output_file(study_row))
})

test_that("extract summary_stats: an unsupported reference build gives a clear error", {
  local_block()
  study_row <- summary_stats_study(summary_stats_gwas(), reference_build = "GRCh36")
  result <- run_summary_stats(study_row)
  expect_step_failed(result, extracted_output_file(study_row), "liftOver combination of GRCh36 GRCh38 not recocognised")
})

# organise ---------------------------------------------------------------------------------------------------------

test_that("organise: extractions are merged into each LD block's extracted_studies.tsv", {
  local_block()
  study_row <- summary_stats_study(summary_stats_gwas())
  empty_row <- summary_stats_study(summary_stats_gwas(), study_name = "unit-empty", p_value_threshold = 1e-300)
  write_studies_to_process(dplyr::bind_rows(study_row, empty_row))
  expect_step_succeeded(run_extraction("extract_regions_from_summary_stats.R", study_row))
  expect_step_succeeded(run_extraction("extract_regions_from_summary_stats.R", empty_row))

  existing <- read_tsv(extracted_output_file(study_row))[1, ]
  existing_studies <- data.frame(
    study = "existing-study", file = "/existing/file.tsv.gz", ancestry = "EUR", chr = 1, bp = 1000000,
    p_value_threshold = 1.5e-4, category = "continuous", sample_size = 1000, cis_trans = NA,
    reference_build = "GRCh38", ld_block = existing$ld_block, variant_type = "common", coverage = "dense"
  )
  vroom::vroom_write(existing_studies, ld_block_file_paths(fixture_block)$extracted_studies)

  output_file <- file.path(pipeline_metadata_dir, "updated_ld_blocks_to_colocalise.tsv")
  for (run in 1:2) {
    result <- run_step("organise_extracted_regions_into_ld_blocks.R", "--output_file", output_file)
    expect_step_succeeded(result, output_file)
  }

  block_1_studies <- read_tsv(ld_block_file_paths(fixture_block)$extracted_studies)
  expect_equal(sort(block_1_studies$study), c("existing-study", "unit-summary-stats"))
  block_2_studies <- read_tsv(ld_block_file_paths(second_block)$extracted_studies)
  expect_equal(block_2_studies$study, "unit-summary-stats")
  expect_false("unit-empty" %in% c(block_1_studies$study, block_2_studies$study))

  updated_blocks <- read_tsv(output_file)
  expect_equal(updated_blocks$ld_block, c(fixture_block, second_block))
  expect_equal(updated_blocks$data_dir, as.character(c(
    ld_block_file_paths(fixture_block)$ld_block_data,
    ld_block_file_paths(second_block)$ld_block_data
  )))
})

# besd -------------------------------------------------------------------------------------------------------------

gtex_eqtl_probe <- "ENSG00000186827.11"

besd_study <- function(besd_name, probe = gtex_eqtl_probe, study_name = NULL, ...) {
  if (is.null(study_name)) study_name <- paste0("unit-besd-", besd_name)
  return(studies_to_process_row(
    study_name, data_formats$besd, fixture_path("input_studies/besd", besd_name),
    probe = probe, data_type = data_types$gene_expression, ...
  ))
}

run_besd <- function(study_row) {
  write_studies_to_process(study_row)
  return(run_extraction("extract_regions_from_besd.R", study_row))
}

test_that("extract besd: a cis study extracts the LD block of the top cis hit", {
  local_block()
  study_row <- besd_study("GTEx-eQTL-tiny-cis", p_value_threshold = 5e-8)

  result <- run_besd(study_row)
  expect_step_succeeded(result, extracted_output_file(study_row))

  extracted_snps <- read_extracted_snps(study_row)
  expect_equal(nrow(extracted_snps), 1)
  expect_equal(extracted_snps$cis_trans, cis_trans$cis_only)
  expect_equal(extracted_snps$ld_block, second_block)
  expect_extractions_inside_blocks(extracted_snps)

  extraction <- read_tsv(extracted_snps$file)
  expect_true(all(c("RSID", "CHR", "BP", "EA", "OA", "EAF", "BETA", "SE", "P") %in% colnames(extraction)))
  expect_false(any(c("Probe", "Probe_Chr", "Probe_bp", "Gene", "Orientation") %in% colnames(extraction)))
  expect_equal(extracted_snps$log_p, -log10(min(extraction$P)), tolerance = 1e-6)
  expect_true(file.exists(file.path(study_row$extracted_location, "clumped_snps.tsv")))
})

# the tiny GTEx BESD only has hits outside the cis block below p = 1e-2
test_that("extract besd: cis_trans studies extract trans regions outside the cis block", {
  local_block()
  study_row <- besd_study("GTEx-eQTL-tiny-cis_trans", p_value_threshold = 1e-2)

  result <- run_besd(study_row)
  expect_step_succeeded(result, extracted_output_file(study_row))

  extracted_snps <- read_extracted_snps(study_row)
  cis_rows <- extracted_snps[extracted_snps$cis_trans == cis_trans$cis_only, ]
  trans_rows <- extracted_snps[extracted_snps$cis_trans == cis_trans$trans_only, ]
  expect_equal(nrow(cis_rows), 1)
  expect_gte(nrow(trans_rows), 1)
  expect_false(any(trans_rows$ld_block == cis_rows$ld_block))
  expect_equal(anyDuplicated(extracted_snps$ld_block), 0)
  expect_extractions_inside_blocks(extracted_snps)
  clumped_snps <- read_tsv(file.path(study_row$extracted_location, "clumped_snps.tsv"))
  expect_true(all(trans_rows$ld_block %in% clumped_snps$ld_block_string))
})

test_that("extract besd: trans only studies do not extract a cis region", {
  local_block()
  study_row <- besd_study("GTEx-eQTL-tiny-trans", p_value_threshold = 1e-2)

  result <- run_besd(study_row)
  expect_step_succeeded(result, extracted_output_file(study_row))
  extracted_snps <- read_extracted_snps(study_row)
  expect_gte(nrow(extracted_snps), 1)
  expect_true(all(extracted_snps$cis_trans == cis_trans$trans_only))
})

test_that("extract besd: thresholds and missing probes", {
  local_block()
  nothing <- besd_study("GTEx-eQTL-tiny-cis", study_name = "unit-besd-nothing", p_value_threshold = 1e-300)
  expect_step_succeeded(run_besd(nothing), extracted_output_file(nothing))
  expect_equal(nrow(read_extracted_snps(nothing)), 0)

  default_threshold <- besd_study("GTEx-eQTL-tiny-cis", study_name = "unit-besd-default", p_value_threshold = NA)
  expect_step_succeeded(run_besd(default_threshold), extracted_output_file(default_threshold))
  expect_equal(nrow(read_extracted_snps(default_threshold)), 1)

  missing_probe <- besd_study("GTEx-eQTL-tiny-cis", probe = "ENSG00000000000.1", study_name = "unit-besd-no-probe")
  expect_step_succeeded(run_besd(missing_probe), extracted_output_file(missing_probe))
  expect_equal(nrow(read_extracted_snps(missing_probe)), 0)
})

test_that("extract besd: studies with more than one probe are extracted separately", {
  local_block()
  probes <- c(
    "chr1:1204236:1205370:clu_63_-:ENSG00000186891.14",
    "chr1:1204316:1205370:clu_63_-:ENSG00000186891.14"
  )
  study_rows <- dplyr::bind_rows(lapply(seq_along(probes), function(i) {
    return(besd_study("GTEx-sQTL-tiny-cis", probe = probes[i], study_name = paste0("unit-sqtl-", i)))
  }))
  write_studies_to_process(study_rows)

  for (i in seq_len(nrow(study_rows))) {
    result <- run_extraction("extract_regions_from_besd.R", study_rows[i, ])
    expect_step_succeeded(result, extracted_output_file(study_rows[i, ]))
    extracted_snps <- read_extracted_snps(study_rows[i, ])
    expect_equal(nrow(extracted_snps), 1)
    expect_true(startsWith(extracted_snps$file, study_rows$extracted_location[i]))
  }
})

test_that("extract besd: only GRCh38 BESD files are supported", {
  local_block()
  study_row <- besd_study("GTEx-eQTL-tiny-cis", reference_build = reference_builds$GRCh37)
  result <- run_besd(study_row)
  expect_step_failed(result, extracted_output_file(study_row), "Only BESD files using GRCh38")
})

# opengwas ---------------------------------------------------------------------------------------------------------

opengwas_study <- function(...) {
  return(studies_to_process_row(
    "ebi-a-GCST90028992", data_formats$opengwas, fixture_path("input_studies/opengwas/ebi-a-GCST90028992"),
    reference_build = reference_builds$GRCh37, ...
  ))
}

run_opengwas <- function(study_row) {
  write_studies_to_process(study_row)
  return(run_extraction("extract_regions_from_opengwas.R", study_row))
}

test_that("extract opengwas: GRCh37 VCFs are lifted over, clumped, and extracted once per LD block", {
  local_block()
  study_row <- opengwas_study(p_value_threshold = 5e-8)

  result <- run_opengwas(study_row)
  expect_step_succeeded(result, extracted_output_file(study_row))
  expect_true(file.exists(file.path(study_row$extracted_location, "vcf/hg38.vcf.gz")))

  extracted_snps <- read_extracted_snps(study_row)
  expect_gte(nrow(extracted_snps), 1)
  expect_equal(anyDuplicated(extracted_snps$ld_block), 0)
  expect_extractions_inside_blocks(extracted_snps)
  clumped_snps <- read_tsv(file.path(study_row$extracted_location, "clumped_snps.tsv"))
  expect_gte(nrow(clumped_snps), nrow(extracted_snps))
  expect_true(all(clumped_snps$P <= 5e-8))
  expect_true(file.exists(file.path(study_row$extracted_location, "svgs/full.zip")))

  # a second run reuses the lifted over VCF
  rerun <- run_opengwas(study_row)
  expect_step_succeeded(rerun, extracted_output_file(study_row))
  expect_true(any(grepl("already converted to hg38", rerun$output)))
  expect_equal(read_extracted_snps(study_row)$ld_block, extracted_snps$ld_block)
})

test_that("extract opengwas: a study with no significant hits extracts nothing", {
  local_block()
  study_row <- opengwas_study(p_value_threshold = 1e-300)

  result <- run_opengwas(study_row)
  expect_step_succeeded(result, extracted_output_file(study_row))
  expect_equal(nrow(read_extracted_snps(study_row)), 0)
})

# rare tsv ---------------------------------------------------------------------------------------------------------

#' Rare variant GWAS: block 1 has a rare 1e-6 hit (plus a rare variant on the other EAF tail, a common 1e-20 hit and
#' non-autosomal rows that are all filtered out), block 2 has only rare variants with no hit
rare_gwas <- function() {
  panel_1 <- read_tsv(ld_block_file_paths(fixture_block)$ld_matrix_tsv)[1:200, ]
  block_1 <- dplyr::transmute(panel_1, CHR, BP, EA, OA, EAF = 0.001, BETA = 0.1, SE = 0.2, P = 0.5)
  block_1$P[50] <- 1e-6
  block_1$EAF[60] <- 0.9995
  block_1$EAF[70] <- 0.2
  block_1$P[70] <- 1e-20
  panel_2 <- read_tsv(ld_block_file_paths(second_block)$ld_matrix_tsv)[1:100, ]
  block_2 <- dplyr::transmute(panel_2, CHR, BP, EA, OA, EAF = 0.002, BETA = 0.1, SE = 0.2, P = 0.5)
  non_autosomal <- data.frame(CHR = "X", BP = 1500000, EA = "A", OA = "G", EAF = 0.001, BETA = 1, SE = 0.1, P = 1e-30)
  autosomal <- dplyr::bind_rows(block_1, block_2) |> dplyr::mutate(CHR = as.character(CHR))
  return(dplyr::bind_rows(autosomal, non_autosomal))
}

rare_study <- function(gwas, study_name = "unit-rare") {
  study_location <- file.path(input_dir(), paste0(study_name, ".tsv.gz"))
  vroom::vroom_write(gwas, study_location)
  return(studies_to_process_row(study_name, data_formats$tsv, study_location, variant_type = variant_types$rare_exome))
}

run_rare <- function(study_row) {
  write_studies_to_process(study_row)
  return(run_extraction("extract_regions_from_rare_tsv.R", study_row))
}

test_that("extract rare tsv: only rare variants are kept, and blocks with a rare hit are extracted", {
  local_block()
  gwas <- rare_gwas()
  study_row <- rare_study(gwas)

  result <- run_rare(study_row)
  expect_step_succeeded(result, extracted_output_file(study_row))

  extracted_snps <- read_extracted_snps(study_row)
  expect_equal(extracted_snps$ld_block, fixture_block)
  expect_equal(extracted_snps$bp, gwas$BP[50])
  extraction <- read_tsv(extracted_snps$file)
  expect_equal(nrow(extraction), 199)
  expect_false(gwas$BP[70] %in% extraction$BP)
  expect_true(gwas$BP[60] %in% extraction$BP)
  expect_true(all(extraction$CHR == 1))
})

test_that("extract rare tsv: odds ratios are converted to BETA and SE", {
  local_block()
  gwas <- rare_gwas() |>
    dplyr::mutate(OR = exp(BETA), CI_LOWER = exp(BETA - 1.96 * SE), CI_UPPER = exp(BETA + 1.96 * SE)) |>
    dplyr::select(-BETA, -SE)
  gwas$OR[40] <- 0
  study_row <- rare_study(gwas)

  result <- run_rare(study_row)
  expect_step_succeeded(result, extracted_output_file(study_row))

  extraction <- read_tsv(read_extracted_snps(study_row)$file)
  expect_equal(extraction$BETA[extraction$BP != gwas$BP[40]][1:5], rep(0.1, 5), tolerance = 1e-6)
  expect_equal(extraction$SE[extraction$BP != gwas$BP[40]][1:5], rep(0.2, 5), tolerance = 1e-6)
  expect_equal(extraction$BETA[extraction$BP == gwas$BP[40]], log(0.01))
  expect_equal(extraction$SE[extraction$BP == gwas$BP[40]], 0)
})

test_that("extract rare tsv: missing columns fail, and a fully filtered GWAS extracts nothing", {
  local_block()
  missing_columns <- rare_study(dplyr::select(rare_gwas(), -SE), study_name = "unit-rare-missing")
  result <- run_rare(missing_columns)
  expect_step_failed(result, extracted_output_file(missing_columns), "Study must have BETA and SE")

  common_only <- rare_study(dplyr::mutate(rare_gwas(), EAF = 0.3), study_name = "unit-rare-common")
  result <- run_rare(common_only)
  expect_step_succeeded(result, extracted_output_file(common_only))
  extracted_snps <- read_extracted_snps(common_only)
  expect_equal(nrow(extracted_snps), 0)
  expect_true(all(c("chr", "bp", "log_p", "ld_block", "file", "cis_trans") %in% colnames(extracted_snps)))
})

test_that("extract rare tsv: significant variants outside every LD block are dropped and logged", {
  local_block()
  gwas <- dplyr::bind_rows(
    rare_gwas(),
    data.frame(CHR = "1", BP = 10, EA = "A", OA = "G", EAF = 0.001, BETA = 1, SE = 0.1, P = 1e-30)
  )
  study_row <- rare_study(gwas)

  result <- run_rare(study_row)
  expect_step_succeeded(result, extracted_output_file(study_row))
  expect_equal(read_extracted_snps(study_row)$ld_block, fixture_block)

  missing_ld_blocks <- read_tsv(
    file.path(pipeline_metadata_dir, "missing_ld_blocks.tsv"),
    col_names = c("study", "chr", "bp"),
    delim = "\t"
  )
  expect_equal(missing_ld_blocks$study, "unit-rare")
  expect_equal(missing_ld_blocks$bp, 10)
})
