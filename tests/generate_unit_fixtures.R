# One-shot generator for the unit (per-step) test fixtures.
#
# Run once, inside the apptainer image, from the repo root:
#   Rscript tests/generate_unit_fixtures.R
#
# It builds /local-scratch/projects/genotype-phenotype-map/test/unit/:
#   fixtures/  read-only inputs used by tests/testthat/unit/*.R
#   data/      DATA_DIR for the unit tests (reference panel links, liftover, variant annotations)
#   results/   RESULTS_DIR for the unit tests
# It needs the e2e test data (and a completed e2e run, for the create_db fixture) to already exist.
# It is not run as part of the test suite.

source("pipeline_steps/common_extraction_functions.R")
source("pipeline_steps/gwas_calculations.R")

study_categories <- list(continuous = "continuous", categorical = "categorical")

project_dir <- "/local-scratch/projects/genotype-phenotype-map"
unit_dir <- file.path(project_dir, "test/unit")
e2e_data_dir <- file.path(project_dir, "test/e2e/data")
e2e_results_dir <- file.path(project_dir, "test/e2e/results/current")
real_panel_dir <- file.path(project_dir, "data/ld_reference_panel_hg38/EUR")

fixtures_dir <- file.path(unit_dir, "fixtures")
data_dir <- file.path(unit_dir, "data")
results_dir <- file.path(unit_dir, "results")
panel_dir <- file.path(data_dir, "ld_reference_panel_hg38/EUR")

fixture_block <- "16103-1170341"
panel_blocks <- c("16103-1170341", "1170341-1730405", "1730405-3355587", "3355587-4320284")
max_bp <- 4400000
panel_sample_size <- 2000
sample_size <- 100000

main <- function() {
  for (dir in c(
    fixtures_dir, file.path(fixtures_dir, "gwas"), file.path(fixtures_dir, "finemapped"),
    file.path(fixtures_dir, "reference_panel"), file.path(fixtures_dir, "input_studies/besd"),
    file.path(fixtures_dir, "input_studies/opengwas"), file.path(panel_dir, "1"),
    file.path(data_dir, "variant_annotation"), file.path(data_dir, "study"), file.path(data_dir, "ld_blocks"),
    file.path(data_dir, "pipeline_metadata"), file.path(results_dir, "current"), file.path(results_dir, "latest")
  )) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  }

  build_reference_panel()
  build_variant_annotations()
  simulate_gwas_fixtures()
  build_besd_fixtures()
  build_opengwas_fixture()
  build_create_db_fixture()
  message("Unit fixtures written to ", unit_dir)
  return(invisible(NULL))
}

#' Links the LD matrices of the test blocks, and writes small (chr1 < 4.4Mb, 2000 samples) bfiles so that
#' clumping and DENTIST run in seconds
build_reference_panel <- function() {
  for (block in panel_blocks) {
    for (suffix in c(".tsv", ".unphased.vcor1.gz", ".unphased.vcor1.vars", ".ldeig.rds")) {
      link_file(file.path(real_panel_dir, "1", paste0(block, suffix)), file.path(panel_dir, "1", paste0(block, suffix)))
    }
  }
  link_file(file.path(project_dir, "data/liftover"), file.path(data_dir, "liftover"))

  keep_file <- file.path(fixtures_dir, "reference_panel/keep_samples.txt")
  fam <- vroom::vroom(file.path(real_panel_dir, "full.fam"), col_names = FALSE, show_col_types = FALSE)
  vroom::vroom_write(fam[seq_len(panel_sample_size), 1:2], keep_file, col_names = FALSE)

  subsets <- list(
    full = file.path(real_panel_dir, "full"),
    full_rsid = file.path(real_panel_dir, "full_rsid"),
    fixture_block = file.path(real_panel_dir, "1", fixture_block)
  )
  for (name in names(subsets)) {
    out_prefix <- file.path(fixtures_dir, "reference_panel", name)
    if (!file.exists(paste0(out_prefix, ".bed"))) {
      plink_command <- glue::glue(
        "plink2 --bfile {subsets[[name]]} --chr 1 --to-bp {max_bp} --keep {keep_file} ",
        "--make-bed --out {out_prefix}"
      )
      run_command(plink_command)
    }
  }
  for (ext in c(".bed", ".bim", ".fam")) {
    link_file(
      file.path(fixtures_dir, "reference_panel", paste0("full", ext)),
      file.path(panel_dir, paste0("full", ext))
    )
    link_file(
      file.path(fixtures_dir, "reference_panel", paste0("full_rsid", ext)),
      file.path(panel_dir, paste0("full_rsid", ext))
    )
    link_file(
      file.path(fixtures_dir, "reference_panel", paste0("fixture_block", ext)),
      file.path(panel_dir, "1", paste0(fixture_block, ext))
    )
  }
  return(invisible(NULL))
}

#' All columns of the real VEP annotations, for chr1 < 4.4Mb only
build_variant_annotations <- function() {
  vep_file <- file.path(data_dir, "variant_annotation/vep_annotations_hg38.tsv.gz")
  if (!file.exists(vep_file)) {
    run_command(glue::glue(
      "zcat {e2e_data_dir}/variant_annotation/vep_annotations_hg38.tsv.gz | ",
      "awk -F'\\t' 'NR == 1 || ($1 == 1 && $2 <= {max_bp})' | gzip > {vep_file}"
    ))
  }
  link_file(
    file.path(e2e_data_dir, "variant_annotation/gene_info.tsv"),
    file.path(data_dir, "variant_annotation/gene_info.tsv")
  )
  return(invisible(NULL))
}

#' Simulates GWASes over the fixture LD block (z = R * lambda + MVN(0, R) noise), and checks SuSiE finds the
#' expected number of credible sets for each one, so each fixture always hits its intended code path
simulate_gwas_fixtures <- function() {
  ld_info <- vroom::vroom(file.path(panel_dir, "1", paste0(fixture_block, ".tsv")), show_col_types = FALSE)
  ld_matrix <- as.matrix(vroom::vroom(
    file.path(panel_dir, "1", paste0(fixture_block, ".unphased.vcor1.gz")),
    col_names = FALSE, show_col_types = FALSE, altrep = FALSE
  ))
  dimnames(ld_matrix) <- NULL
  eigen_ld <- eigen(ld_matrix, symmetric = TRUE)
  vep <- vroom::vroom(
    file.path(data_dir, "variant_annotation/vep_annotations_hg38.tsv.gz"),
    col_select = c("snp", "rsid"), show_col_types = FALSE
  )

  candidates <- which(ld_info$EAF > 0.1 & ld_info$EAF < 0.9 & nchar(ld_info$EA) == 1 & nchar(ld_info$OA) == 1)
  low_ld <- function(i, chosen) all(abs(ld_matrix[i, chosen]) < 0.05)

  # study A: 2 signals, far apart; study C: a signal within 50kb of A's first signal, but not in LD with it
  causal_1 <- candidates[round(length(candidates) / 3)]
  causal_2 <- candidates[candidates != causal_1 & abs(ld_info$BP[candidates] - ld_info$BP[causal_1]) > 150000]
  causal_2 <- causal_2[vapply(causal_2, low_ld, logical(1), chosen = causal_1)][1]
  causal_3 <- candidates[abs(ld_info$BP[candidates] - ld_info$BP[causal_1]) < 40000 & candidates != causal_1]
  causal_3 <- causal_3[vapply(causal_3, low_ld, logical(1), chosen = c(causal_1, causal_2))][1]
  four_signals <- causal_1
  for (i in candidates) {
    far_enough <- all(abs(ld_info$BP[i] - ld_info$BP[four_signals]) > 60000)
    if (length(four_signals) < 4 && low_ld(i, four_signals) && far_enough) {
      four_signals <- c(four_signals, i)
    }
  }
  stopifnot(!is.na(causal_2), !is.na(causal_3), length(four_signals) == 4)

  study_designs <- list(
    study_a = list(causal = c(causal_1, causal_2), z = c(10, -8), expected_cs = 2),
    study_b = list(causal = causal_1, z = 9, expected_cs = 1),
    study_c = list(causal = causal_3, z = 9, expected_cs = 1),
    study_d = list(causal = four_signals, z = c(11, -10, 9, -9), expected_cs = 4),
    study_null = list(causal = integer(0), z = numeric(0), expected_cs = 0)
  )

  for (seed in 1:200) {
    set.seed(seed)
    simulated <- lapply(study_designs, function(design) {
      lambda <- rep(0, nrow(ld_matrix))
      lambda[design$causal] <- design$z
      noise <- as.vector(eigen_ld$vectors %*% (sqrt(pmax(eigen_ld$values, 0)) * rnorm(nrow(ld_matrix))))
      z <- as.vector(ld_matrix %*% lambda) + noise
      susie_result <- susieR::susie_rss(z = z, R = ld_matrix, n = sample_size)
      return(list(z = z, susie = susie_result))
    })
    cs_counts <- vapply(simulated, function(s) length(s$susie$sets$cs_index), numeric(1))
    expected <- vapply(study_designs, function(d) d$expected_cs, numeric(1))
    cs_ok <- all(cs_counts[names(expected) != "study_d"] == expected[names(expected) != "study_d"]) &&
      cs_counts[["study_d"]] >= 4
    if (!cs_ok) next

    lbfs <- lapply(simulated[c("study_a", "study_b", "study_c")], function(s) {
      lbf <- s$susie$lbf_variable[s$susie$sets$cs_index, , drop = FALSE]
      colnames(lbf) <- ld_info$SNP
      return(lbf)
    })
    cs_1_of_a <- which.max(abs(lbfs$study_a[, causal_1]))
    h4_shared <- coloc::coloc.bf_bf(lbfs$study_a[cs_1_of_a, ], lbfs$study_b[1, ], p1 = 1e-4, p2 = 1e-4, p12 = 5e-6)
    h4_distinct <- coloc::coloc.bf_bf(lbfs$study_a[cs_1_of_a, ], lbfs$study_c[1, ], p1 = 1e-4, p2 = 1e-4, p12 = 5e-6)
    if (h4_shared$summary$PP.H4.abf > 0.9 && h4_distinct$summary$PP.H4.abf < 0.2) break
  }
  message(glue::glue(
    "Seed {seed}: CS counts {paste(names(cs_counts), cs_counts, collapse = ', ')}; ",
    "H4 shared {signif(h4_shared$summary$PP.H4.abf, 3)}, distinct {signif(h4_distinct$summary$PP.H4.abf, 3)}"
  ))
  stopifnot(cs_ok, h4_shared$summary$PP.H4.abf > 0.9, h4_distinct$summary$PP.H4.abf < 0.2)

  rsids <- vep$rsid[match(ld_info$SNP, vep$snp)]
  gwases <- lapply(simulated, function(s) {
    se <- 1 / sqrt(2 * ld_info$EAF * (1 - ld_info$EAF) * sample_size)
    return(data.frame(
      CHR = ld_info$CHR, SNP = ld_info$SNP, OA = ld_info$OA, EA = ld_info$EA, BP = ld_info$BP,
      RSID = ifelse(is.na(rsids), paste0("rs_unit_", seq_along(rsids)), sub(",.*", "", rsids)),
      EAF = ld_info$EAF, BETA = s$z * se, SE = se, P = 2 * pnorm(-abs(s$z)), Z = s$z
    ))
  })

  for (name in names(gwases)) {
    imputed <- dplyr::mutate(gwases[[name]], IMPUTED = FALSE)
    vroom::vroom_write(imputed, file.path(fixtures_dir, "gwas", glue::glue("{name}_imputed.tsv.gz")))
  }

  # study A also gets a standardised form (25% of SNPs missing, causal SNPs kept), and an extracted form
  # (about half the alleles swapped, plus SNPs not in the reference panel)
  set.seed(1)
  study_a <- gwases$study_a
  droppable <- setdiff(seq_len(nrow(study_a)), c(causal_1, causal_2))
  dropped <- sort(sample(droppable, round(0.25 * nrow(study_a))))
  vroom::vroom_write(study_a[-dropped, ], file.path(fixtures_dir, "gwas", "study_a_standardised.tsv.gz"))
  vroom::vroom_write(
    data.frame(SNP = study_a$SNP[dropped]),
    file.path(fixtures_dir, "gwas", "study_a_standardised_dropped_snps.tsv")
  )

  extracted <- dplyr::select(study_a, RSID, CHR, BP, EA, OA, EAF, BETA, SE, P)
  to_swap <- seq_len(nrow(extracted)) %% 2 == 0
  extracted <- flip_alleles(extracted, to_swap)
  not_in_panel <- extracted[c(10, 20, 30, 40, 50), ] |>
    dplyr::mutate(BP = BP + 1, RSID = paste0("rs_not_in_panel_", dplyr::row_number()))
  extracted <- dplyr::bind_rows(extracted, not_in_panel) |> dplyr::arrange(BP)
  vroom::vroom_write(extracted, file.path(fixtures_dir, "gwas", "study_a_extracted.tsv.gz"))

  # finemapped (LBF) files and a finemapped_studies template for coloc / cluster
  finemapped_rows <- list()
  for (name in c("study_a", "study_b", "study_c")) {
    susie_result <- simulated[[name]]$susie
    gwas <- gwases[[name]]
    for (cs_number in seq_along(susie_result$sets$cs_index)) {
      lbf <- susie_result$lbf_variable[susie_result$sets$cs_index[cs_number], ]
      log_p <- convert_lbf_to_log_p_value(lbf, gwas$SE)
      lead <- which.min(log_p)
      file_name <- glue::glue("{name}_{cs_number}.tsv.gz")
      vroom::vroom_write(
        dplyr::select(gwas, SNP, CHR, BP, BETA, SE, EAF) |> dplyr::mutate(IMPUTED = FALSE, LBF = lbf),
        file.path(fixtures_dir, "finemapped", file_name)
      )
      finemapped_rows[[length(finemapped_rows) + 1]] <- data.frame(
        study = name, cs = cs_number, file_name = file_name, snp = gwas$SNP[lead], bp = gwas$BP[lead],
        min_p = exp(min(log_p))
      )
    }
  }
  vroom::vroom_write(
    dplyr::bind_rows(finemapped_rows),
    file.path(fixtures_dir, "finemapped", "finemapped_studies_template.tsv")
  )

  vroom::vroom_write(
    data.frame(
      name = c("causal_1", "causal_2", "causal_3", paste0("four_signals_", 1:4)),
      snp = ld_info$SNP[c(causal_1, causal_2, causal_3, four_signals)],
      bp = ld_info$BP[c(causal_1, causal_2, causal_3, four_signals)]
    ),
    file.path(fixtures_dir, "gwas", "causal_snps.tsv")
  )
  vroom::vroom_write(
    data.frame(seed = seed, sample_size = sample_size, ld_block = paste0("EUR/1/", fixture_block)),
    file.path(fixtures_dir, "gwas", "simulation_info.tsv")
  )
  return(invisible(NULL))
}

#' Cuts the GTEx eQTL (1 probe) and sQTL (2 probes) BESDs down to chr1 < 4.4Mb
build_besd_fixtures <- function() {
  qtl_dir <- file.path(e2e_data_dir, "input_studies/qtl")
  besd_dir <- file.path(fixtures_dir, "input_studies/besd")
  besds <- list(
    "GTEx-eQTL-tiny" = "GTEx-eQTL-v10-Brain_Cortex",
    "GTEx-sQTL-tiny" = "GTEx-sQTL-v10-Brain_Cortex"
  )
  for (name in names(besds)) {
    source_prefix <- file.path(qtl_dir, besds[[name]])
    snp_list <- file.path(besd_dir, glue::glue("{name}_snps.list"))
    probe_list <- file.path(besd_dir, glue::glue("{name}_probes.list"))
    run_command(glue::glue("awk '$1 == 1 && $4 <= {max_bp} {{print $2}}' {source_prefix}.esi > {snp_list}"))
    run_command(glue::glue("awk '{{print $2}}' {source_prefix}.epi > {probe_list}"))
    run_command(glue::glue(
      "smr --beqtl-summary {source_prefix} --extract-probe {probe_list} --extract-snp {snp_list} ",
      "--make-besd --out {file.path(besd_dir, name)}"
    ))

    metadata <- jsonlite::fromJSON(paste0(source_prefix, ".json"))
    for (cis_trans_value in c("cis", "trans", "cis_trans")) {
      variant_prefix <- file.path(besd_dir, glue::glue("{name}-{cis_trans_value}"))
      for (ext in c(".besd", ".esi", ".epi")) {
        file.copy(file.path(besd_dir, paste0(name, ext)), paste0(variant_prefix, ext), overwrite = TRUE)
      }
      metadata$cis_trans <- cis_trans_value
      jsonlite::write_json(metadata, paste0(variant_prefix, ".json"), auto_unbox = TRUE, pretty = TRUE)
    }
  }
  return(invisible(NULL))
}

build_opengwas_fixture <- function() {
  source_dir <- file.path(e2e_data_dir, "input_studies/opengwas/ebi-a-GCST90028992")
  out_dir <- file.path(fixtures_dir, "input_studies/opengwas/ebi-a-GCST90028992")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  for (file in c("ebi-a-GCST90028992.vcf.gz", "ebi-a-GCST90028992.vcf.gz.tbi", "ebi-a-GCST90028992.json")) {
    file.copy(file.path(source_dir, file), file.path(out_dir, file), overwrite = TRUE)
  }
  return(invisible(NULL))
}

#' Compiled results of the last e2e run, plus the LD block index files and study files they reference,
#' with paths rewritten to point at the unit data dir
build_create_db_fixture <- function() {
  out_dir <- file.path(fixtures_dir, "create_db")
  dir.create(file.path(out_dir, "results"), recursive = TRUE, showWarnings = FALSE)
  rewrite_paths <- function(df) {
    return(dplyr::mutate(df, dplyr::across(
      dplyr::where(is.character),
      \(x) gsub(paste0(project_dir, "/test/(e2e/)?data/"), paste0(data_dir, "/"), x)
    )))
  }

  result_files <- c(
    "studies_processed", "traits_processed", "study_extractions",
    "coloc_clustered_results", "coloc_pairwise_results", "rare_results"
  )
  for (result_file in result_files) {
    results <- vroom::vroom(
      file.path(e2e_results_dir, glue::glue("{result_file}.tsv.gz")),
      show_col_types = FALSE, col_types = vroom::cols(.default = "c")
    ) |> rewrite_paths()
    vroom::vroom_write(results, file.path(out_dir, "results", glue::glue("{result_file}.tsv.gz")), na = "NA")
  }

  for (block_dir in list.dirs(file.path(e2e_data_dir, "ld_blocks/EUR/1"), recursive = FALSE)) {
    out_block_dir <- file.path(out_dir, "data/ld_blocks/EUR/1", basename(block_dir))
    dir.create(out_block_dir, recursive = TRUE, showWarnings = FALSE)
    for (index_file in c("standardised_studies.tsv", "imputed_studies.tsv", "finemapped_studies.tsv")) {
      if (!file.exists(file.path(block_dir, index_file))) next
      index <- vroom::vroom(
        file.path(block_dir, index_file),
        show_col_types = FALSE, col_types = vroom::cols(.default = "c")
      ) |> rewrite_paths()
      vroom::vroom_write(index, file.path(out_block_dir, index_file), na = "NA")
    }
  }

  study_dirs <- list.dirs(file.path(e2e_data_dir, "study"), recursive = FALSE)
  for (study_dir in study_dirs) {
    for (sub_dir in c("standardised", "imputed", "finemapped")) {
      out_study_dir <- file.path(out_dir, "data/study", basename(study_dir), sub_dir)
      dir.create(out_study_dir, recursive = TRUE, showWarnings = FALSE)
      files <- list.files(file.path(study_dir, sub_dir), full.names = TRUE)
      file.copy(files, out_study_dir, overwrite = TRUE)
    }
  }
  return(invisible(NULL))
}

link_file <- function(target, link) {
  if (!file.exists(target)) {
    message("Skipping missing ", target)
    return(invisible(NULL))
  }
  if (file.exists(link) || nzchar(Sys.readlink(link))) file.remove(link)
  file.symlink(target, link)
  return(invisible(NULL))
}

run_command <- function(command) {
  message(command)
  status <- system(command)
  if (status != 0) stop("Command failed: ", command)
  return(invisible(NULL))
}

main()
