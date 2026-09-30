#!/usr/bin/env Rscript
# Repairs the SE=0.00001 sentinel bug introduced in
# standardise_studies_in_ld_block.R:293 (fixed in the accompanying code
# change), in two parts:
#
#   1. Patches associations.se directly in an already-built associations
#      DuckDB file (associations_full.db / associations_specific.db). This
#      alone fixes what the API/gpmapr serve - beta/se/p there are never
#      touched by finemapping, so no pipeline rerun is needed for it.
#
#   2. For (study, ld_block) pairs whose SuSiE input Z (= BETA / SE) was
#      corrupted by the same sentinel, patches SE in the per-study imputed
#      GWAS files (what finemap_rule actually reads) and strips those
#      studies out of the finemapped_studies.tsv / coloc_pairwise_results
#      tracking files, so a subsequent pipeline run only recomputes the
#      affected studies instead of skipping them as "already done".
#
# This script only edits/deletes files - it never invokes snakemake,
# run_pipeline.sh, or any Makefile target. Triggering an actual pipeline
# run is a separate, deliberate step you take afterwards.
#
# Usage:
#   Rscript scripts/fix_se_sentinel.R <associations_db_file> [<studies_db_file>] [--dry-run]
#
# associations_db_file: associations_full.db or associations_specific.db
# studies_db_file (optional): studies.db - used to resolve affected
#                              (study_name, ld_block) pairs for step 2. If
#                              omitted, step 2 is skipped.

source("../pipeline_steps/constants.R")

SENTINEL <- 1e-05
TOLERANCE <- 1e-9 # float32 round-trip noise around the sentinel

back_derive_se <- function(beta, p) {
  # Invert the two-tailed Wald p-value formula: se = |beta| / qnorm(1 - p/2)
  se <- abs(beta) / qnorm(1 - p / 2)
  se[!is.finite(se) | se <= 0] <- NA
  return(se)
}

# --- Step 1: patch associations.se in the final DuckDB file -----------------

repair_associations_db <- function(associations_db_file, dry_run) {
  con <- duckdb::dbConnect(duckdb::duckdb(), associations_db_file)
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE))

  candidates <- DBI::dbGetQuery(con, glue::glue(
    "SELECT variant_id, study_id, beta, se, p FROM associations
     WHERE abs(se - {SENTINEL}) < {TOLERANCE}"
  ))
  message(nrow(candidates), " rows carry the SE sentinel in ", associations_db_file)
  if (nrow(candidates) == 0) {
    return(candidates)
  }

  candidates$new_se <- back_derive_se(candidates$beta, candidates$p)
  unfixable <- is.na(candidates$new_se)
  if (any(unfixable)) {
    message(
      sum(unfixable),
      " rows can't be repaired from beta/p (beta==0 and/or p in {0,1,NA}) - left as-is"
    )
    vroom::vroom_write(candidates[unfixable, ], "unfixable_se_sentinel_rows.tsv")
    candidates <- candidates[!unfixable, ]
  }

  vroom::vroom_write(candidates, "se_sentinel_repairs_audit.tsv")
  message("Repairing ", nrow(candidates), " rows; audit written to se_sentinel_repairs_audit.tsv")

  if (!dry_run) {
    DBI::dbExecute(
      con,
      "CREATE OR REPLACE TEMP TABLE se_fixes (variant_id INTEGER, study_id INTEGER, new_se DOUBLE)"
    )
    DBI::dbAppendTable(con, "se_fixes", candidates[, c("variant_id", "study_id", "new_se")])
    DBI::dbExecute(con, "
      UPDATE associations
      SET se = se_fixes.new_se
      FROM se_fixes
      WHERE associations.variant_id = se_fixes.variant_id
        AND associations.study_id = se_fixes.study_id
    ")
    message("associations.se repaired in place in ", associations_db_file)
  } else {
    message("--dry-run: no changes written to ", associations_db_file)
  }

  return(candidates)
}

# --- Resolve which (study_name, ld_block) pairs need a finemap/coloc redo ---
# min_p in study_extractions/coloc_groups comes from feeding Z into SuSiE and
# converting LBF back to a p-value (gwas_calculations.R:94-103,
# finemap_studies_in_ld_block.R:321,467) - not invertible by formula.

resolve_affected_study_ld_blocks <- function(studies_db_file, affected_study_ids) {
  con <- duckdb::dbConnect(duckdb::duckdb(), studies_db_file, read_only = TRUE)
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE))

  affected <- DBI::dbGetQuery(con, glue::glue(
    "SELECT DISTINCT study AS study_name, ld_block
     FROM study_extractions
     WHERE study_id IN ({paste(affected_study_ids, collapse = ',')})"
  ))
  vroom::vroom_write(affected, "affected_study_ld_blocks_for_rerun.tsv")
  message(
    nrow(affected),
    " (study_name, ld_block) pairs need finemap+coloc rerun -> affected_study_ld_blocks_for_rerun.tsv"
  )
  return(affected)
}

# --- Step 2: patch imputed data + strip resume-tracking rows ----------------

prepare_ld_blocks_for_refinemap_recoloc <- function(affected_pairs, dry_run) {
  if (nrow(affected_pairs) == 0) {
    return(invisible())
  }

  for (block in split(affected_pairs, affected_pairs$ld_block)) {
    ld_block <- block$ld_block[[1]]
    affected_studies <- unique(block$study_name)
    paths <- ld_block_file_paths(ld_block)
    message(glue::glue("{ld_block}: preparing {length(affected_studies)} studies for refinemap/recoloc"))

    # 1. Patch SE in the per-study imputed GWAS files - this is what
    #    finemap_rule actually reads (finemap_studies_in_ld_block.R:62,121).
    if (file.exists(paths$imputed_studies)) {
      imputed_studies <- vroom::vroom(paths$imputed_studies, show_col_types = FALSE)
      relevant <- imputed_studies[imputed_studies$study %in% affected_studies, , drop = FALSE]
      for (i in seq_len(nrow(relevant))) {
        study_file <- relevant$file[[i]]
        if (!file.exists(study_file)) next
        gwas <- vroom::vroom(study_file, show_col_types = FALSE)
        zero_se <- abs(gwas$SE - SENTINEL) < TOLERANCE
        if (!any(zero_se)) next
        derived_se <- back_derive_se(gwas$BETA[zero_se], gwas$P[zero_se])
        derived_se[is.na(derived_se)] <- SENTINEL # leave genuinely unfixable rows as-is
        if (!dry_run) {
          gwas$SE[zero_se] <- derived_se
          vroom::vroom_write(gwas, study_file)
        }
        message(glue::glue("  {relevant$study[[i]]}: patched {sum(zero_se)} SE values in {study_file}"))
      }
    }

    # 2. Strip affected studies out of finemapped_studies.tsv, so
    #    finemap_studies_in_ld_block.R:79-85 no longer treats them as done,
    #    and remember which unique_study_id values that removes.
    removed_unique_study_ids <- character(0)
    if (file.exists(paths$finemapped_studies)) {
      finemapped <- vroom::vroom(
        paths$finemapped_studies,
        show_col_types = FALSE, col_types = finemapped_column_types
      )
      to_remove <- finemapped[finemapped$study %in% affected_studies, , drop = FALSE]
      removed_unique_study_ids <- unique(to_remove$unique_study_id)
      if (!dry_run && nrow(to_remove) > 0) {
        remaining <- finemapped[!finemapped$study %in% affected_studies, , drop = FALSE]
        vroom::vroom_write(remaining, paths$finemapped_studies)
      }
      message(
        glue::glue(
          "  removed {nrow(to_remove)} rows ({length(removed_unique_study_ids)} unique_study_ids) from finemapped_studies.tsv"
        )
      )
    }

    # 3. Delete the now-stale per-study finemap output files.
    flattened_block_name <- flattened_ld_block_name(ld_block)
    for (study in affected_studies) {
      stale_files <- Sys.glob(glue::glue(
        "{extracted_study_dir}/{study}/finemapped/{flattened_block_name}*"
      ))
      if (length(stale_files) > 0 && !dry_run) {
        file.remove(stale_files)
      }
      if (length(stale_files) > 0) {
        message(glue::glue("  removed {length(stale_files)} stale finemap output files for {study}"))
      }
    }

    # 4. Strip any coloc pair touching a removed unique_study_id, so
    #    coloc_studies_in_ld_block.R:117-126 treats those pairs as new again.
    if (length(removed_unique_study_ids) > 0 && file.exists(paths$coloc_pairwise)) {
      coloc_results <- vroom::vroom(
        paths$coloc_pairwise,
        show_col_types = FALSE, col_types = coloc_pairwise_results_column_types
      )
      stale_pairs <- coloc_results$unique_study_a %in% removed_unique_study_ids |
        coloc_results$unique_study_b %in% removed_unique_study_ids
      if (!dry_run && any(stale_pairs)) {
        vroom::vroom_write(coloc_results[!stale_pairs, , drop = FALSE], paths$coloc_pairwise)
      }
      message(glue::glue("  removed {sum(stale_pairs)} stale pairs from coloc_pairwise_results.tsv.gz"))
    }

    # 5. Delete completion markers so Snakemake's DAG sees this ld_block as
    #    needing finemap_rule/coloc_rule/clustering_rule/compare_rare_rule
    #    again. standardisation_complete/imputation_complete are left alone -
    #    those stages don't need to be redone.
    markers <- c(
      paths$finemapping_complete, paths$coloc_complete,
      paths$clustering_complete, paths$compare_rare_complete
    )
    markers <- markers[file.exists(markers)]
    if (!dry_run) {
      file.remove(markers)
    }
    message(glue::glue("  removed {length(markers)} completion markers"))
  }

  message(
    "Done. This script did not run the pipeline - use run_pipeline.sh/snakemake ",
    "yourself once you've reviewed the changes above."
  )
  return(invisible())
}

# --- Entry point --------------------------------------------------------

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  dry_run <- "--dry-run" %in% args
  positional <- args[!args %in% c("--dry-run")]
  associations_db_file <- positional[[1]]
  studies_db_file <- if (length(positional) > 1) positional[[2]] else NA

  stopifnot(file.exists(associations_db_file))

  candidates <- repair_associations_db(associations_db_file, dry_run)
  if (nrow(candidates) == 0) {
    return(invisible())
  }

  affected_study_ids <- unique(candidates$study_id)
  if (is.na(studies_db_file) || !file.exists(studies_db_file)) {
    writeLines(as.character(affected_study_ids), "affected_study_ids_for_rerun.txt")
    message(
      length(affected_study_ids),
      " affected study_ids written to affected_study_ids_for_rerun.txt",
      " (pass studies_db_file to also prepare ld_blocks for refinemap/recoloc)"
    )
    return(invisible())
  }

  affected_pairs <- resolve_affected_study_ld_blocks(studies_db_file, affected_study_ids)
  prepare_ld_blocks_for_refinemap_recoloc(affected_pairs, dry_run)
  return(invisible())
}

main()
