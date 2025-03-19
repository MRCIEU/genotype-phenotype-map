source('constants.R')

library(dplyr)
library(duckdb)
library(data.table)
library(validate)
library(tidyr)
library(R.utils)
library(purrr)
library(furrr)
library(parallel)

parser <- argparser::arg_parser('Create DuckDB from pipeline results')
#INPUT
parser <- argparser::add_argument(parser, '--results_dir', help = 'Results directory of pipeline', type = 'character')
#OUTPUT
parser <- argparser::add_argument(parser, '--studies_db_file', help = 'Existing DuckDB studies file', type = 'character')
parser <- argparser::add_argument(parser, '--associations_db_file', help = 'Existing DuckDB associations file', type = 'character')

args <- argparser::parse_args(parser)


main <- function() {
  study_extractions <- vroom::vroom(file.path(args$results_dir, "all_study_blocks.tsv"), show_col_types = F)
  raw_coloc_results <- vroom::vroom(file.path(args$results_dir, "raw_coloc_results.tsv"), show_col_types = F)
  rare_results <- vroom::vroom(file.path(args$results_dir, "rare_results.tsv"), show_col_types = F)
  results_metadata <- vroom::vroom(file.path(args$results_dir, "results_metadata.tsv"), show_col_types = F)
  studies <- vroom::vroom(file.path(args$results_dir, "studies_processed.tsv"), show_col_types = F)
  variant_annotations_full <- vroom::vroom(file.path(variant_annotation_dir, "vep_annotations_hg38.tsv.gz"), show_col_types =  F) |>
   dplyr::mutate(id=1:dplyr::n())

  results_metadata <- results_metadata |>
    tidyr::separate(ld_block, into=c("pop", "chr", "start", "end"), sep="[/-]", remove=FALSE) |>
    dplyr::mutate(id=1:dplyr::n(), chr=as.integer(chr), start=as.integer(start), end=as.integer(end))
  
  vl <- parallel::mclapply(1:nrow(results_metadata), mc.cores=50, \(i) {
    chr <- results_metadata$chr[i]
    start <- results_metadata$start[i]
    end <- results_metadata$end[i]
    ind <- variant_annotations_full$CHR == chr & variant_annotations_full$BP >= start & variant_annotations_full$BP < end
    tibble::tibble(SNP=variant_annotations_full$SNP[ind], ld_block=results_metadata$ld_block[i])
  }) |> dplyr::bind_rows()

  variant_annotations_full <- dplyr::left_join(variant_annotations_full, vl, by="SNP")

  coloc <- raw_coloc_results |>
    dplyr::filter(posterior_prob > 0.5) |>
    dplyr::mutate(coloc_group_id=1:dplyr::n()) |>
    tidyr::separate_longer_delim(cols=traits, delim=", ") |>
    dplyr::rename(unique_study_id=traits)
  
  #this is needed to deduplicate study / snp pairs.
  coloc <- coloc |>
    dplyr::group_by(unique_study_id, candidate_snp) |>
    dplyr::slice_max(posterior_prob, n = 1, with_ties = FALSE) |>
    dplyr::ungroup()

  rare_results <- rare_results |>
    dplyr::mutate(id=1:dplyr::n()) |>
    tidyr::separate_longer_delim(cols=traits, delim=", ") |>
    dplyr::rename(unique_study_id=traits)

  coloc <- coloc |> dplyr::left_join(dplyr::select(study_extractions, study, unique_study_id, chr, bp, min_p, cis_trans, ld_block, known_gene), by="unique_study_id")

  # study_extractions includes variants that don't colocalise with anything. Need to get those, but they don't have variant IDs
  study_extractions <- study_extractions |> 
    dplyr::left_join(dplyr::select(variant_annotations_full, SNP, CHR, BP), by=c("chr"="CHR", "bp"="BP")) |>
    dplyr::filter(!duplicated(unique_study_id))
  dim(study_extractions)

  # Find finemapped results that colocalise with nothing
  no_coloc <- dplyr::filter(study_extractions, !unique_study_id %in% coloc$unique_study_id) |>
    dplyr::rename(candidate_snp=SNP, traits=unique_study_id) |>
    dplyr::mutate(id=(1:dplyr::n())+nrow(coloc))

  coloc <- dplyr::bind_rows(coloc, no_coloc)

  # Add ld_block to variant_annotations
  snp_and_ld_blocks <- coloc |>
    dplyr::filter(!duplicated(candidate_snp)) |>
    dplyr::select(candidate_snp, ld_block)

  variant_annotations_of_candidate_snps <- dplyr::left_join(variant_annotations_full, snp_and_ld_blocks, by=c("SNP"="candidate_snp"))

  populate_associations_db(ld_blocks, variant_annotations_of_candidate_snps$SNP)


  studies_db <- duckdb::dbConnect(duckdb::duckdb(), args$studies_db_file)
  


  tables_config <- list(
    "study_extractions" = list(
      data = study_extractions,
      keys = c("id"),
      foreign_keys = list(
        list(columns = c("study"), references = list(table = "studies", columns = c("id")))
      )
    ),
    "results_metadata" = list(
      data = results_metadata
    ),
    "studies" = list(
      data = studies,
      keys = c("id")
    ),
    "variant_annotations" = list(
      data = variant_annotations_full,
      keys = c("id"),
      unique = c("SNP")
    ),
    "colocalisations" = list(
      data = coloc,
      keys = c("study_extraction_id", "candidate_snp"),
      foreign_keys = list(
        list(columns = c("study_extraction_id"), references = list(table = "study_extractions", columns = c("id"))),
        list(columns = c("candidate_snp"), references = list(table = "variant_annotations", columns = c("SNP")))
      )
    ),
    "ld" = list(
      data = ldl,
      keys = c("lead", "variant", "ld_block")
    ),
    "rare_results" = list(
      data = rare_results,
      keys = c("id")
    ),
    # "assocs" = list(
    #   data = assocs,
    #   keys = c("SNP", "study")
    # )
  )


  # Process each table
  for (table_name in names(tables_config)) {
    append_unique_rows(
      studies_db, 
      table_name, 
      tables_config[[table_name]]$data,
      tables_config[[table_name]]$keys,
      tables_config[[table_name]]$unique
    )

    # Add foreign key constraints if defined
    if (!is.null(tables_config[[table_name]]$foreign_keys)) {
      for (fk in tables_config[[table_name]]$foreign_keys) {
        fk_cols <- paste(fk$columns, collapse = ", ")
        ref_cols <- paste(fk$references$columns, collapse = ", ")
        sql <- sprintf(
          "ALTER TABLE %s ADD CONSTRAINT fk_%s_%s FOREIGN KEY (%s) REFERENCES %s(%s)",
          table_name,
          table_name,
          fk$references$table,
          fk_cols,
          fk$references$table,
          ref_cols
        )
        tryCatch({
          duckdb::dbExecute(studies_db, sql)
        }, error = function(e) {
          message(sprintf("Could not add foreign key constraint to %s: %s", table_name, e$message))
        })
      }
    }
  }

  duckdb::dbDisconnect(studies_db, shutdown=TRUE)

  # Extract summary statistics
  varids <- unique(study_extractions$SNP)
  length(varids)
  varids_info <- tibble(varids) |> separate(varids, into=c("chr", "other"), sep=":", remove=FALSE)
  varids_list <- lapply(1:22, \(x) {
    subset(varids_info, chr==x)$varids
  })

  sources <- unique(studies$source)

  assocs <- lapply(sources, \(x) assocs_source(studies, varids_list, x, 30))
  associations_db <- duckdb::dbConnect(duckdb::duckdb(), args$associations_db_file)
  
  if (dbExistsTable(associations_db, "assocs")) {
    # Get existing combinations of SNP and study
    existing_keys <- duckdb::dbGetQuery(associations_db, "SELECT DISTINCT SNP, study FROM assocs")
    
    # Process each source's associations
    for (i in seq_along(assocs)) {
      if (nrow(assocs[[i]]) > 0) {
        # Filter out existing combinations
        new_rows <- dplyr::anti_join(assocs[[i]], existing_keys, by = c("SNP", "study"))
        if (nrow(new_rows) > 0) {
          duckdb::dbAppendTable(associations_db, "assocs", new_rows)
        }
      }
    }
  } else {
    duckdb::dbWriteTable(associations_db, "assocs", assocs[[1]])
    for(i in 2:length(assocs)) {
      if(nrow(assocs[[i]]) > 0) {
        duckdb::dbAppendTable(associations_db, "assocs", assocs[[i]])
      }
    }
  }
  
  duckdb::dbDisconnect(associations_db, shutdown=TRUE)

  ensure_dbs_are_valid()
}

populate_associations_db <- function(ld_blocks, snps) {
  associations_db <- duckdb::dbConnect(duckdb::duckdb(), args$associations_db_file)
  existing_keys <- duckdb::dbGetQuery(associations_db, "SELECT DISTINCT SNP FROM assocs")

  existing_keys <- as.data.frame(existing_keys)
  new_rows <- dplyr::anti_join(ld_blocks, existing_keys, by = "ld_block")
  
  if (nrow(new_rows) > 0) {
    return()
  }
  ldl <- parallel::mclapply(ld_blocks, \(x) generate_ld_obj(x, snps), mc.cores=30) |> bind_rows()

  duckdb::dbAppendTable(associations_db, "assocs", ldl)
}

generate_ld_obj <- function(ld_block, snps) {
  message(ld_block)
  ld <- suppressMessages(vroom::vroom(file.path(ld_reference_panel_dir, paste0(ld_block, ".unphased.vcor1")), show_col_types = F))
  ldvars <- suppressMessages(scan(file.path(ld_reference_panel_dir, paste0(ld_block, ".unphased.vcor1.vars")), character()))
  names(ld) <- ldvars
  ld$lead <- ldvars
  ind <- which(ldvars %in% snps)
  ld <- ld[ind,]
  ldl <- tidyr::pivot_longer(ld, cols=-lead, names_to="variant", values_to="r") |>
    dplyr::filter(r^2 > 0.8 | variant %in% ld$lead) |> 
    dplyr::filter(lead != variant) |>
    dplyr::mutate(ld_block=ld_block)
  return(ldl)
}

# Helper function to append new rows while avoiding duplicates
append_unique_rows <- function(conn, table_name, new_data, key_columns, unique_columns = NULL) {
  if (duckdb::dbExistsTable(conn, table_name)) {
    # Get existing keys
    key_cols_str <- paste(key_columns, collapse = ", ")
    existing_keys <- duckdb::dbGetQuery(conn, sprintf("SELECT DISTINCT %s FROM %s", key_cols_str, table_name))
    
    # Filter out existing rows
    existing_keys <- as.data.frame(existing_keys)
    new_rows <- dplyr::anti_join(new_data, existing_keys, by = key_columns)
    
    if (nrow(new_rows) > 0) {
      duckdb::dbAppendTable(conn, table_name, new_rows)
    }
    message(sprintf("Added %d new rows to %s", nrow(new_rows), table_name))
  } else {
    # Create table with primary key and unique constraints
    cols <- paste(names(new_data), sapply(new_data, function(x) typeof(x)), collapse=", ")
    key_cols_str <- paste(key_columns, collapse = ", ")
    
    # Add PRIMARY KEY and UNIQUE constraints
    constraints <- sprintf("PRIMARY KEY (%s)", key_cols_str)
    if (!is.null(unique_columns)) {
      unique_constraints <- sapply(unique_columns, function(col) {
        sprintf("UNIQUE (%s)", col)
      })
      constraints <- paste(c(constraints, unique_constraints), collapse=", ")
    }
    
    sql <- sprintf("CREATE TABLE %s (%s, %s)", table_name, cols, constraints)
    duckdb::dbExecute(conn, sql)
    
    # Insert the data
    duckdb::dbAppendTable(conn, table_name, new_data)
    message(sprintf("Created new table %s with %d rows", table_name, nrow(new_data)))
  }
}

extract_variants <- function(varids_list, path, study="study") {
  file_list <- list.files(path) |>
    file.path(path, .) |>
    grep("pre_filter", ., value=TRUE, invert=TRUE) |>
    grep("dentist", ., value=TRUE, invert=TRUE)

  if(length(file_list) == 0) {
    return(NULL)
  }

  ext <- lapply(file_list, \(y) {
    chr <- strsplit(basename(y), "_")[[1]][2] |> as.numeric()
    tryCatch({
      vroom::vroom(y, show_col_types = F) |>
        dplyr::filter(SNP %in% varids_list[[chr]]) |>
        dplyr::select(SNP, BETA, SE, IMPUTED, P, EAF) |>
        dplyr::mutate(study=study)
    }, error=function(e) {
      message(e)
      return(NULL)
    })
  }) 

  ext <- ext[!sapply(ext, is.null)] |> dplyr::bind_rows()
  return(ext)
}

assocs_source <- function(studies, varids_list, source, mc.cores=10) {
  assocs <- parallel::mclapply(studies$study_name[studies$source == source], \(x) {
    message(x)
    path <- file.path(data_dir, "/study/", x, "imputed")
    tryCatch({
      extract_variants(varids_list, path, x)
    }, error=function(e) {
      message(e)
      return(NULL)
    })
  }, mc.cores=mc.cores)
  assocs <- assocs[!sapply(assocs, \(x) inherits(x, "try-error"))]
  return(assocs |> dplyr::bind_rows())
}

ensure_dbs_are_valid <- function() {
  #TODO: Check that the databases are valid
}

main()