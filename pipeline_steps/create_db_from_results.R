source('constants.R')

parser <- argparser::arg_parser('Create DuckDB from pipeline results')
parser <- argparser::add_argument(parser, '--results_dir', help = 'Results directory of pipeline', type = 'character')
parser <- argparser::add_argument(parser, '--gpm_db_file', help = 'Existing DuckDB GPM file', type = 'character')

args <- argparser::parse_args(parser)

#args <- list(results_dir='/local-scratch/projects/genotype-phenotype-map/results/2025_02_19-09_08', gpm_db_file='/local-scratch/projects/genotype-phenotype-map/results/gpm_db.db')
simple_db_tables <- list(
  study_sources = list(
    name = "study_sources",
    query = "CREATE TABLE study_sources (
      id INTEGER PRIMARY KEY,
      source TEXT NOT NULL,
      name TEXT NOT NULL,
      url TEXT NOT NULL,
      doi TEXT NOT NULL
    )"
  ),
  ld_blocks = list(
    name = "ld_blocks",
    query = "CREATE TABLE ld_blocks (
      id INTEGER PRIMARY KEY,
      chr INTEGER NOT NULL,
      start INTEGER NOT NULL,
      stop INTEGER NOT NULL,
      ancestry TEXT NOT NULL,
      ld_block TEXT NOT NULL
    )"
  ),
  studies = list(
    name = "studies",
    query = "CREATE TABLE studies (
      id INTEGER PRIMARY KEY,
      data_type TEXT NOT NULL,
      data_format TEXT NOT NULL,
      study_name TEXT NOT NULL,
      trait TEXT NOT NULL,
      ancestry TEXT NOT NULL,
      sample_size INTEGER NOT NULL,
      category TEXT,
      study_location TEXT,
      extracted_location TEXT,
      probe TEXT,
      tissue TEXT,
      source_id INTEGER ,
      variant_type TEXT,
      p_value_threshold REAL,
      gene TEXT,
      FOREIGN KEY (source_id) REFERENCES study_sources(id)
    )"
  ),
  study_extractions = list(
    name = "study_extractions",
    query = "CREATE TABLE study_extractions (
      id INTEGER PRIMARY KEY,
      study_id INTEGER,
      ld_block_id INTEGER,
      unique_study_id TEXT,
      study TEXT,
      file TEXT,
      chr INTEGER,
      bp INTEGER,
      min_p REAL CHECK (min_p BETWEEN 0 AND 1),
      cis_trans TEXT,
      ld_block TEXT,
      known_gene TEXT,
      FOREIGN KEY (study_id) REFERENCES studies(id),
      FOREIGN KEY (ld_block_id) REFERENCES ld_blocks(id)
    )"
  ),
  results_metadata = list(
    name = "results_metadata",
    query = "CREATE TABLE results_metadata (
      ld_block_id INTEGER,
      ld_block TEXT,
      number_extracted INTEGER,
      number_standardised INTEGER,
      mean_snps_removed_by_reference_panel REAL,
      number_imputed INTEGER,
      significant_snps_imputed INTEGER,
      significant_imputed_snps_filtered INTEGER,
      number_finemapped INTEGER,
      finemapped_per_imputed REAL,
      num_finemap_failed INTEGER,
      standardised_time_taken REAL,
      imputed_time_taken REAL,
      finemapped_time_taken REAL,
      FOREIGN KEY (ld_block_id) REFERENCES ld_blocks(id)
    )"
  ),
  snp_annotations = list(
    name = "snp_annotations",
    query = "CREATE TABLE snp_annotations (
      id INTEGER PRIMARY KEY,
      SNP TEXT,
      CHR INTEGER,
      BP INTEGER,
      EA TEXT,
      OA TEXT,
      EAF REAL,
      Feature_type TEXT,
      Consequence TEXT,
      cDNA_position TEXT,
      CDS_position TEXT,
      Protein_position TEXT,
      Amino_acids TEXT,
      Codons TEXT,
      RSID TEXT,
      impact TEXT,
      symbol TEXT,
      biotype TEXT,
      strand TEXT,
      canonical TEXT
    )"
  ),
  colocalisations = list(
    name = "colocalisations",
    query = "CREATE TABLE colocalisations (
      study_extraction_id INTEGER,
      snp_id INTEGER,
      ld_block_id INTEGER,
      coloc_group_id INTEGER,
      iteration INTEGER,
      unique_study_id TEXT,
      posterior_prob REAL CHECK (posterior_prob BETWEEN 0 AND 1),
      regional_prob REAL CHECK (regional_prob BETWEEN 0 AND 1),
      posterior_prob_explained_by_snp REAL CHECK (posterior_prob_explained_by_snp BETWEEN 0 AND 1),
      candidate_snp TEXT,
      study_id INTEGER,
      chr INTEGER,
      bp INTEGER,
      min_p REAL CHECK (min_p BETWEEN 0 AND 1),
      cis_trans TEXT,
      known_gene TEXT,
      PRIMARY KEY (study_extraction_id, snp_id),
      FOREIGN KEY (study_extraction_id) REFERENCES study_extractions(id),
      FOREIGN KEY (snp_id) REFERENCES snp_annotations(id),
      FOREIGN KEY (ld_block_id) REFERENCES ld_blocks(id)
    )"
  ),
  rare_results = list(
    name = "rare_results",
    query = "CREATE TABLE rare_results (
      study_extraction_id INTEGER,
      snp_id INTEGER,
      unique_study_id TEXT,
      candidate_snp TEXT,
      rare_result_group_id INTEGER,
      study_id INTEGER,
      chr INTEGER,
      bp INTEGER,
      min_p REAL CHECK (min_p BETWEEN 0 AND 1),
      cis_trans TEXT,
      ld_block_id INTEGER,
      known_gene TEXT,
      FOREIGN KEY (study_extraction_id) REFERENCES study_extractions(id),
      FOREIGN KEY (snp_id) REFERENCES snp_annotations(id),
      FOREIGN KEY (ld_block_id) REFERENCES ld_blocks(id)
    )"
  )
)

complex_db_tables <- list(
  ld = list(
    name = "ld",
    query = "CREATE TABLE ld (
      lead_snp_id INTEGER,
      variant_snp_id INTEGER,
      ld_block_id INTEGER,
      r2 REAL CHECK (r2 BETWEEN 0 AND 1),
      PRIMARY KEY (lead_snp_id, variant_snp_id),
      FOREIGN KEY (lead_snp_id) REFERENCES snp_annotations(id),
      FOREIGN KEY (variant_snp_id) REFERENCES snp_annotations(id),
      FOREIGN KEY (ld_block_id) REFERENCES ld_blocks(id)
    )"
  ),
  assocs = list(
    name = "assocs",
    query = "CREATE TABLE assocs (
      snp_id INTEGER,
      study_id INTEGER,
      BETA REAL CHECK (BETA BETWEEN -1 AND 1),
      SE REAL CHECK (SE > 0),
      IMPUTED BOOLEAN,
      P REAL CHECK (P BETWEEN 0 AND 1),
      EAF REAL CHECK (EAF BETWEEN 0 AND 1),
      PRIMARY KEY (snp_id, study_id),
      FOREIGN KEY (snp_id) REFERENCES snp_annotations(id),
      FOREIGN KEY (study_id) REFERENCES studies(id)
    )"
  )
)

main <- function() {
  if (file.exists(args$gpm_db_file)) {
    file.copy(args$gpm_db_file, args$results_dir)
  }
  new_gpm_db <- file.path(args$results_dir, "gpm_db.db")

  simple_db_tables$study_sources$data <- vroom::vroom(file.path("data/study_sources.csv"), show_col_types = F) |>
    dplyr::mutate(id=1:dplyr::n())
  simple_db_tables$ld_blocks$data <- vroom::vroom(file.path("data/ld_blocks.tsv"), show_col_types = F) |>
    dplyr::mutate(id=1:dplyr::n(), ld_block=paste0(ancestry, "/", chr, "/", start, "-", stop))

  simple_db_tables$studies$data <- vroom::vroom(file.path(args$results_dir, "studies_processed.tsv"), show_col_types = F) |>
    dplyr::left_join(dplyr::select(simple_db_tables$study_sources$data, name, id) |> dplyr::rename(source_id=id), by=c("source"="name"))
  simple_db_tables$study_extractions$data <- vroom::vroom(file.path(args$results_dir, "all_study_blocks.tsv"), show_col_types = F)

  # Remove the studies that don't have any study extractions
  simple_db_tables$studies$data <- simple_db_tables$studies$data |>
    dplyr::filter(study_name %in% simple_db_tables$study_extractions$data$study) |>
    dplyr::mutate(id=1:dplyr::n()) |>
    dplyr::select(-reference_build, -source)

  simple_db_tables$study_extractions$data <- simple_db_tables$study_extractions$data |>
    dplyr::mutate(id=1:dplyr::n()) |>
    dplyr::left_join(dplyr::select(simple_db_tables$studies$data, study_name, id) |> dplyr::rename(study_id=id), by=c("study"="study_name")) |>
    dplyr::left_join(dplyr::select(simple_db_tables$ld_blocks$data, ld_block, id) |> dplyr::rename(ld_block_id=id), by="ld_block")

  gpm_db <- duckdb::dbConnect(duckdb::duckdb(), new_gpm_db)
  if (!file.exists(args$gpm_db_file)) {
    lapply(simple_db_tables, \(table) DBI::dbExecute(gpm_db, table$query))
  } 

  lapply(simple_db_tables, \(table) append_unique_rows(gpm_db, table))
  q()

  simple_db_tables$snp_annotations$data <- vroom::vroom(file.path(variant_annotation_dir, "vep_annotations_hg38.tsv.gz"), show_col_types =  F) |>
    dplyr::mutate(id=1:dplyr::n())

  colocalisations <- vroom::vroom(file.path(args$results_dir, "raw_coloc_results.tsv"), show_col_types = F) |> 
    format_colocalisations(simple_db_tables$study_extractions$data, simple_db_tables$snp_annotations$data)

  rare_results <- vroom::vroom(file.path(args$results_dir, "rare_results.tsv"), show_col_types = F) |>
    format_rare_results(simple_db_tables$study_extractions$data, simple_db_tables$snp_annotations$data)

  results_metadata <- vroom::vroom(file.path(args$results_dir, "results_metadata.tsv"), show_col_types = F) |>
    dplyr::left_join(dplyr::select(simple_db_tables$ld_blocks$data, ld_block, id) |> dplyr::rename(ld_block_id=id), by="ld_block")


  # study_extractions includes variants that don't colocalise with anything. Need to get those, but they don't have variant IDs
  study_extractions <- study_extractions |> 
    dplyr::left_join(dplyr::select(snp_annotations, SNP, CHR, BP), by=c("chr"="CHR", "bp"="BP")) |>
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

  variant_annotations_of_candidate_snps <- dplyr::left_join(snp_annotations, snp_and_ld_blocks, by=c("SNP"="candidate_snp"))

  populate_associations_db(ld_blocks, variant_annotations_of_candidate_snps$SNP)


  studies_db <- duckdb::dbConnect(duckdb::duckdb(), args$studies_db_file)
  
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
  
  if (DBI::dbExistsTable(associations_db, "assocs")) {
    # Get existing combinations of SNP and study
    existing_keys <- DBI::dbGetQuery(associations_db, "SELECT DISTINCT SNP, study FROM assocs")
    
    # Process each source's associations
    for (i in seq_along(assocs)) {
      if (nrow(assocs[[i]]) > 0) {
        # Filter out existing combinations
        new_rows <- dplyr::anti_join(assocs[[i]], existing_keys, by = c("SNP", "study"))
        if (nrow(new_rows) > 0) {
          DBI::dbAppendTable(associations_db, "assocs", new_rows)
        }
      }
    }
  } else {
    DBI::dbWriteTable(associations_db, "assocs", assocs[[1]])
    for(i in 2:length(assocs)) {
      if(nrow(assocs[[i]]) > 0) {
        DBI::dbAppendTable(associations_db, "assocs", assocs[[i]])
      }
    }
  }
  
  DBI::dbDisconnect(associations_db, shutdown=TRUE)

  ensure_dbs_are_valid()
}

format_colocalisations <- function(colocalisations, study_extractions, snp_annotations) {
  study_extractions_subset <- study_extractions |>
    dplyr::select(id, study_id, unique_study_id, chr, bp, min_p, cis_trans, ld_block_id, known_gene) |>
    dplyr::rename(study_extraction_id=id)

  snp_annotations_subset <- snp_annotations |>
    dplyr::select(SNP, id) |>
    dplyr::rename(snp_id=id)

  colocalisations <- colocalisations |>
    dplyr::filter(posterior_prob > 0.5) |>
    dplyr::mutate(coloc_group_id=1:dplyr::n()) |>
    tidyr::separate_longer_delim(cols=traits, delim=", ") |>
    dplyr::rename(unique_study_id=traits) |>
    dplyr::left_join(study_extractions_subset, by="unique_study_id") |>
    dplyr::left_join(snp_annotations_subset, by=c("candidate_snp"="SNP"), relationship="many-to-many") |>
    dplyr::group_by(study_extraction_id, snp_id) |>
    dplyr::slice_max(posterior_prob, n = 1, with_ties = FALSE) |>
    dplyr::ungroup()

  return(colocalisations)
}

format_rare_results <- function(rare_results, study_extractions, snp_annotations) {
  study_extractions_subset <- study_extractions |>
    dplyr::select(id, study_id, unique_study_id, chr, bp, min_p, cis_trans, ld_block_id, known_gene) |>
    dplyr::rename(study_extraction_id=id)

  snp_annotations_subset <- snp_annotations |>
    dplyr::select(SNP, id) |>
    dplyr::rename(snp_id=id)

  rare_results <- rare_results |>
    dplyr::mutate(rare_result_group_id=1:dplyr::n()) |>
    tidyr::separate_longer_delim(cols=traits, delim=", ") |>
    dplyr::rename(unique_study_id=traits) |>
    dplyr::left_join(study_extractions_subset, by="unique_study_id") |>
    dplyr::left_join(snp_annotations_subset, by=c("candidate_snp"="SNP"), relationship="many-to-many")

  return(rare_results)
}


populate_associations_db <- function(ld_blocks, snps) {
  associations_db <- duckdb::dbConnect(duckdb::duckdb(), args$associations_db_file)
  existing_keys <- DBI::dbGetQuery(associations_db, "SELECT DISTINCT SNP FROM assocs")

  existing_keys <- as.data.frame(existing_keys)
  new_rows <- dplyr::anti_join(ld_blocks, existing_keys, by = "ld_block")
  
  if (nrow(new_rows) > 0) {
    return()
  }

  DBI::dbAppendTable(associations_db, "assocs", ldl)
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

create_table <- function(db, query) {
  DBI::dbExecute(db, query)
}

append_unique_rows <- function(conn, table) {
  if (is.null(table$data) || nrow(table$data) == 0) return()
  # If the id column doesn't exist, we can replace the table, otherwise we need to filter out existing rows
  id_column_exists <- DBI::dbGetQuery(conn, sprintf("SELECT COUNT(*) FROM pragma_table_info('%s') WHERE name = 'id';", table$name))

  if (id_column_exists$count_star== 0) {
    DBI::dbExecute(conn, glue::glue("TRUNCATE TABLE {table$name}"))
  } else {
    existing_ids <- DBI::dbGetQuery(conn, glue::glue("SELECT DISTINCT id FROM {table$name}"))
    table$data <- dplyr::anti_join(table$data, existing_ids, by = "id")
  }

  if (nrow(table$data) > 0) {
    DBI::dbAppendTable(conn, table$name, table$data)
  }
  message(sprintf("Added %d new rows to %s", nrow(table$data), table$name))
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