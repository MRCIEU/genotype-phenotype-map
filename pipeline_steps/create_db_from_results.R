source('constants.R')

parser <- argparser::arg_parser('Create DuckDB from pipeline results')
parser <- argparser::add_argument(parser, '--results_dir', help = 'Results directory of pipeline', type = 'character')
parser <- argparser::add_argument(parser, '--gpm_db_file', help = 'Existing DuckDB GPM file', type = 'character')

args <- argparser::parse_args(parser)

max_cores <- 50

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
  snp_annotations = list(
    name = "snp_annotations",
    query = "CREATE TABLE snp_annotations (
      id INTEGER PRIMARY KEY,
      snp TEXT,
      chr INTEGER,
      bp INTEGER,
      ea TEXT,
      oa TEXT,
      gene TEXT,
      feature_type TEXT,
      consequence TEXT,
      cdna_position TEXT,
      cds_position TEXT,
      protein_position TEXT,
      amino_acids TEXT,
      codons TEXT,
      rsid TEXT,
      impact TEXT,
      symbol TEXT,
      biotype TEXT,
      strand TEXT,
      canonical TEXT,
      all_af REAL CHECK (all_af BETWEEN 0 AND 1),
      eur_af REAL CHECK (eur_af BETWEEN 0 AND 1),
      amr_af REAL CHECK (amr_af BETWEEN 0 AND 1),
      eas_af REAL CHECK (eas_af BETWEEN 0 AND 1),
      sas_af REAL CHECK (sas_af BETWEEN 0 AND 1),
      afr_af REAL CHECK (afr_af BETWEEN 0 AND 1)
    )"
  ),
  study_extractions = list(
    name = "study_extractions",
    query = "CREATE TABLE study_extractions (
      id INTEGER PRIMARY KEY,
      study_id INTEGER,
      snp_id INTEGER,
      snp TEXT NOT NULL,
      ld_block_id INTEGER,
      unique_study_id TEXT NOT NULL,
      study TEXT NOT NULL,
      file TEXT NOT NULL,
      chr INTEGER NOT NULL,
      bp INTEGER NOT NULL,
      min_p REAL CHECK (min_p BETWEEN 0 AND 1),
      cis_trans TEXT,
      ld_block TEXT,
      known_gene TEXT,
      FOREIGN KEY (study_id) REFERENCES studies(id),
      FOREIGN KEY (snp_id) REFERENCES snp_annotations(id),
      FOREIGN KEY (ld_block_id) REFERENCES ld_blocks(id)
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
      posterior_explained_by_snp REAL CHECK (posterior_explained_by_snp BETWEEN 0 AND 1),
      candidate_snp TEXT,
      study_id INTEGER,
      chr INTEGER NOT NULL,
      bp INTEGER NOT NULL,
      min_p REAL CHECK (min_p BETWEEN 0 AND 1),
      cis_trans TEXT,
      ld_block TEXT,
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
      rare_result_group_id INTEGER,
      study_extraction_id INTEGER,
      snp_id INTEGER,
      ld_block_id INTEGER,
      unique_study_id TEXT,
      candidate_snp TEXT,
      study_id INTEGER,
      chr INTEGER,
      bp INTEGER,
      min_p REAL CHECK (min_p BETWEEN 0 AND 1),
      cis_trans TEXT,
      known_gene TEXT,
      FOREIGN KEY (study_extraction_id) REFERENCES study_extractions(id),
      FOREIGN KEY (snp_id) REFERENCES snp_annotations(id),
      FOREIGN KEY (ld_block_id) REFERENCES ld_blocks(id)
    )"
  ),
  ld = list(
    name = "ld",
    query = "CREATE TABLE ld (
      lead_snp_id INTEGER,
      variant_snp_id INTEGER,
      ld_block_id INTEGER,
      r REAL CHECK (r BETWEEN -1 AND 1),
      PRIMARY KEY (lead_snp_id, variant_snp_id),
      FOREIGN KEY (lead_snp_id) REFERENCES snp_annotations(id),
      FOREIGN KEY (variant_snp_id) REFERENCES snp_annotations(id),
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
  )
)

associations_table <- list(
  name = "assocs",
  query = "CREATE TABLE assocs (
    snp_id INTEGER,
    study_id INTEGER,
    beta REAL CHECK (beta BETWEEN -1 AND 1),
    se REAL CHECK (se > 0),
    imputed BOOLEAN,
    p REAL CHECK (p BETWEEN 0 AND 1),
    eaf REAL CHECK (eaf BETWEEN 0 AND 1),
    PRIMARY KEY (snp_id, study_id),
    FOREIGN KEY (snp_id) REFERENCES snp_annotations(id),
    FOREIGN KEY (study_id) REFERENCES studies(id)
  )"
)

gwas_upload_table <- list(
  name = "gwas_upload",
  query = "CREATE TABLE gwas_upload (
    guid TEXT,
    name TEXT,
    sample_size INTEGER,
    p_value_threshold REAL CHECK (p_value_threshold BETWEEN 0 AND 1),
    ancestry TEXT,
    category TEXT,
    is_published BOOLEAN,
    doi TEXT,
    should_be_added BOOLEAN,
    PRIMARY KEY (guid)
  )"
)

main <- function() {
  if (file.exists(args$gpm_db_file)) {
    file.copy(args$gpm_db_file, args$results_dir)
  } else {
    file.remove(file.path(args$results_dir, "gpm.db"))
  }
  new_gpm_db <- file.path(args$results_dir, "gpm.db")

  data <- load_data_for_simple_tables(simple_db_tables)
  simple_db_tables <- data$simple_db_tables
  all_relevant_snps <- data$all_relevant_snps
  
  gpm_db <- duckdb::dbConnect(duckdb::duckdb(), new_gpm_db)
  if (!file.exists(args$gpm_db_file)) {
    lapply(simple_db_tables, \(table) DBI::dbExecute(gpm_db, table$query))
    DBI::dbExecute(gpm_db, associations_table$query)
    DBI::dbExecute(gpm_db, gwas_upload_table$query)
  } 

  assocs <- load_data_for_associations(gpm_db, simple_db_tables, all_relevant_snps)
  DBI::dbAppendTable(gpm_db, associations_table$name, assocs)
  
  DBI::dbDisconnect(gpm_db, shutdown=TRUE)
  file.copy(new_gpm_db, args$gpm_db_file, overwrite=TRUE)
}

load_data_for_simple_tables <- function(simple_db_tables) {
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

  simple_db_tables$snp_annotations$data <- vroom::vroom(file.path(variant_annotation_dir, "vep_annotations_hg38.tsv.gz"), show_col_types =  F) |>
    dplyr::rename_with(tolower) |>
    dplyr::mutate(id=1:dplyr::n(), snp=trimws(snp))

  simple_db_tables$study_extractions$data <- simple_db_tables$study_extractions$data |>
    dplyr::mutate(id=1:dplyr::n(), file=sub(ld_block_data_dir, "", file)) |>
    dplyr::left_join(dplyr::select(simple_db_tables$studies$data, study_name, id) |> dplyr::rename(study_id=id), by=c("study"="study_name")) |>
    dplyr::left_join(dplyr::select(simple_db_tables$ld_blocks$data, ld_block, id) |> dplyr::rename(ld_block_id=id), by="ld_block")
  
  snp_annotation_subset <- simple_db_tables$snp_annotations$data |>
    dplyr::select(snp, id, chr, bp) |>
    dplyr::rename(snp_id=id)

  simple_db_tables$study_extractions$data <- simple_db_tables$study_extractions$data |>
    dplyr::left_join(snp_annotation_subset, by=c("chr"="chr", "bp"="bp")) |>
    dplyr::filter(!duplicated(unique_study_id))
  #TODO: Remove this once we have fixed the SNP annotations
  problematic_extractions <- simple_db_tables$study_extractions$data |>
    dplyr::filter(is.na(snp_id))
  vroom::vroom_write(problematic_extractions, "/home/wt23152/problematic_extractions.tsv", show_col_types = F)
  simple_db_tables$study_extractions$data <- simple_db_tables$study_extractions$data |>
    dplyr::filter(!is.na(snp_id))

  simple_db_tables$colocalisations$data <- vroom::vroom(file.path(args$results_dir, "raw_coloc_results.tsv"), show_col_types = F) |> 
    format_colocalisations(simple_db_tables$study_extractions$data, simple_db_tables$snp_annotations$data)
  
  #TODO: Remove this once we have fixed SNP annotations
  problematic_colocalisations <- simple_db_tables$colocalisations$data |>
    dplyr::filter(is.na(study_extraction_id))
  vroom::vroom_write(problematic_colocalisations, "/home/wt23152/problematic_colocalisations.tsv", show_col_types = F)
  simple_db_tables$colocalisations$data <- simple_db_tables$colocalisations$data |>
    dplyr::filter(!is.na(study_extraction_id))

  simple_db_tables$rare_results$data <- vroom::vroom(file.path(args$results_dir, "rare_results.tsv"), show_col_types = F) |>
    dplyr::rename_with(tolower) |>
    format_rare_results(simple_db_tables$study_extractions$data, simple_db_tables$snp_annotations$data)
  #TODO: Remove this once we have fixed SNP annotations
  problematic_rare_results <- simple_db_tables$rare_results$data |>
    dplyr::filter(is.na(study_extraction_id))
  vroom::vroom_write(problematic_rare_results, "/home/wt23152/problematic_rare_results.tsv", show_col_types = F)
  simple_db_tables$rare_results$data <- simple_db_tables$rare_results$data |>
    dplyr::filter(!is.na(study_extraction_id))

  simple_db_tables$results_metadata$data <- vroom::vroom(file.path(args$results_dir, "results_metadata.tsv"), show_col_types = F) |>
    dplyr::left_join(dplyr::select(simple_db_tables$ld_blocks$data, ld_block, id) |> dplyr::rename(ld_block_id=id), by="ld_block")


  ld_data <- retrieve_ld_data(simple_db_tables$ld_blocks$data,
                              simple_db_tables$snp_annotations$data,
                              simple_db_tables$study_extractions$data,
                              simple_db_tables$colocalisations$data) 

  simple_db_tables$ld$data <- ld_data$ld_data
  all_relevant_snps <- ld_data$all_relevant_snps

  return(list(simple_db_tables = simple_db_tables, all_relevant_snps = all_relevant_snps))
}

format_colocalisations <- function(colocalisations, study_extractions, snp_annotations) {
  study_extractions_subset <- study_extractions |>
    dplyr::select(id, study_id, unique_study_id, chr, bp, min_p, cis_trans, ld_block, ld_block_id, known_gene) |>
    dplyr::rename(study_extraction_id=id)

  snp_annotations_subset <- snp_annotations |>
    dplyr::select(snp, id) |>
    dplyr::rename(snp_id=id)

  colocalisations <- colocalisations |>
    dplyr::filter(posterior_prob > 0.5) |>
    dplyr::mutate(coloc_group_id=1:dplyr::n(), candidate_snp=trimws(candidate_snp)) |>
    dplyr::left_join(snp_annotations_subset, by=c("candidate_snp"="snp"), relationship="many-to-one") |>
    tidyr::separate_longer_delim(cols=traits, delim=", ") |>
    dplyr::rename(unique_study_id=traits) |>
    dplyr::left_join(study_extractions_subset, by="unique_study_id", relationship="many-to-one") |>
    dplyr::group_by(study_extraction_id, snp_id) |>
    dplyr::slice_max(posterior_prob, n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::select(-dropped_trait)

  return(colocalisations)
}

format_rare_results <- function(rare_results, study_extractions, snp_annotations) {
  study_extractions_subset <- study_extractions |>
    dplyr::select(id, study_id, unique_study_id, chr, bp, min_p, cis_trans, ld_block_id, known_gene) |>
    dplyr::rename(study_extraction_id=id)

  snp_annotations_subset <- snp_annotations |>
    dplyr::select(snp, id) |>
    dplyr::rename(snp_id=id)

  rare_results <- rare_results |>
    dplyr::mutate(rare_result_group_id=1:dplyr::n()) |>
    tidyr::separate_longer_delim(cols=traits, delim=", ") |>
    dplyr::rename(unique_study_id=traits) |>
    dplyr::left_join(study_extractions_subset, by="unique_study_id") |>
    dplyr::left_join(snp_annotations_subset, by=c("candidate_snp"="snp"), relationship="many-to-many")

  return(rare_results)
}

retrieve_ld_data <- function(ld_blocks, variant_annotations, study_extractions, colocalisations) {
  colocalising_snps <- dplyr::select(colocalisations, candidate_snp, snp_id, ld_block) |> dplyr::distinct()

  # Find finemapped results that colocalise with nothing
  non_colocalising_snps <- study_extractions |>
    dplyr::filter(!unique_study_id %in% colocalisations$unique_study_id) |>
    dplyr::select(snp, snp_id, ld_block) |>
    dplyr::rename(candidate_snp=snp) |>
    dplyr::distinct()

  all_relevant_snps <- dplyr::bind_rows(colocalising_snps, non_colocalising_snps) |> dplyr::distinct()

  relevant_ld_blocks <- unique(study_extractions$ld_block)
  ld_data <- parallel::mclapply(relevant_ld_blocks, mc.cores=max_cores, \(ld_block) {
    relevant_snps <- all_relevant_snps |> dplyr::filter(ld_block == ld_block)
    generate_ld_obj(ld_block, relevant_snps$candidate_snp)
  }) |> dplyr::bind_rows()

  # swap lead and variant snps and ld_block for their respective ids 
  variant_annotations_subset <- dplyr::select(variant_annotations, snp, id)
  ld_blocks_subset <- dplyr::select(ld_blocks, ld_block, id) |> dplyr::rename(ld_block_id=id)
  ld_data <- ld_data |>
    dplyr::left_join(variant_annotations_subset |> dplyr::rename(lead_snp_id=id), by=c("lead"="snp")) |>
    dplyr::left_join(variant_annotations_subset |> dplyr::rename(variant_snp_id=id), by=c("variant"="snp")) |>
    dplyr::left_join(ld_blocks_subset, by="ld_block") |>
    dplyr::select(-lead, -variant, -ld_block)

  return(list(ld_data = ld_data, all_relevant_snps = all_relevant_snps))
}

generate_ld_obj <- function(ld_block, snps) {
  ld_file <- file.path(ld_reference_panel_dir, glue::glue("{ld_block}.unphased.vcor1"))
  ld <- vroom::vroom(ld_file, col_names = FALSE, show_col_types = F)
  ldvars <- vroom::vroom(glue::glue("{ld_file}.vars")), col_names = FALSE, show_col_types = F)

  names(ld) <- ldvars$X1
  ld$lead <- ldvars$X1
  ind <- which(ldvars$X1 %in% snps)
  ld <- ld[ind,]

  ldl <- tidyr::pivot_longer(ld, cols=-lead, names_to="variant", values_to="r") |>
    dplyr::filter(r^2 > 0.8 | variant %in% ld$lead) |> 
    dplyr::filter(lead != variant) |>
    dplyr::mutate(ld_block=ld_block)
  return(ldl)
}

append_unique_rows <- function(conn, table) {
  if (is.null(table$data) || nrow(table$data) == 0) return()
  # If the id column doesn't exist, we can replace the table, otherwise we need to filter out existing rows
  id_column_exists <- DBI::dbGetQuery(conn, glue::glue("SELECT COUNT(*) FROM pragma_table_info('{table$name}') WHERE name = 'id';"))

  if (id_column_exists$count_star== 0) {
    DBI::dbExecute(conn, glue::glue("TRUNCATE TABLE {table$name}"))
  } else {
    existing_ids <- DBI::dbGetQuery(conn, glue::glue("SELECT DISTINCT id FROM {table$name}"))
    table$data <- dplyr::anti_join(table$data, existing_ids, by = "id")
  }

  if (nrow(table$data) > 0) {
    DBI::dbAppendTable(conn, table$name, table$data)
  }
  message(glue::glue("Added {nrow(table$data)} new rows to {table$name}"))
}

load_data_for_associations <- function(gpm_db, simple_db_tables, all_relevant_snps) {
  existing_assocs <- DBI::dbGetQuery(gpm_db, "SELECT study_id, snp_id FROM assocs")
  message(nrow(existing_assocs), ' existing associations found for ', length(unique(existing_assocs$study_id)), ' studies')

  all_study_extractions_by_study <- split(simple_db_tables$study_extractions$data, simple_db_tables$study_extractions$data$study_id)

  new_assocs <- purrr::imap(all_study_extractions_by_study, \(study_extractions, study_id) {
    relevant_ld_blocks <- study_extractions |>
      dplyr::select(ld_block) |>
      dplyr::distinct()

    existing_assocs <- existing_assocs |>
      dplyr::filter(study_id == study_id)

    relevant_snps_subset <- all_relevant_snps |>
      dplyr::filter(ld_block %in% relevant_ld_blocks) |>
      dplyr::filter(!snp_id %in% existing_assocs$snp_id)

    if (nrow(relevant_snps_subset) == 0) return(data.frame())
    return(data.frame(study_name=study_extractions$study[1],
                      candidate_snp=relevant_snps_subset$candidate_snp,
                      ld_block=relevant_snps_subset$ld_block
    ))
  }) |> dplyr::bind_rows()

  message('Retrieving ', nrow(new_assocs), ' new associations for ', length(unique(new_assocs$study_name)), ' studies')

  assocs_by_ld_block <- split(new_assocs, new_assocs$ld_block)
  assocs <- parallel::mclapply(assocs_by_ld_block, mc.cores=max_cores, \(assocs) get_associations_for_ld_block(assocs)) |>
    dplyr::bind_rows()
  
  snp_annotations_subset <- simple_db_tables$snp_annotations$data |>
    dplyr::select(snp, id) |>
    dplyr::rename(snp_id=id)

  studies_subset <- simple_db_tables$studies$data |>
    dplyr::select(study_name, id) |>
    dplyr::rename(study_id=id)

  assocs <- assocs |> 
    dplyr::left_join(snp_annotations_subset, by="snp") |>
    dplyr::left_join(studies_subset, by=c("study"="study_name")) |>
    dplyr::select(snp_id, study_id, beta, se, imputed, p, eaf)

  return(assocs)
}

get_associations_for_ld_block <- function(assocs) {
  imputed_studies <- vroom::vroom(file.path(ld_block_data_dir, assocs$ld_block[1], "imputed_studies.tsv"), show_col_types = F) |>
    dplyr::filter(study %in% assocs$study)

  associations <- apply(imputed_studies, 1, \(study) {
    vroom::vroom(study[['file']], show_col_types = F) |>
      dplyr::filter(SNP %in% assocs$candidate_snp) |>
      dplyr::select(SNP, BETA, SE, IMPUTED, P, EAF) |>
      dplyr::mutate(study=study[['study']]) |>
      dplyr::rename_with(tolower)
  })

  associations <- associations[!sapply(associations, is.null)] |> dplyr::bind_rows()
  return(associations)
}


main()