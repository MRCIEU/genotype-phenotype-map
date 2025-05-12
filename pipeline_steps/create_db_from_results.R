source('constants.R')
source('database_definitions.R')

parser <- argparser::arg_parser('Create DuckDB from pipeline results')
parser <- argparser::add_argument(parser, '--results_dir', help = 'Results directory of pipeline', type = 'character')
parser <- argparser::add_argument(parser, '--studies_db_file', help = 'Studies database file', type = 'character')
parser <- argparser::add_argument(parser, '--associations_db_file', help = 'Associations database file', type = 'character')
parser <- argparser::add_argument(parser, '--ld_db_file', help = 'LD database file', type = 'character')
parser <- argparser::add_argument(parser, '--gwas_upload_db_file', help = 'GWAS upload database file', type = 'character')

args <- argparser::parse_args(parser)

max_cores <- 40

main <- function() {
  if (file.exists(args$studies_db_file)) file.remove(args$studies_db_file)
  if (file.exists(args$associations_db_file)) file.remove(args$associations_db_file)
  if (file.exists(args$ld_db_file)) file.remove(args$ld_db_file)
  if (file.exists(args$gwas_upload_db_file)) file.remove(args$gwas_upload_db_file)

  studies_conn <- duckdb::dbConnect(duckdb::duckdb(), args$studies_db_file)
  associations_conn <- duckdb::dbConnect(duckdb::duckdb(), args$associations_db_file)
  ld_conn <- duckdb::dbConnect(duckdb::duckdb(), args$ld_db_file)
  gwas_upload_conn <- duckdb::dbConnect(duckdb::duckdb(), args$gwas_upload_db_file)

  lapply(studies_db, \(table) DBI::dbExecute(studies_conn, table$query))
  lapply(gwas_upload_db, \(table) DBI::dbExecute(gwas_upload_conn, table$query))
  DBI::dbExecute(associations_conn, associations_table$query)
  DBI::dbExecute(ld_conn, ld_table$query)

  studies_db <- load_data_for_studies_db(studies_db)

  lapply(studies_db, \(table) append_unique_rows(studies_conn, table))

  all_relevant_snps <- find_relevant_snps(studies_db)
  message("Found ", nrow(all_relevant_snps), " relevant SNPs")

  load_data_into_ld_db(ld_conn, studies_db, all_relevant_snps)
  load_data_into_associations_db(associations_conn, studies_db, all_relevant_snps)

  DBI::dbDisconnect(studies_conn, shutdown=TRUE)
  DBI::dbDisconnect(associations_conn, shutdown=TRUE)
  DBI::dbDisconnect(ld_conn, shutdown=TRUE)
  DBI::dbDisconnect(gwas_upload_conn, shutdown=TRUE)

  file.copy(args$studies_db_file, file.path(results_dir, "latest/studies.db"))
  file.copy(args$associations_db_file, file.path(results_dir, "latest/associations.db"))
  file.copy(args$ld_db_file, file.path(results_dir, "latest/ld.db"))
  file.copy(args$gwas_upload_db_file, file.path(results_dir, "latest/gwas_upload.db"))
}

load_data_for_studies_db <- function(studies_db) {
  studies_db$study_sources$data <- vroom::vroom(file.path("data/study_sources.csv"), show_col_types = F) |>
    dplyr::mutate(id=1:dplyr::n())
  studies_db$ld_blocks$data <- vroom::vroom(file.path("data/ld_blocks.tsv"), show_col_types = F) |>
    dplyr::mutate(id=1:dplyr::n(), ld_block=paste0(ancestry, "/", chr, "/", start, "-", stop))

  studies_db$study_extractions$data <- vroom::vroom(file.path(args$results_dir, "study_extractions.tsv"), show_col_types = F)

  sources_subset <- studies_db$study_sources$data |>
    dplyr::select(source, id) |>
    dplyr::rename(source_id=id)

  # Remove the studies that don't have any study extractions
  studies_db$studies$data <- vroom::vroom(file.path(args$results_dir, "studies_processed.tsv.gz"), show_col_types = F) |>
    dplyr::left_join(sources_subset, by=c("source"="source")) |>
    dplyr::filter(study_name %in% studies_db$study_extractions$data$study & trait %in% studies_db$study_extractions$data$study) |>
    dplyr::mutate(id=1:dplyr::n()) |>
    dplyr::select(-reference_build, -source, -trait_name)

  # Remove the traits that don't have any study extractions
  studies_db$traits$data <- vroom::vroom(file.path(args$results_dir, "traits_processed.tsv.gz"), show_col_types = F) |>
    dplyr::filter(study_name %in% studies_db$studies$data$study_name) |>
    dplyr::mutate(id=1:dplyr::n()) |>
    dplyr::rename(trait_name=trait, trait=study_name)

  traits_subset <- studies_db$traits$data |>
    dplyr::select(trait, id) |>
    dplyr::rename(trait_id=id)
  
  studies_db$studies$data <- studies_db$studies$data |>
    dplyr::left_join(traits_subset, by=c("trait"="trait")) |>
    dplyr::select(-trait)

  studies_db$snp_annotations$data <- vroom::vroom(file.path(variant_annotation_dir, "vep_annotations_hg38.tsv.gz"), show_col_types =  F) |>
    dplyr::rename_with(tolower) |>
    dplyr::mutate(id=1:dplyr::n(), snp=trimws(snp))

  studies_db$gene_annotations$data <- vroom::vroom(file.path(liftover_dir, "gene_name_map.tsv"), show_col_types = F) |>
    dplyr::mutate(id=1:dplyr::n()) |>
    dplyr::rename_with(tolower) |>
    dplyr::rename(start=bp_start, end=bp_end, symbol=gene_name)

  studies_db$study_extractions$data <- studies_db$study_extractions$data |>
    dplyr::mutate(id=1:dplyr::n(), file=sub(ld_block_data_dir, "", file)) |>
    dplyr::left_join(dplyr::select(studies_db$studies$data, study_name, id) |> dplyr::rename(study_id=id), by=c("study"="study_name")) |>
    dplyr::left_join(dplyr::select(studies_db$ld_blocks$data, ld_block, id) |> dplyr::rename(ld_block_id=id), by="ld_block")
  
  snp_annotation_subset <- studies_db$snp_annotations$data |>
    dplyr::select(snp, id, chr, bp) |>
    dplyr::rename(snp_id=id)

  # TODO: Fix this is.na(snp_id) filter once it's fixed
  studies_db$study_extractions$data <- studies_db$study_extractions$data |>
    dplyr::left_join(snp_annotation_subset, by=c("chr"="chr", "bp"="bp")) |>
    dplyr::filter(!duplicated(unique_study_id)) |>
    dplyr::filter(!is.na(snp_id))

  # TODO: Fix this is.na(study_extraction_id) filter once it's fixed
  studies_db$colocalisations$data <- vroom::vroom(file.path(args$results_dir, "raw_coloc_results.tsv"), show_col_types = F) |> 
    format_colocalisations(studies_db$study_extractions$data, studies_db$snp_annotations$data) |>
    dplyr::filter(!is.na(study_extraction_id))

  # TODO: Fix this is.na(study_extraction_id) filter once it's fixed
  studies_db$rare_results$data <- vroom::vroom(file.path(args$results_dir, "raw_rare_results.tsv"), show_col_types = F) |>
    dplyr::rename_with(tolower) |>
    format_rare_results(studies_db$study_extractions$data, studies_db$snp_annotations$data) |>
    dplyr::filter(!is.na(study_extraction_id))

  studies_db$results_metadata$data <- vroom::vroom(file.path(args$results_dir, "results_metadata.tsv"), show_col_types = F) |>
    dplyr::left_join(dplyr::select(studies_db$ld_blocks$data, ld_block, id) |> dplyr::rename(ld_block_id=id), by="ld_block")

  return(studies_db)
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
    dplyr::select(id, unique_study_id, ld_block_id, study_id) |>
    dplyr::rename(study_extraction_id=id)

  snp_annotations_subset <- snp_annotations |>
    dplyr::select(snp, id) |>
    dplyr::rename(snp_id=id)
  
  rare_results <- rare_results |>
    dplyr::mutate(rare_result_group_id=1:dplyr::n()) |>
    dplyr::mutate(candidate_snp=trimws(candidate_snp)) |>
    tidyr::separate_rows(traits, min_ps, genes, files, sep=", ") |>
    dplyr::rename(unique_study_id=traits, min_p=min_ps, known_gene=genes, file=files) |>
    dplyr::mutate(min_p = as.numeric(min_p)) |>
    tidyr::separate(unique_study_id, into = c("study", "ancestry", "chr", "bp"), sep = "_", remove = F) |>
    dplyr::left_join(study_extractions_subset, by="unique_study_id") |>
    dplyr::left_join(snp_annotations_subset, by=c("candidate_snp"="snp"), relationship="many-to-many") |>
    dplyr::select(rare_result_group_id, study_extraction_id, snp_id, ld_block_id, unique_study_id, candidate_snp, study_id, file, chr, bp, min_p, known_gene)

  return(rare_results)
}

# Find relevant snps: all coloc SNPs and all finemapped SNPs that colocalise with nothing (at genome wide significance)
find_relevant_snps <- function(studies_db) {
  colocalising_snps <- dplyr::select(studies_db$colocalisations$data, candidate_snp, snp_id, ld_block) |> dplyr::distinct()

  non_colocalising_snps <- studies_db$study_extractions$data |>
    dplyr::filter(!unique_study_id %in% studies_db$colocalisations$data$unique_study_id & min_p < genome_wide_p_value_threshold) |>
    dplyr::select(snp, snp_id, ld_block) |>
    dplyr::rename(candidate_snp=snp) |>
    dplyr::distinct()
  
  # all_relevant_snps <- dplyr::bind_rows(colocalising_snps, non_colocalising_snps) |> dplyr::distinct()
  return(colocalising_snps)
}

load_data_into_ld_db <- function(ld_conn, studies_db, all_relevant_snps) {
  relevant_ld_blocks <- unique(studies_db$study_extractions$data$ld_block)

  ld_data <- parallel::mclapply(relevant_ld_blocks, mc.cores=max_cores, \(ld_block) {
    tryCatch({
      relevant_snps <- all_relevant_snps |> dplyr::filter(ld_block == ld_block)
      generate_ld_obj(ld_block, relevant_snps$candidate_snp)
    }, error = function(e) {
      message("Error generating LD object for ", ld_block, " - ", e)
      throw(e)
    })
  }) |> dplyr::bind_rows()

  variant_annotations_subset_lead <- dplyr::select(studies_db$snp_annotations$data, snp, id) |> dplyr::rename(lead_snp_id=id)
  variant_annotations_subset_variant <- dplyr::select(studies_db$snp_annotations$data, snp, id) |> dplyr::rename(variant_snp_id=id)
  ld_blocks_subset <- dplyr::select(studies_db$ld_blocks$data, ld_block, id) |> dplyr::rename(ld_block_id=id)

  ld_data <- ld_data |>
    dplyr::left_join(variant_annotations_subset_lead, by=c("lead"="snp")) |>
    dplyr::left_join(variant_annotations_subset_variant, by=c("variant"="snp")) |>
    dplyr::left_join(ld_blocks_subset, by="ld_block") |>
    dplyr::select(-lead, -variant, -ld_block)

  DBI::dbAppendTable(ld_conn, ld_table$name, ld_data)
}

generate_ld_obj <- function(ld_block, snps) {
  message("Generating LD object for ", ld_block)
  ld_file <- file.path(ld_reference_panel_dir, glue::glue("{ld_block}.unphased.vcor1"))
  ld <- vroom::vroom(ld_file, col_names = FALSE, show_col_types = F)
  ldvars <- vroom::vroom(glue::glue("{ld_file}.vars"), delim = " ", col_names = FALSE, show_col_types = F)

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

load_data_into_associations_db <- function(conn, studies_db, all_relevant_snps) {
  message('Retrieving ', nrow(all_relevant_snps), ' SNPs for ', nrow(studies_db$studies$data), ' studies')

  snp_annotations_subset <- studies_db$snp_annotations$data |>
    dplyr::select(snp, id) |>
    dplyr::rename(snp_id=id)

  studies_subset <- studies_db$studies$data |>
    dplyr::select(study_name, id) |>
    dplyr::rename(study_id=id)

  relevant_snps_per_ld_block <- split(all_relevant_snps, all_relevant_snps$ld_block)
  associations <- parallel::mclapply(names(relevant_snps_per_ld_block), mc.cores=max_cores, \(ld_block) {
    gc()
    message("Processing ld_block: ", ld_block)
    tryCatch({
      relevant_snps <- relevant_snps_per_ld_block[[ld_block]]
      imputed_studies_file <- file.path(ld_block_data_dir, ld_block, "imputed_studies.tsv")
      if (!file.exists(imputed_studies_file)) return(NULL)
      
      imputed_studies <- vroom::vroom(imputed_studies_file, show_col_types = F) |>
        dplyr::filter(study %in% studies_db$studies$data$study_name)
      if (nrow(imputed_studies) == 0) return(NULL)

      associations <- apply(imputed_studies, 1, \(study) {
        tryCatch({
          extractions <- vroom::vroom(study[['file']],
            show_col_types = F,
            altrep = FALSE,
            col_select = dplyr::any_of(c("SNP", "BETA", "SE", "P", "EAF", "IMPUTED"))
          )
          if (!"IMPUTED" %in% names(extractions)) {
            extractions$IMPUTED <- FALSE
          }
          extractions <- extractions |>
            dplyr::filter(SNP %in% relevant_snps$candidate_snp) |>
            dplyr::mutate(study = study[['study']]) |>
            dplyr::rename_with(tolower)
          return(extractions)
        }, error = function(e) {
          message('Error processing file: ', study[['file']], ' - ', e)
          return(NULL)
        })
      })
      associations <- associations[!sapply(associations, is.null)]
      if (length(associations) == 0) return(NULL)
      associations <- do.call(rbind, associations)

      associations <- associations |> 
        dplyr::left_join(snp_annotations_subset, by="snp") |>
        dplyr::left_join(studies_subset, by=c("study"="study_name")) |>
        dplyr::select(snp_id, study_id, beta, se, imputed, p, eaf)

      message('Extracted ', nrow(associations), ' associations for ', ld_block)

      # output_file <- glue::glue("{args$results_dir}/{gsub('[/:]', '_', ld_block)}_associations.tsv.gz")
      # vroom::vroom_write(associations, output_file)
      return(associations)
    }, error = function(e) {
      message('Error processing ld_block: ', ld_block, ' - ', e)
      return(NULL)
    })
  }) 

  gc()
  associations <- associations[!sapply(associations, is.null)]
  associations <- do.call(rbind, associations)

  DBI::dbAppendTable(conn, associations_table$name, associations)
  num_rows <- DBI::dbGetQuery(conn, glue::glue("SELECT COUNT(*) FROM {associations_table$name}"))
  message('Added ', num_rows$count_star, ' rows to ', associations_table$name)
  if (num_rows$count_star > 1000000) {
    association_files <- Sys.glob(glue::glue("{args$results_dir}/*_associations.tsv.gz"))
    file.remove(association_files)
  }
}


main()
