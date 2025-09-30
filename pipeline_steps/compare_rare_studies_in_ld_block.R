source('constants.R')

# In the absence of a complex LD structure surrounding rare variant associations, compare WES and WGS studies by variant identifer alone
# for all studies with a significant hit within a given LD block

parser <- argparser::arg_parser('Compare rare variant sequencing studies per LD block')
parser <- argparser::add_argument(parser, '--ld_block', help = 'LD block', type = 'character')
parser <- argparser::add_argument(parser, '--completed_output_file', help = 'Result file to save', type = 'character')
args <- argparser::parse_args(parser)

main <- function() {
  ld_info <- ld_block_dirs(args$ld_block)
  block <- vroom::vroom(glue::glue('{pipeline_metadata_dir}updated_ld_blocks_to_colocalise.tsv'), show_col_types=F) |>
    dplyr::filter(data_dir == ld_info$ld_block_data)

  compare_rare_cache_file <- glue::glue('{ld_info$ld_block_data}/compare_rare_cached_studies.tsv')
  standardised_file <- glue::glue('{ld_info$ld_block_data}/standardised_studies.tsv')
  if (file.exists(standardised_file)) {
  standardised_studies <- vroom::vroom(standardised_file, delim = '\t', show_col_types = F) |>
    dplyr::filter(variant_type != variant_types$common)
  }

  if (!file.exists(standardised_file) || nrow(block) == 0 || nrow(standardised_studies) == 0) {
    message(glue::glue('No rare studies to compare in LD region {ld_info$ld_block_data}, skipping.'))
    vroom::vroom_write(data.frame(), args$completed_output_file)
    return()
  }

  already_run <- already_run_comparison(compare_rare_cache_file, standardised_studies)
  if (already_run) {
    vroom::vroom_write(data.frame(), args$completed_output_file)
    return()
  }

  standardised_studies$unique_study_id <- glue::glue('{standardised_studies$study}_{file_prefix(standardised_studies$file)}')
  
  studies_to_compare <- lapply(standardised_studies$file, function(file) {
    study <- vroom::vroom(file, show_col_types = F) |>
      dplyr::select(dplyr::any_of(standardised_gwas_columns))
    if (!"GENE" %in% names(study)) {
      study$GENE <- NA
    }
    study$GENE <- gsub("'", "", study$GENE)
    return(study)
  })
  names(studies_to_compare) <- standardised_studies$unique_study_id

  # Add unique study id as column
  studies_to_compare <- lapply(names(studies_to_compare), function(name){
    study_name <- strsplit(name, '_')[[1]][1]
    study <- studies_to_compare[[name]]
    study$unique_study_id <- name
    study$file <- standardised_studies$file[match(study_name, standardised_studies$study)]
    return(study)
  })
  # Re-assign names
  names(studies_to_compare) <- standardised_studies$unique_study_id

  message(glue::glue('Comparing {length(studies_to_compare)} rare variants in LD region {ld_info$ld_block_data}'))

  # All rare variants available in ld_block
  variants <- do.call("rbind", lapply(studies_to_compare, function(x) {
    x |> dplyr::select(SNP, P)})) |> 
    dplyr::group_by(SNP) |> 
    dplyr::summarise(min_P = min(P))  
  
  message(glue::glue('Found {nrow(variants)} rare variants in LD region {ld_info$ld_block_data}'))

  # Compare rare variant hits across studies at given p-value threshold
  compare_results <- compare_by_variant(studies = studies_to_compare, variants = variants, P_thresh = lowest_p_value_threshold)

  if(nrow(compare_results) == 0){
    message(glue::glue('No rare variants to compare in LD region {ld_info$ld_block_data}, skipping.'))
     vroom::vroom_write(data.frame(), args$completed_output_file)
     return() 
  }

  studies_to_cache <- data.frame(study=standardised_studies$study)

  compare_results_file <- glue::glue('{ld_info$ld_block_data}/compare_rare_results.tsv')
  vroom::vroom_write(compare_results, compare_results_file)
  vroom::vroom_write(studies_to_cache, compare_rare_cache_file)
  vroom::vroom_write(data.frame(), args$completed_output_file)
}

already_run_comparison <- function(compare_rare_cache_file, standardised_studies) {
  if (file.exists(compare_rare_cache_file)) {
    compare_rare_cache <- vroom::vroom(compare_rare_cache_file, delim = '\t', show_col_types = F)
    already_run <- setdiff(standardised_studies$study, compare_rare_cache$study)
    if (length(already_run) == 0) {
      message('No new rare studies to compare, skipping.')
      return(TRUE)
    }
  }
  return(FALSE)
}

compare_by_variant <- function(studies, variants, P_thresh) {
  variants_keep <- variants |> dplyr::filter(min_P <= P_thresh)|> dplyr::pull(SNP)
  compare_wide <- data.frame()

  for (var in variants_keep) {
    found_studies <- lapply(studies, function(study) {
      study <- dplyr::filter(study, SNP == var & P <= P_thresh) |>
        dplyr::select(unique_study_id, BP, P, GENE, file) |>
        dplyr::mutate(GENE = as.character(GENE))
    }) |> dplyr::bind_rows()

    if (nrow(found_studies) <= 1) {
      next
    }

    #Update unique study id to include the new BP value, to avoid duplicates
    new_bp <- strsplit(var, ':')[[1]][2]
    new_bp <- gsub('_', '-', new_bp)
    found_studies$unique_study_id <- sub('(.*)_\\d+', paste0('\\1_', new_bp), found_studies$unique_study_id)

    compare_wide <- rbind(compare_wide, data.frame(
      traits = paste(found_studies$unique_study_id, collapse = ", "),
      candidate_snp = var,
      min_ps = paste(found_studies$P, collapse = ", "),
      genes = paste(found_studies$GENE, collapse = ", "),
      ld_block = args$ld_block,
      files = paste(found_studies$file, collapse = ", "))
    )
  }

  return(compare_wide)
}

main()