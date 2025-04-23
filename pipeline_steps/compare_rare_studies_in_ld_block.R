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
  
  studies_to_compare <- lapply(standardised_studies$file, 
    function(file) vroom::vroom(file, show_col_types = F) |>
      dplyr::select(dplyr::all_of(standardised_gwas_columns))
  )
  names(studies_to_compare) <- standardised_studies$unique_study_id

  # Add unique study id as column
  studies_to_compare <- lapply(names(studies_to_compare), function(name){
    study <- studies_to_compare[[name]]
    study$unique_study_id <- name
    return(study)
  })
  # Re-assign names
  names(studies_to_compare) <- standardised_studies$unique_study_id

  # All rare variants available in ld_block
  variants <- do.call("rbind", lapply(studies_to_compare, function(x) {
    x |> dplyr::select(SNP, P)})) |> 
    dplyr::group_by(SNP) |> 
    dplyr::summarise(min_P = min(P))  

  # Compare rare variant hits across studies at given p-value threshold
  compare_results <- compare_by_variant(studies = studies_to_compare, variants = variants, P_thresh = lowest_rare_p_value_threshold)

  if(nrow(compare_results) == 0){
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
  # Retain variants under the p-value threshold to compare
  variants_keep <- variants |> dplyr::filter(min_P <= P_thresh)|> dplyr::pull(SNP)
  compare_wide <- data.frame()
  rare_studies <- data.frame()

  # Pull studies with shared varaiants
  for (var in variants_keep) {
    found_studies <- sapply(studies, function(study) {
      study <- dplyr::filter(study, SNP == var & P <= P_thresh) |>
        dplyr::pull(unique_study_id, P, GENE)
    })

    if (length(found_studies$unique_study_id) == 0){
      next
    }
    found_study_ids <- unlist(found_studies$unique_study_id)
    found_study_ids <- found_study_ids[order(found_study_ids)]

    if (length(found_study_ids) == 1){
      next
    }
    
    rare_studies <- rbind(rare_studies, found_studies)
    compare_wide <- rbind(compare_wide, data.frame(
      traits = paste(found_study_ids, collapse = ", "), candidate_snp = var))
  }
  
  if(nrow(compare_wide) == 0){
    return(data.frame())
  }

  colnames(compare_wide) <- c("traits","candidate_snp")
  return(compare_wide)
}

invisible(main())