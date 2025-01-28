source('constants.R')

# In the absence of a complex LD structure surrounding rare variant associations, compare WES and WGS studies by variant identifer alone
# for all studies with a significant hit within a given LD block

parser <- argparser::arg_parser('Compare rare variant sequencing studies per LD block')
parser <- argparser::add_argument(parser, '--ld_block', help = 'LD block', type = 'character')
parser <- argparser::add_argument(parser, '--completed_output_file', help = 'Result file to save', type = 'character')
args <- argparser::parse_args(parser)

p_value_threshold = 0.00015

main <- function() {
    ld_info <- ld_block_dirs(args$ld_block)
    block <- vroom::vroom(glue::glue('{pipeline_metadata_dir}updated_ld_blocks_to_colocalise.tsv'), show_col_types=F) |>
        dplyr::filter(data_dir == ld_info$ld_block_data)

    standardised_file <- glue::glue('{ld_info$ld_block_data}/standardised_studies.tsv')
    
    if (file.exists(standardised_file)) {
    standardised_studies <- vroom::vroom(standardised_file, delim = '\t', show_col_types = F) |>
        dplyr::filter(varaint_type == "rare")
    }

    if (!file.exists(standardised_file) || nrow(block) == 0 || nrow(standardised_studies) == 0) {
        message(glue::glue('Nothing to process for LD region {ld_info$ld_block_data}, skipping.'))
        vroom::vroom_write(data.frame(), args$completed_output_file)
        return()
    }

    standardised_studies$unique_study_id <- glue::glue('{standardised_studies$study}_{file_prefix(standardised_studies$file)}')
    
    studies_to_compare <- lapply(standardised_studies$file, function(file) vroom::vroom(file, delim = '\t', show_col_types = F))
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
    variants <- do.call("rbind", lapply(studies_to_compare, function(x,y){x |> dplyr::select(SNP, P)})) |> 
        dplyr::group_by(SNP) |> 
        dplyr::summarise(min_P = min(P))    

    # Compare rare variant hits across studies at given p-value threshold
    compare_results <- compare_by_variant(studies = studies_to_compare, variants = variants, P_thresh = p_value_threshold)

    if(nrow(compare_results) == 0){
       vroom::vroom_write(data.frame(), args$completed_output_file)
       return() 
    }

    compare_results_file <- glue::glue('{ld_info$ld_block_data}/compare_rare_results.tsv')
    vroom::vroom_write(compare_results, compare_results_file)
    vroom::vroom_write(data.frame(), args$completed_output_file)
}

compare_by_variant <- function(studies, variants, P_thresh){
    
    # Retain variants under the p-value threshold to compare
    variants_keep <- variants |> dplyr::filter(min_P <= P_thresh)|> dplyr::pull(SNP)

    # Output data frame
    compare_wide <- data.frame()

    # Pull studies with shared varaiants
    for(var in variants_keep){
        found_studies <- do.call("rbind",lapply(studies, function(x){x |> dplyr::filter(SNP == var & P <= P_thresh)}))

        if(nrow(found_studies) == 0){
            next
        }

        found_study_ids <- found_studies$unique_study_id
        found_study_ids <- found_study_ids[order(found_study_ids)]

        if(length(found_study_ids) == 1){
            next
        }
      
        # Wide format output
        compare_wide <- rbind(compare_wide, data.frame(traits = paste(found_study_ids, collapse = ","), candidate_snp = var))    
    }
    
    if(nrow(compare_wide) == 0){
        return(data.frame())
    }

    colnames(compare_wide) <- c("traits","candidate_snp")
    return(compare_wide)
}

main()