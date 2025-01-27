source('constants.R')

# In the absence of a complex LD structure surrounding rare variant associations, compare WES and WGS studies by variant identifer alone
# for all studies with a significant hit within each LD block

parser <- argparser::arg_parser('Compare rare variant sequencing studies per LD block')
parser <- argparser::add_argument(parser, '--ld_block', help = 'LD block', type = 'character')
parser <- argparser::add_argument(parser, '--completed_output_file', help = 'Result file to save', type = 'character')
args <- argparser::parse_args(parser)

lp_thresh =  -log10(0.00015)
test_dir = "/local-scratch/projects/genotype-phenotype-map/test/results/rare_variant_compare"

main <- function() {
    ## ----- FOR TESTING ONLY
    args <- list(ld_block="EUR/1/95684841-97259500")
    ## -----

    ld_info <- ld_block_dirs(args$ld_block)
    #block <- vroom::vroom(glue::glue('{pipeline_metadata_dir}updated_ld_blocks_to_colocalise.tsv'), show_col_types=F) |>
    #    dplyr::filter(data_dir == ld_info$ld_block_data)

    ## ----- FOR TESTING ONLY
    standardised_file <- "data/example_studiesinblock.txt"
    standardised_studies <- vroom::vroom(standardised_file, delim = '\t', show_col_types = F)
    colnames(standardised_studies) <- "file"
    standardised_studies$study <- lapply(stringr::str_split(standardised_studies$file, pattern = "/"), function(x){x[7]}) |> unlist()
    ## -----

    #standardised_file <- glue::glue('{ld_info$ld_block_data}/standardised_studies.tsv')
    
    #if (file.exists(standardised_file)) {
    #standardised_studies <- vroom::vroom(standardised_file, delim = '\t', show_col_types = F) |>
    #    dplyr::filter(varaint_type == "rare")
    #}

    #if (!file.exists(standardised_file) || nrow(block) == 0 || nrow(standardised_studies) == 0) {
    if (!file.exists(standardised_file) || nrow(standardised_studies) == 0) {
        message(glue::glue('Nothing to process for LD region {ld_info$ld_block_data}, skipping.'))
        #vroom::vroom_write(data.frame(), args$completed_output_file)
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
    names(studies_to_compare) <- standardised_studies$unique_study_id

    # Test output directory
    compare_results_dir <- glue::glue('{test_dir}/{gsub("/","_",args$ld_block)}')
    dir.create(compare_results_dir, recursive = T, showWarnings = F)

    # Make study metadata
    all_studies <- make_metadata(studies = studies_to_compare)
    all_studies <- all_studies |> 
        dplyr::inner_join(standardised_studies) |>
        dplyr::select(study, unique_study_id, file, chr, bp, min_p)

    # Keep only the file named by the top hit in the region (some duplicate files exists if multiple hits are found in a ld block)
    all_studies <- all_studies |> dplyr::arrange(study, min_p) |> dplyr::filter(duplicated(study) == FALSE)
    vroom::vroom_write(all_studies, glue::glue('{compare_results_dir}/all_study_blocks.tsv'))

    # Remove duplicate files from list of study data
    studies_to_compare <- studies_to_compare[all_studies$unique_study_id]

    # All rare variants available in ld_block
    variants <- do.call("rbind", lapply(studies_to_compare, function(x,y){x |> dplyr::select(NAME, LP)})) |> 
        dplyr::group_by(NAME) |> 
        dplyr::summarise(max_LP = max(LP))    

    # Compare rare variant hits across studies
    compare_results <- compare_by_variant(studies = studies_to_compare, variants = variants, LP_thresh = lp_thresh)

    if(length(compare_results) == 0){
       #vroom::vroom_write(data.frame(), args$completed_output_file)
       return() 
    }

    vroom::vroom_write(compare_results$compare_long, glue::glue('{compare_results_dir}/compare_results.tsv'))
    vroom::vroom_write(compare_results$compare_wide, glue::glue('{compare_results_dir}/raw_compare_results.tsv'))
    #vroom::vroom_write(data.frame(), args$completed_output_file)

    # Some metadata on number of comparisons at different -log(P) thresholds
    lp_results <- compare_across_lp(lp_list = list(5,4,3.82,3,2), studies = studies_to_compare, variants = variants)
    vroom::vroom_write(lp_results, glue::glue('{compare_results_dir}/shared_by_lp.tsv'))
}

make_metadata <- function(studies){

    studies_tophit <- do.call("rbind", lapply(studies, function(x){
        # Report P-value for variant in filename
        top <- x|> dplyr::filter(BP == unique(gsub(".*_","",unique_study_id)))
        return(data.frame(unique_study_id = top$unique_study_id, chr = top$CHR, bp = top$BP, min_p = 10^-top$LP))
    }))

    return(studies_tophit)
}

compare_by_variant <- function(studies, variants, LP_thresh){
    
    # Retain variants passing -log10(P) threshold
    variants_keep <- variants |> dplyr::filter(max_LP >= LP_thresh)|> dplyr::pull(NAME)

    # Output data frames
    compare_long <- data.frame()
    compare_wide <- data.frame()

    # For variants exceeding theshold -log(P) in at least one study
    for(var in variants_keep){
        found_studies <- do.call("rbind",lapply(studies, function(x){x |> dplyr::filter(NAME == var & LP >= LP_thresh)}))

        if(nrow(found_studies) == 0){
            next
        }

        found_study_ids <- found_studies$unique_study_id
        found_study_ids <- found_study_ids[order(found_study_ids)]

        if(length(found_study_ids) == 1){
            next
        }

        # Long format output
        compare_long <- rbind(compare_long, as.data.frame(cbind(t(combn(found_study_ids, 2)), var)))        
        # Wide format output
        compare_wide <- rbind(compare_wide, data.frame(traits = paste(found_study_ids, collapse = ","), candidate_snp = var))    
    }
    
    if(length(compare_long) == 0){
        return(list())
    }

    colnames(compare_long) <- c("unique_study_a","unique_study_b","candidate_snp")
    colnames(compare_wide) <- c("traits","candidate_snp")
    return(list(compare_long = compare_long, compare_wide = compare_wide))
}

compare_across_lp <- function(lp_list, studies, variants){
    
    # Pairwise variant comparison at differing -log10(P) thresholds 
    lp_results <- lapply(lp_list, compare_by_variant, studies = studies_to_compare, variants = variants)
    lp_results <- lapply(lp_results, function(x){
        if(length(x) > 0){
            return(x[[1]])
        }else{return(NULL)}})
    names(lp_results) <- paste0("lp_",unlist(lp_list))
    
    # Count number of shared hits between pairs of studies

    # Dummy dataframe of all pairs of studies
    df <- data.frame(t(combn(names(studies), 2)))
    names(df) <- c("unique_study_a", "unique_study_b")

    lp_results_count <- lapply(lp_results, function(x){
        if(is.null(x)){
            df <- data.frame(t(combn(names(studies), 2)), NA)
            names(df) <- c("unique_study_a", "unique_study_b", "count")
            return(df)
        }else{
            x |> dplyr::group_by(unique_study_a,unique_study_b) |> dplyr::summarise(count = dplyr::n())}})

    joined_results <- purrr::reduce(names(lp_results_count), function(acc, name){
        if (is.null(acc)) {
            lp_results_count[[name]]
        } else {
            dplyr::left_join(acc, lp_results_count[[name]], by = c("unique_study_a","unique_study_b"), suffix = c("", paste0("_", name)))
        }}, .init = df)

    names(joined_results) <- c("unique_study_a", "unique_study_b", paste0("count_lp_", unlist(lp_list)))

    joined_results[is.na(joined_results)] <- 0

    # Remove studies with no shared variants at any threshold
    joined_results <- joined_results[rowSums(joined_results[,paste0("count_lp_", unlist(lp_list))]) > 0,]
    return(joined_results)
}

main()