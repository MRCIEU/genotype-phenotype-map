source('../pipeline_steps/constants.R')

update_studies_processed <- function() {
  gene_name_map <- vroom::vroom(paste0('/Users/wt23152/Documents/Projects/genotype-phenotype-map/pipeline_steps/scratch/data/1000genomes/gene_name_map.tsv'), show_col_types=F)
  results_dir <- '/Users/wt23152/Documents/Projects/genotype-phenotype-map/pipeline_steps/scratch/results'
  studies_processed_file <- paste0(results_dir, '/studies_processed.tsv')
  studies_processed <- vroom::vroom(studies_processed_file)
  gene_names <- gene_name_map$GENE_NAME[match(studies_processed$gene, gene_name_map$ENSEMBL_ID)]
  studies_processed$gene <- gene_names
  vroom::vroom_write(studies_processed, studies_processed_file)

}

update_imputed_studies <- function() {
  ld_regions <- vroom::vroom('../pipeline_steps/data/ld_regions.tsv')
  ld_info <- construct_ld_block(ld_regions$ancestry, ld_regions$chr, ld_regions$start, ld_regions$stop)

  lapply(ld_info$ld_block_data, function(ld_block) {
    extracted_studies_file <- paste0(ld_block, '/extracted_studies.tsv')
    imputed_studies_file <- paste0(ld_block, '/imputed_studies.tsv')
    if (!file.exists(extracted_studies_file)) return()

    extracted_studies <- vroom::vroom(extracted_studies_file, show_col_types = F)
    imputed_studies <- vroom::vroom(imputed_studies_file, col_types = vroom::cols(
      chr = vroom::col_character(),
      bp = vroom::col_number(),
      p_value_threshold = vroom::col_number(),
      eaf_from_reference_panel = vroom::col_logical(),
      sample_size = vroom::col_number()
      ), show_col_types = F)
    if (nrow(extracted_studies) == nrow(imputed_studies)) return()
    print(paste(ld_block, 'sizes are different...', nrow(extracted_studies) - nrow(imputed_studies)))

    add_to_imputed <- apply(extracted_studies, 1, function(extracted_study) {
      already_in_imputed <- dplyr::filter(imputed_studies, study == extracted_study['study'])
      if (nrow(already_in_imputed) > 0) return()

      imputed_study <- sub('original', 'imputed', extracted_study[['file']])
      if (!file.exists(imputed_study)) return()

      eaf_from_reference_panel <- ifelse(grepl('GTEx', extracted_study['study']), TRUE, NA)

      return(data.frame(study=extracted_study['study'], file=extracted_study['file'], chr=extracted_study['chr'], bp=as.numeric(extracted_study['bp']),
                        p_value_threshold=as.numeric(extracted_study['p_value_threshold']), sample_size=as.numeric(extracted_study['sample_size']), cis_trans=extracted_study['cis_trans'],
                        rows_imputed=NA, eaf_from_reference_panel=eaf_from_reference_panel
                        )
      )

    }) |> dplyr::bind_rows()

    imputed_studies <- dplyr::bind_rows(imputed_studies, add_to_imputed) |> dplyr::distinct()
    vroom::vroom_write(imputed_studies, imputed_studies_file)
  })

}

update_finemapped_studies <- function() {
  ld_regions <- vroom::vroom('../pipeline_steps/data/ld_regions.tsv')
  ld_regions <- ld_regions[ld_regions$chr >= 12, ]
  ld_info <- construct_ld_block(ld_regions$ancestry, ld_regions$chr, ld_regions$start, ld_regions$stop)

  lapply(ld_info$ld_block_data, function(ld_block) {
    finemapped_studies_file <- paste0(ld_block, '/finemapped_studies.tsv')
    imputed_studies_file <- paste0(ld_block, '/imputed_studies.tsv')
    if (!file.exists(finemapped_studies_file)) return()
    print(ld_block)

    imputed_studies <- vroom::vroom(imputed_studies_file, show_col_types=F) 
    finemapped_studies <- vroom::vroom(finemapped_studies_file, col_types = vroom::cols(
      chr = vroom::col_character(),
      bp = vroom::col_number(),
      min_p = vroom::col_number(),
      p_value_threshold = vroom::col_number(),
      sample_size = vroom::col_number()
      ), show_col_types = F)

    if (nrow(imputed_studies) == 0) return()
    add_to_finemapped <- apply(imputed_studies, 1, function(imputed_study) {
      finemap_file_prefix <- sub('imputed', 'finemapped', imputed_study[['file']])
      finemap_file_prefix <- sub('\\..*', '', finemap_file_prefix)
      finemap_files <- Sys.glob(paste0(finemap_file_prefix, '*'))

      if (length(finemap_files) == 0) return()
      if (nrow(dplyr::filter(finemapped_studies, grepl(finemap_file_prefix, file))) > 0) return()

      unique_study_ids <- file_prefix(finemap_files)
      message <- if(length(finemap_files) > 1) 'success' else NA

      bps <- lapply(finemap_files, function(file) {
        finemap <- vroom::vroom(file, show_col_types=F)
        min <- finemap[which.min(finemap$P), ]
        return(data.frame(bp=min$BP, min_p=min$P))
      }) |> dplyr::bind_rows()

      print(unique_study_ids)
      print(finemap_files)
      print(bps$bp)
      print(bps$min_p)

      return(data.frame(study=imputed_study['study'], unique_study_id=unique_study_ids, file=finemap_files, chr=imputed_study['chr'], bp=as.numeric(bps$bp),
                        p_value_threshold=as.numeric(imputed_study['p_value_threshold']), sample_size=as.numeric(imputed_study['sample_size']), cis_trans=imputed_study['cis_trans'],
                        message=message, min_p=bps$min_p
                        )
      )


    }) |> dplyr::bind_rows()

    if (nrow(add_to_finemapped) == 0) return()
    print(paste('missing imputed from finemapped',nrow(add_to_finemapped)))

    finemapped_studies <- dplyr::bind_rows(finemapped_studies, add_to_finemapped) |> dplyr::distinct()
    vroom::vroom_write(finemapped_studies, finemapped_studies_file)
  })
}

populate_json_for_besd <- function() {
  gtex <- vroom::vroom('/local-scratch/projects/genotype-phenotype-map/GTEx_portal.csv')

  gtex$Tissue <- gsub('[ -]', '_', gtex$Tissue)
  gtex$Tissue <- gsub('_+', '_', gtex$Tissue)
  gtex$Tissue <- gsub('[()]', '', gtex$Tissue)

  apply(gtex, 1, function(tissue) {
    metadata <- list(sample_size=as.numeric(tissue['# RNASeq and Genotyped samples']), ancestry='EUR', cis_trans='cis', category='continuous')
    metadata_json <- jsonlite::toJSON(metadata, auto_unbox = T, pretty = T)
    file_name <- paste0('/local-scratch/data/GTEx-cis/', tissue['Tissue'], '.json')
    write(metadata_json, file_name)
  })

}

update_finemapped_studies()
