source('constants.R')

parser <- argparser::arg_parser('Standardise GWAS for pipeline')
parser <- argparser::add_argument(parser, '--ld_block', help = 'LD block that the ', type = 'character')
parser <- argparser::add_argument(parser, '--completed_output_file', help = 'Completed output file', type = 'character')
args <- argparser::parse_args(parser)

main <- function(args) {
  ld_info <- ld_block_dirs(args$ld_block)
  ld_region <- vroom::vroom(paste0(ld_info$ld_matrix_prefix, '.tsv'), show_col_types = F)

  extracted_studies_file <- paste0(ld_info$ld_block_data, '/extracted_studies.tsv')
  extracted_studies  <- vroom::vroom(extracted_studies_file , show_col_types = F)

  standardised_studies_file <- paste0(ld_info$ld_block_data, '/standardised_studies.tsv')
  if (file.exists(standardised_studies_file)) {
    existing_standardised_studies <- vroom::vroom(standardised_studies_file,
                                                  show_col_types = F,
                                                  col_types = vroom::cols(
                                                    chr = vroom::col_character(),
                                                    bp = vroom::col_number(),
                                                    sample_size = vroom::col_number(),
                                                    p_value_threshold = vroom::col_number(),
                                                    snps_removed_by_reference_panel=vroom::col_number(),
                                                    snps_removed_by_dentist=vroom::col_number()
                                                  )
    )
  } else {
    existing_standardised_studies <- empty_standardised_studies()
  }

  if (nrow(extracted_studies) > 0) {
    standardised_studies <- apply(extracted_studies, 1, function (study) {
      result <- perform_standardisation(study, ld_region)
      result <- perform_qc(result$gwas, result$study)
      vroom::vroom_write(result$gwas)

      return(result$study)
    }) |> dplyr::bind_rows()
    if (nrow(standardised_studies) > 0) dplyr::mutate_at(standardised_studies, c('bp', 'p_value_threshold', 'sample_size',
                                                                                 'snps_removed_by_reference_panel', 'snps_removed_by_dentist'), as.numeric)

  }

  standardised_studies <- dplyr::bind_rows(existing_standardised_studies, standardised_studies) |>
    dplyr::distinct()

  vroom::vroom_write(standardised_studies, standardised_studies_file)
  vroom::vroom_write(data.frame(), args$completed_output_file)
}

perform_qc <- function(gwas, study) {
  dentist_gwas <- dplyr::rename(gwas, REAL_SNP='SNP', SNP='RSID', A1='EA', A2='OA', freq='EAF', beta='BETA', se='SE', p='P') |>
    dplyr::mutate(N = as.numeric(study[['sample_size']])) |>
    dplyr::select(SNP, A1, A2, freq, beta, se, p, N)

  dentist_tmp_file <- tempfile(basename(study[['file']]))
  vroom::vroom_write(dentist_gwas, dentist_tmp_file)

  dentist_command <- paste('DENTIST --bfile', paste0(thousand_genomes_dir, study[['ancestry']]),
                           '--gwas-summary', dentist_tmp_file,
                           '--chrID', study[['chr']],
                           '--out', dentist_tmp_file
  )
  system(dentist_command, wait = T)

  dentist_file_to_remove <- paste0(dentist_tmp_file, '.DENTIST.short.txt')
  if (!file.exists(dentist_file_to_remove)) {
    stop('DENTIST command failed')
  }
  dentist_to_remove <- vroom::vroom(dentist_file_to_remove, col_names = F, show_col_types = F)
  gwas <- dplyr::filter(gwas, RSID %in% dentist_to_remove$X1)

  study['snps_removed_by_dentist'] <- nrow(dentist_to_remove)
  return(list(gwas=gwas, study=study))
}

perform_standardisation <- function(study, ld_region) {
  standardised_file <- sub('original', 'standardised', study[['file']])
  if (standardised_file %in% existing_standardised_studies$file) {
    return()
  }

  gwas <- vroom::vroom(study[['file']], show_col_types = F)

  response <- convert_reference_build_via_liftover(gwas, study[['reference_build']], reference_builds$GRCh37) |>
    standardise_alleles() |>
    standardise_extracted_gwas(ld_region)

  study['file'] <- standardised_file
  study['eaf_from_reference_panel'] <- response$eaf_from_reference_panel
  study['snps_removed_by_reference_panel'] <- response$snps_removed_by_reference_panel
  study <- study[-match('reference_build', names(study))]

  #maybe cant save it as gzipped file?
  vroom::vroom_write(response$gwas, study[['file']])

  return(list(gwas=response$gwas, study=as.list(study)))
}

empty_standardised_studies <- function() {
  return(data.frame(study=character(),
                     file=character(),
                     ancestry=character(),
                     chr=character(),
                     bp=numeric(),
                     p_value_threshold=numeric(),
                     category=character(),
                     sample_size=numeric(),
                     cis_trans=character(),
                     eaf_from_reference_panel=logical(),
                     snps_removed_by_reference_panel=numeric(),
                     snps_removed_by_dentist=numeric(),
  ))

}

standardise_extracted_gwas <- function(gwas, ld_region) {
  eaf_from_reference_panel <- FALSE
  original_gwas_size <- nrow(gwas)
  gwas <- dplyr::distinct(gwas, CHR, BP, EA, OA, .keep_all = TRUE)
  
  if (!"Z" %in% colnames(gwas)) {
    gwas <- dplyr::mutate(gwas, SE = replace(SE, SE == 0, 0.00001), Z = BETA / SE)
  }
  
  if (!"P" %in% colnames(gwas) && "LP" %in% colnames(gwas)) {
    gwas <- dplyr::mutate(gwas, P = 10^(-LP)) |>
      dplyr::select(-LP)
  }

  gwas <- dplyr::filter(gwas, SNP %in% ld_region$SNP)
  ld_region <- dplyr::filter(ld_region, SNP %in% gwas$SNP)

  columns_to_coerce <- c("EAF") # Add BETA and SE if needed
  gwas <- dplyr::mutate(gwas, dplyr::across(dplyr::all_of(columns_to_coerce), as.numeric))
  
  if (all(is.na(gwas$EAF))) {
    gwas <- dplyr::select(gwas, -EAF) |>
      dplyr::left_join(ld_region |> dplyr::select(SNP, EAF), by = "SNP")
    eaf_from_reference_panel <- TRUE
  }

  gwas <- tidyr::drop_na(gwas, dplyr::all_of(columns_to_coerce)) |>
    dplyr::arrange(match(SNP, ld_region$SNP))
  
  return(list(gwas = gwas,
              eaf_from_reference_panel = eaf_from_reference_panel,
              snps_removed_by_reference_panel = nrow(gwas) - original_gwas_size))
}

standardise_alleles <- function(gwas) {
  gwas$EA <- toupper(gwas$EA)
  gwas$OA <- toupper(gwas$OA)

  to_flip <- (gwas$EA > gwas$OA) & (!gwas$EA %in% c("D", "I"))
  if (any(to_flip)) {
    gwas$EAF[to_flip] <- 1 - gwas$EAF[to_flip]
    if ('BETA' %in% names(gwas)) {
      gwas$BETA[to_flip] <- -1 * gwas$BETA[to_flip]
    }
    if ('Z' %in% names(gwas)) {
      gwas$Z[to_flip] <- -1 * gwas$Z[to_flip]
    }

    temp <- gwas$OA[to_flip]
    gwas$OA[to_flip] <- gwas$EA[to_flip]
    gwas$EA[to_flip] <- temp
  }

  gwas$SNP <- toupper(paste0(gwas$CHR, ":", format(gwas$BP, scientific = F, trim = T), "_", gwas$EA, "_", gwas$OA))

  return(gwas)
}


#' convert_reference_build_via_liftover: Change reference build of BP marker from allow list of liftOver conversions
#'
#' @param gwas: GWAS (file or dataframe) of standardised GWAS
#' @param input_reference_build: string reference build, found in reference_builds list
#' @param output_reference_build: string reference build that GWAS is to change to, found in reference_builds list
#' @return gwas input is altered and returned
convert_reference_build_via_liftover <- function(gwas,
                                                 input_reference_build=reference_builds$GRCh37,
                                                 output_reference_build=reference_builds$GRCh38) {
  if (input_reference_build == output_reference_build) {
    return(gwas)
  }

  liftover_conversion <- available_liftover_conversions[[paste0(input_reference_build, output_reference_build)]]
  if (is.null(liftover_conversion)) {
    stop(paste(c("Error: liftOver combination of", input_build, output_build, "not recocognised.",
                 "Reference builds must be one of:", reference_builds), collapse = " "))
  }

  original_gwas_size <- nrow(gwas)

  bed_file_input <- tempfile(fileext = ".bed")
  bed_file_output <- tempfile(fileext = ".bed")
  unmapped <- tempfile(fileext = ".unmapped")

  create_bed_file_from_gwas(gwas, bed_file_input)
  run_liftover(bed_file_input, bed_file_output, input_reference_build, output_reference_build, unmapped)
  gwas <- use_bed_file_to_update_gwas(gwas, bed_file_output)

  updated_gwas_size <- nrow(gwas)
  if (updated_gwas_size < original_gwas_size) {
    message(paste("Warning: During liftover conversion, the GWAS lost", original_gwas_size - updated_gwas_size,
                  "rows out of ", original_gwas_size))
  }

  return(gwas)
}

create_bed_file_from_gwas <- function(gwas, output_file) {
  bed_format <- tibble::tibble(
    CHR = paste0("chr", gwas$CHR),
    BP1 = gwas$BP,
    BP2 = gwas$BP+1,
    CHRBP = paste0(CHR, ":", BP1, "-", BP2)
  )

  vroom::vroom_write(bed_format, output_file, col_names=F, delim=" ")
  return(bed_format)
}

run_liftover <- function(bed_file_input, bed_file_output, input_build, output_build, unmapped) {
  lifover_binary <- paste0(liftover_dir, "liftOver")
  liftover_conversion <- available_liftover_conversions[[paste0(input_build, output_build)]]

  chain_file <- paste0(liftover_dir, liftover_conversion)
  liftover_command <- paste(lifover_binary, bed_file_input, chain_file, bed_file_output, unmapped)
  system(liftover_command, wait=T)
}

use_bed_file_to_update_gwas <- function(gwas, bed_file) {
  liftover_bed <- vroom::vroom(bed_file, col_names=F, show_col_types=F)
  liftover_bed$X1 <- gsub("chr", "", liftover_bed$X1)
  liftover_bed$X4 <- sub(".*:(\\d+)-.*", "\\1", liftover_bed$X4, perl=T)

  bed_map <- tibble::tibble(
    NEW_BP = liftover_bed$X2,
    ORIGINAL_CHRBP = paste(liftover_bed$X1, liftover_bed$X4, sep=":")
  )

  matching <- match(gwas$CHRBP, bed_map$ORIGINAL_CHRBP)
  gwas$BP <- bed_map$NEW_BP[matching]
  gwas <- gwas[!is.na(gwas$BP), !names(gwas) %in% "CHRBP", drop = FALSE]

  return(gwas)
}

main(args)
