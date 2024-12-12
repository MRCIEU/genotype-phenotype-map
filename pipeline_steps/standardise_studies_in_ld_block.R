source('constants.R')

parser <- argparser::arg_parser('Standardise GWAS for pipeline')
parser <- argparser::add_argument(parser, '--ld_block', help = 'LD block that the ', type = 'character')
parser <- argparser::add_argument(parser, '--completed_output_file', help = 'Completed output file', type = 'character')
args <- argparser::parse_args(parser)

main <- function() {
  ld_info <- ld_block_dirs(args$ld_block)
  ld_matrix_info <- vroom::vroom(glue::glue('{ld_info$ld_reference_panel_prefix}.tsv'), show_col_types = F)

  extracted_studies_file <- glue::glue('{ld_info$ld_block_data}/extracted_studies.tsv')
  extracted_studies  <- vroom::vroom(extracted_studies_file , show_col_types = F)

  standardised_studies_file <- glue::glue('{ld_info$ld_block_data}/standardised_studies.tsv')
  if (file.exists(standardised_studies_file)) {
    existing_standardised_studies <- vroom::vroom(standardised_studies_file,
                                                  show_col_types = F,
                                                  col_types = standardised_column_types 
    )
  } else {
    existing_standardised_studies <- empty_standardised_studies()
  }

  if (nrow(extracted_studies) > 0) {
    standardised_studies <- apply(extracted_studies, 1, function (study) {
      start_time <- Sys.time()
      standardised_file <- sub('extracted', 'standardised', study[['file']])

      if (standardised_file %in% existing_standardised_studies$file) {
        return()
      }

      result <- perform_standardisation(study, ld_matrix_info)

      if (nrow(result$gwas) < 100) {
        return()
      }
      vroom::vroom_write(result$gwas, result$study$file)

      result$study$time_taken <- hms::as_hms(difftime(Sys.time(), start_time)) 
      return(result$study)
    }) |>
      dplyr::bind_rows() |>
      type.convert(as.is=T)

    standardised_studies$chr <- as.character(standardised_studies$chr)
  }

  if (nrow(standardised_studies) > 0) {
    standardised_studies <- dplyr::bind_rows(existing_standardised_studies, standardised_studies) |>
      dplyr::distinct()

    vroom::vroom_write(standardised_studies, standardised_studies_file)
  }

  vroom::vroom_write(data.frame(), args$completed_output_file)
}

perform_standardisation <- function(study, ld_matrix_info) {
  standardised_file <- sub('extracted', 'standardised', study[['file']])
  gwas <- vroom::vroom(study[['file']], show_col_types = F)

  response <- standardise_alleles(gwas) |>
    standardise_extracted_gwas(ld_matrix_info)

  study['ld_block'] <- args$ld_block
  study['file'] <- standardised_file
  study['eaf_from_reference_panel'] <- response$eaf_from_reference_panel
  study['snps_removed_by_reference_panel'] <- response$snps_removed_by_reference_panel
  study <- study[-match('reference_build', names(study))]

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
                    time_taken=character()
  ))
}

standardise_extracted_gwas <- function(gwas, ld_matrix_info) {
  eaf_from_reference_panel <- FALSE
  original_gwas_size <- nrow(gwas)
  gwas <- dplyr::distinct(gwas, CHR, BP, EA, OA, .keep_all = TRUE)
  
  if (!"Z" %in% colnames(gwas)) {
    gwas <- dplyr::mutate(gwas, SE = replace(SE, SE == 0, 0.00001), Z = BETA / SE)
  }
  
  if (!"P" %in% colnames(gwas) && "LP" %in% colnames(gwas)) {
    gwas$LP <- as.numeric(gwas$LP)
    gwas <- dplyr::mutate(gwas, P = 10^(-LP)) |>
      dplyr::select(-LP)
  }

  gwas <- dplyr::filter(gwas, SNP %in% ld_matrix_info$SNP)
  ld_matrix_info <- dplyr::filter(ld_matrix_info, SNP %in% gwas$SNP)

  if (all(is.na(gwas$EAF))) {
    gwas <- dplyr::select(gwas, -EAF) |>
      dplyr::left_join(ld_matrix_info |> dplyr::select(SNP, EAF), by = 'SNP')
    eaf_from_reference_panel <- TRUE
  }

  columns_to_coerce <- c("EAF") # Add BETA and SE if needed
  gwas <- tidyr::drop_na(gwas, dplyr::all_of(columns_to_coerce)) |>
    dplyr::arrange(match(SNP, ld_matrix_info$SNP))
  
  return(list(gwas = gwas,
              eaf_from_reference_panel = eaf_from_reference_panel,
              snps_removed_by_reference_panel = original_gwas_size - nrow(gwas)))
}

standardise_alleles <- function(gwas) {
  gwas$EA <- toupper(gwas$EA)
  gwas$OA <- toupper(gwas$OA)
  columns_to_coerce <- c("EAF", "BETA", "SE")
  gwas <- dplyr::mutate(gwas, dplyr::across(dplyr::all_of(columns_to_coerce), as.numeric))

  to_flip <- (gwas$EA > gwas$OA) & (!gwas$EA %in% c("D", "I"))
  if (any(to_flip)) {
    if ('EAF' %in% names(gwas)) {
      gwas$EAF[to_flip] <- 1 - gwas$EAF[to_flip]
    }
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

  compressed_ea <- compress_alleles(gwas$EA)
  compressed_oa <- compress_alleles(gwas$OA)
  formatted_bp <- format(gwas$BP, scientific = F, trim = T)

  gwas$SNP <- glue::glue('{gwas$CHR}:{formatted_bp}_{compressed_ea}_{compressed_oa}')
  return(gwas)
}

compress_alleles <- function(alleles) {
  sapply(alleles, function(allele) if(nchar(allele) > 10) digest::digest(allele, algo='murmur32') else allele)
}

# TODO: if we need this in the future, maybe move it to the extract step
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

  chain_file <- available_liftover_conversions[[glue::glue('{input_reference_build}{output_reference_build}')]]
  if (is.null(chain_file)) {
    stop(paste(c("Error: liftOver combination of", input_build, output_build, "not recocognised.",
                 "Reference builds must be one of:", reference_builds), collapse = " "))
  }

  original_gwas_size <- nrow(gwas)

  bed_file_input <- tempfile(fileext = ".bed")
  bed_file_output <- tempfile(fileext = ".bed")
  unmapped <- tempfile(fileext = ".unmapped")

  create_bed_file_from_gwas(gwas, bed_file_input)
  run_liftover(bed_file_input, bed_file_output, chain_file, unmapped)
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
    CHR = glue::glue('chr{gwas$CHR}'),
    BP1 = gwas$BP,
    BP2 = gwas$BP+1,
    CHRBP = glue::glue('{CHR}:{BP1}-{BP2}')
  )

  vroom::vroom_write(bed_format, output_file, col_names=F, delim=" ")
  return(bed_format)
}

run_liftover <- function(bed_file_input, bed_file_output, chain_file, unmapped) {
  lifover_binary <- glue::glue('{liftover_dir}liftOver')

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

main()
