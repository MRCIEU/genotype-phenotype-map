format_unique_snp_string <- function(chr, bp, ea, oa) {
  # Was going to use compress_alleles but the reference panel doesn't have compressed alleles
  compressed_ea <- trimws(as.character(ea))
  compressed_oa <- trimws(as.character(oa))
  formatted_bp <- format(as.numeric(bp), scientific = FALSE, trim = TRUE)
  snp_string <- glue::glue("{chr}:{formatted_bp}_{compressed_ea}_{compressed_oa}")
  return(snp_string)
}

format_compressed_allele_snp_string <- function(chr, bp, ea, oa) {
  compressed_ea <- compress_alleles(as.character(ea))
  compressed_oa <- compress_alleles(as.character(oa))
  formatted_bp <- trimws(as.character(bp))
  snp_string <- glue::glue("{chr}:{formatted_bp}_{compressed_ea}_{compressed_oa}")
  return(snp_string)
}

compress_alleles <- function(alleles) {
  return(sapply(
    alleles,
    function(allele) if (nchar(allele) > 10) digest::digest(allele, algo = "murmur32") else as.character(allele)
  ))
}


convert_reference_build <- function(study,
                                    file_type,
                                    file,
                                    input_reference_build = reference_builds$GRCh37,
                                    output_reference_build = reference_builds$GRCh38) {
  if (file_type == "vcf") {
    return(convert_vcf_reference_build(study, file, input_reference_build, output_reference_build))
  } else if (file_type == "csv") {
    return(convert_csv_reference_build(study, file, input_reference_build, output_reference_build))
  } else {
    stop(paste(c(
      "Error: file type", file_type, "not recocognised.",
      "File types must be one of:", extraction_file_types
    ), collapse = " "))
  }
}

convert_vcf_reference_build <- function(study, vcf_file) {
  if (input_reference_build == output_reference_build) {
    return(vcf_file)
  }

  dir.create(glue::glue("{study$extracted_location}/vcf"), showWarnings = F, recursive = T)
  liftover_conversion <- available_liftover_conversions[[glue::glue("{input_reference_build}{output_reference_build}")]]
  if (is.null(liftover_conversion)) {
    stop(paste(c(
      "Error: liftOver combination of", input_build, output_build, "not recocognised.",
      "Reference builds must be one of:", reference_builds
    ), collapse = " "))
  }
  output_file <- glue::glue("{study$extracted_location}vcf/hg38.vcf.gz")
  rejected_file <- glue::glue("{study$extracted_location}vcf/hg38_rejected.vcf")
  fasta_file <- glue::glue("{liftover_dir}/hg38.fa")

  if (file.exists(output_file)) {
    message(glue::glue("{vcf_file} already converted to hg38."))
    return(output_file)
  }

  bcf_liftover_command <- glue::glue(
    "/home/bcftools/bcftools annotate --rename-chrs {liftover_dir}/num_to_chr.txt {vcf_file} | ",
    "/home/bcftools/bcftools +liftover --no-version -Ou -- ",
    "-s {liftover_dir}/hg19.fa ",
    "-f {liftover_dir}/hg38.fa ",
    "-c {liftover_conversion} ",
    "--reject {rejected_file} | ",
    "/home/bcftools/bcftools annotate --rename-chrs {liftover_dir}/chr_to_num.txt | ",
    "/home/bcftools/bcftools sort -Oz -o {output_file} -W=tbi"
  )
  system(bcf_liftover_command, wait = T, ignore.stdout = T)

  return(output_file)
}
#' convert_reference_build_via_liftover: Change reference build of BP marker from allow list of liftOver conversions
#'
#' @param gwas: GWAS (file or dataframe) of standardised GWAS
#' @param input_reference_build: string reference build, found in reference_builds list
#' @param output_reference_build: string reference build that GWAS is to change to, found in reference_builds list
#' @return gwas input is altered and returned
#' @export
convert_dataframe_reference_build <- function(gwas,
                                              input_reference_build = reference_builds$GRCh37,
                                              output_reference_build = reference_builds$GRCh38) {
  if (input_reference_build == output_reference_build) {
    return(gwas)
  }

  liftover_conversion <- available_liftover_conversions[[paste0(input_reference_build, output_reference_build)]]
  if (is.null(liftover_conversion)) {
    stop(paste(c(
      "Error: liftOver combination of", input_build, output_build, "not recocognised.",
      "Reference builds must be one of:", reference_builds
    ), collapse = " "))
  }

  original_gwas_size <- nrow(gwas)

  gwas$CHRBP <- paste(gwas$CHR, gwas$BP, sep = ":")
  if (input_reference_build == reference_builds$GRCh37) {
    gwas$BP37 <- gwas$BP
  } else if (input_reference_build == reference_builds$GRCh38) {
    gwas$BP38 <- gwas$BP
  } else if (input_reference_build == reference_builds$GRCh36) {
    gwas$BP36 <- gwas$BP
  }

  bed_file_input <- withr::local_tempfile(fileext = ".bed")
  bed_file_output <- withr::local_tempfile(fileext = ".bed")
  unmapped <- withr::local_tempfile(fileext = ".unmapped")

  create_bed_file_from_gwas(gwas, bed_file_input)
  run_liftover(bed_file_input, bed_file_output, input_reference_build, output_reference_build, unmapped)
  gwas <- use_bed_file_to_update_gwas(gwas, bed_file_output)

  updated_gwas_size <- nrow(gwas)
  if (updated_gwas_size < original_gwas_size) {
    message(paste(
      "Warning: During liftover conversion, the GWAS lost", original_gwas_size - updated_gwas_size,
      "rows out of ", original_gwas_size
    ))
  }

  return(gwas)
}

#' @import tibble
#' @import vroom
create_bed_file_from_gwas <- function(gwas, output_file) {
  bed_format <- tibble::tibble(
    CHR = paste0("chr", gwas$CHR),
    BP1 = gwas$BP,
    BP2 = gwas$BP + 1,
    CHRBP = paste0(CHR, ":", BP1, "-", BP2)
  )

  vroom::vroom_write(bed_format, output_file, col_names = F, delim = " ")
  return(bed_format)
}


run_liftover <- function(bed_file_input, bed_file_output, input_build, output_build, unmapped) {
  lifover_binary <- paste0(liftover_dir, "liftOver")
  liftover_conversion <- available_liftover_conversions[[paste0(input_build, output_build)]]

  liftover_command <- paste(lifover_binary, bed_file_input, liftover_conversion, bed_file_output, unmapped)
  system(liftover_command, wait = T)
  return()
}


use_bed_file_to_update_gwas <- function(gwas, bed_file) {
  liftover_bed <- vroom::vroom(bed_file, col_names = F, show_col_types = F)
  liftover_bed$X1 <- gsub("chr", "", liftover_bed$X1)
  liftover_bed$X4 <- sub(".*:(\\d+)-.*", "\\1", liftover_bed$X4, perl = T)

  bed_map <- tibble::tibble(
    NEW_BP = liftover_bed$X2,
    ORIGINAL_CHRBP = paste(liftover_bed$X1, liftover_bed$X4, sep = ":")
  )

  matching <- match(gwas$CHRBP, bed_map$ORIGINAL_CHRBP)
  gwas$BP <- bed_map$NEW_BP[matching]
  gwas <- gwas[!is.na(gwas$BP), !names(gwas) %in% "CHRBP", drop = FALSE]
  return(gwas)
}

standardise_columns <- function(gwas, N) {
  gwas_columns <- colnames(gwas)

  if (!all(c("CHR", "BP") %in% gwas_columns)) {
    if (all(grepl("\\d:\\d", gwas$SNP))) {
      gwas <- tidyr::separate(data = gwas, col = "SNP", into = c("CHR", "BP"), sep = "[:_]", remove = F)
      gwas$BP <- as.numeric(gwas$BP)
    }
  }

  if (all(c("OR", "OR_LB", "OR_UB") %in% gwas_columns) && !all(c("BETA", "SE") %in% colnames(gwas))) {
    gwas <- convert_or_to_beta(gwas)
  }

  if ("LOG_P" %in% gwas_columns && !"P" %in% gwas_columns) {
    gwas <- convert_negative_log_p_to_p(gwas)
  }

  if ("Z" %in% gwas_columns && !"BETA" %in% gwas_columns) {
    gwas <- convert_z_score_to_beta(gwas)
  }

  if ("BP" %in% gwas_columns) gwas$BP <- as.numeric(gwas$BP)
  if ("P" %in% gwas_columns) {
    gwas$P <- as.numeric(gwas$P)
    gwas$P[gwas$P == 0] <- .Machine$double.xmin
  }
  if ("BETA" %in% gwas_columns) {
    gwas$BETA <- as.numeric(gwas$BETA)
  }

  if ("BETA" %in% gwas_columns) {
    gwas$BETA <- as.numeric(gwas$BETA)
  }

  return(gwas)
}
standardise_alleles <- function(gwas) {
  gwas$EA <- toupper(gwas$EA)
  gwas$OA <- toupper(gwas$OA)

  to_flip <- (gwas$EA > gwas$OA) & (!gwas$EA %in% c("D", "I"))
  if (any(to_flip)) {
    gwas$EAF[to_flip] <- 1 - gwas$EAF[to_flip]
    gwas$BETA[to_flip] <- -1 * gwas$BETA[to_flip]

    temp <- gwas$OA[to_flip]
    gwas$OA[to_flip] <- gwas$EA[to_flip]
    gwas$EA[to_flip] <- temp
  }

  gwas$SNP <- toupper(paste0(gwas$CHR, ":", format(gwas$BP, scientific = F, trim = T), "_", gwas$EA, "_", gwas$OA))
  gwas <- dplyr::select(gwas, SNP, CHR, BP, EA, OA, dplyr::everything())

  return(gwas)
}

#' change_column_names: private function that takes a named list of column names
#'  and changes the supplied data frame's column names accordingly
#'
#' @param: gwas dataframe of gwas to standardise column names
#' @param: columns named list for
#' @param: opposite_mapping logical flag on if we are mapping from key to value or vice verca
change_column_names <- function(gwas, columns = list(), remove_extra_columns = F) {
  if (!is.null(columns$N) && is.numeric(columns$N)) {
    gwas$N <- columns$N
    # TODO: remove N from columns
  }
  for (name in names(columns)) {
    # Skip if the mapping value is empty, NULL, or zero-length
    if (is.null(columns[[name]]) || length(columns[[name]]) == 0 ||
          (is.character(columns[[name]]) && nchar(columns[[name]]) == 0)) {
      next
    }

    # this deletes an existing column that we're about to rename, so we don't have 2 columns
    already <- name != columns[[name]]
    already <- length(already) > 0 && any(already, na.rm = TRUE)
    if (name %in% names(gwas) && already) {
      gwas <- gwas[, -which(names(gwas) %in% c(name))]
    }
    names(gwas)[names(gwas) == columns[name]] <- name
  }

  if (remove_extra_columns) {
    columns_to_remove <- setdiff(colnames(gwas), names(columns))
    gwas <- gwas[, -which(colnames(gwas) %in% columns_to_remove)]
  }

  return(gwas)
}

split_into_regions <- function(gwas, ld_blocks, study_metadata, p_value_threshold) {
  ld_blocks_prepared <- ld_blocks |>
    dplyr::mutate(
      ld_block = ld_block_string(study_metadata$ancestry, chr, start, stop),
      chr = as.numeric(chr),
      start = as.numeric(start),
      stop = as.numeric(stop)
    )

  gwas <- data.table::as.data.table(gwas)
  ld_blocks <- data.table::as.data.table(ld_blocks_prepared)

  data.table::setkey(ld_blocks, chr, start, stop)
  gwas[, `:=`(start = BP, end = BP)]
  data.table::setkey(gwas, CHR, start, end)

  gwas_with_blocks <- data.table::foverlaps(
    gwas,
    ld_blocks,
    by.x = c("CHR", "start", "end"),
    by.y = c("chr", "start", "stop"),
    type = "within",
    nomatch = NULL
  )
  gwas_with_blocks <- gwas_with_blocks[!is.na(ld_block)]
  gwas_with_blocks[, `:=`(start = NULL, end = NULL)]
  gwas_with_blocks <- tibble::as_tibble(gwas_with_blocks) |>
    dplyr::select(-dplyr::any_of(c("chr", "start", "stop", "ancestry")))
  gwas_split <- split(gwas_with_blocks, gwas_with_blocks$ld_block)

  extracted_regions <- parallel::mclapply(gwas_split, mc.cores = 10, function(snps_in_block) {
    top_hit <- dplyr::arrange(snps_in_block, P) |>
      dplyr::slice(1)

    variant_p <- as.numeric(top_hit$P)
    variant_bp <- as.numeric(top_hit$BP)
    variant_chr <- as.numeric(top_hit$CHR)

    if (variant_p > p_value_threshold) {
      return(NULL)
    }

    extracted_file <- glue::glue(
      "{extracted_study_dir}/extracted/{study_metadata$ancestry}_{variant_chr}_{variant_bp}.tsv.gz"
    )
    vroom::vroom_write(snps_in_block, extracted_file)

    extraction_info <- data.frame(
      chr = variant_chr,
      bp = variant_bp,
      log_p = -log10(variant_p),
      ld_block = snps_in_block$ld_block[1],
      file = extracted_file,
      cis_trans = NA
    )

    return(extraction_info)
  })

  extracted_regions <- extracted_regions[!sapply(extracted_regions, is.null)] |>
    dplyr::bind_rows()

  message(glue::glue("Extracted {nrow(extracted_regions)} regions from {study_metadata$study_name}"))

  if (nrow(extracted_regions) > 0) {
    extracted_regions <- extracted_regions |>
      dplyr::arrange(dplyr::desc(log_p))
    return(extracted_regions)
  } else {
    return(data.frame())
  }
}

#' convert_or_to_beta: Given an OR and lower and upper bounds,
#'   calculates the BETA, and SE.
#'   based on this answer: https://stats.stackexchange.com/a/327684
#'
#' @param gwas: dataframe with the following columns: OR, LB (lower bound), UB (upper bound)
#' @return gwas with new columns BETA and SE
#' @import stats
#' @export
convert_or_to_beta <- function(gwas) {
  gwas <- get_file_or_dataframe(gwas)
  if (!all(c("OR", "OR_LB", "OR_UB") %in% colnames(gwas))) {
    stop("Need OR, OR_LB + OR_UB to complete conversion")
  }

  z_score <- stats::qnorm(.975, mean = 0, sd = 1) #1.96
  gwas$BETA <- log(gwas$OR)
  gwas$SE <- (log(gwas$OR_LB) - gwas$BETA) / -z_score

  return(gwas)
}