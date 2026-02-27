source("constants.R")

# Input study_location and extracted_location from study metadata file
parser <- argparser::arg_parser("Extract genomic regions from sequencing GWAS summary statistics")
parser <- argparser::add_argument(
  parser,
  "--extracted_study_location",
  help = "Summary stats file",
  type = "character"
)
parser <- argparser::add_argument(
  parser,
  "--bespoke_parsing",
  "Flag to signal bespoke parsing of besd format",
  type = "character",
  default = "none"
)
parser <- argparser::add_argument(
  parser,
  "--extracted_output_file",
  help = "Extracted output location",
  type = "character"
)
args <- argparser::parse_args(parser)

ld_blocks <- vroom::vroom("data/ld_blocks.tsv", show_col_types = F)
ld_blocks$ld_block_string <- ld_block_string(ld_blocks$ancestry, ld_blocks$chr, ld_blocks$start, ld_blocks$stop)

# Lower MAF theshold (for 10 expected observations of minor allele in 431,000 Europeans (approx n UKB))
eaf_min_threshold <- 0.000012
eaf_max_threshold <- 0.01
required_columns <- c("CHR", "BP", "EA", "OA", "EAF", "BETA", "SE", "P")
beneficial_columns <- c("GENE", "ANNOTATION")

#' extract_regions_from_rare_tsv:
#' 1. Identify GRCh38 ld regions containing rare variant top hits passing threshold (no clumping performed here)
#' 2. Extract SNPs across region and store
main <- function() {
  study <- vroom::vroom(glue::glue("{pipeline_metadata_dir}/studies_to_process.tsv"), show_col_types = F) |>
    dplyr::filter(extracted_location == args$extracted_study_location)
  if (nrow(study) != 1) stop("Error: cant find study to process")

  p_value_threshold <- ifelse(is.na(study$p_value_threshold), lowest_p_value_threshold, study$p_value_threshold)

  dir.create(glue::glue("{study$extracted_location}/svgs"), showWarnings = F, recursive = T)
  dir.create(glue::glue("{study$extracted_location}/extracted"), showWarnings = F, recursive = T)
  dir.create(glue::glue("{study$extracted_location}/standardised"), showWarnings = F, recursive = T)
  dir.create(glue::glue("{study$extracted_location}/imputed"), showWarnings = F, recursive = T)
  dir.create(glue::glue("{study$extracted_location}/finemapped"), showWarnings = F, recursive = T)

  gwas <- vroom::vroom(study$study_location, show_col_types = F)
  gwas <- check_gwas(gwas)

  gwas <- filter_snps(gwas)
  if (nrow(gwas) == 0) {
    empty_extracted_snps <- data.frame(
      chr = character(),
      bp = numeric(),
      log_p = numeric(),
      ld_block = character(),
      file = character(),
      cis_trans = character()
    )
    vroom::vroom_write(empty_extracted_snps, args$extracted_output_file)
  } else {
    # Identify regions top variants fall in and write to file
    vars_regions <- split_into_regions(gwas, study, p_value_threshold)
    vroom::vroom_write(vars_regions, args$extracted_output_file)
  }
  return()
}

check_gwas <- function(gwas) {
  study_cols <- colnames(gwas)

  if (all(required_columns %in% study_cols)) {
    return(gwas)
  }
  if (all(c("OR", "CI_UPPER", "CI_LOWER") %in% study_cols) && all(required_columns[c(-6, -7)] %in% study_cols)) {
    gwas$OR <- ifelse(gwas$OR == 0, 0.01, gwas$OR)
    gwas$BETA <- log(gwas$OR)

    UCI <- log(gwas$CI_UPPER)
    gwas$SE <- (UCI - gwas$BETA) / 1.96
    gwas$SE <- ifelse(gwas$OR == 0.01, 0, gwas$SE)

    return(gwas)
  } else {
    stop("Study must have BETA and SE or OR and CI_UPPER/CI_LOWER plus CHR, BP, EA, OA, EAF, P")
  }
}

filter_snps <- function(gwas) {
  study_cols <- colnames(gwas)

  if (!all(c("CHR", "P", "EAF") %in% study_cols)) {
    stop("Missing CHR, P or EAF column")
  }

  gwas_filt <- gwas |>
    dplyr::filter(
      CHR %in% seq(1, 22),
      ((EAF <= eaf_max_threshold & EAF >= eaf_min_threshold) |
         (1 - EAF <= eaf_max_threshold & 1 - EAF >= eaf_min_threshold))
    ) |>
    dplyr::arrange(CHR, BP)
  return(gwas_filt)
}

split_into_regions <- function(gwas, study, p_value_threshold) {
  gwas <- dplyr::select(gwas, dplyr::any_of(c(required_columns, beneficial_columns)))

  ld_block_strings <- apply(gwas, 1, function(extracted) {
    bp <- as.numeric(extracted["BP"])
    ld_block <- dplyr::filter(
      ld_blocks,
      chr == as.numeric(extracted["CHR"]) & start <= bp & stop > bp & ancestry == study$ancestry
    )
    return(ld_block_string(ld_block$ancestry, ld_block$chr, ld_block$start, ld_block$stop))
  })
  gwas$ld_block_string <- ld_block_strings

  extracted_regions <- lapply(unique(ld_block_strings), function(ld_block_identifier) {
    snps_in_block <- dplyr::filter(gwas, ld_block_string == ld_block_identifier) |>
      dplyr::select(-ld_block_string)

    top_hit <- dplyr::arrange(snps_in_block, P) |>
      dplyr::slice(1)

    variant_p <- as.numeric(top_hit["P"])
    variant_bp <- as.numeric(top_hit["BP"])
    variant_chr <- as.numeric(top_hit["CHR"])

    if (variant_p > p_value_threshold) {
      return()
    }

    extracted_file <- glue::glue(
      "{study$extracted_location}extracted/{study$ancestry}_{variant_chr}_{variant_bp}.tsv.gz"
    )
    vroom::vroom_write(snps_in_block, extracted_file)

    extraction_info <- data.frame(
      chr = variant_chr,
      bp = variant_bp,
      log_p = -log10(as.numeric(top_hit["P"])),
      ld_block = ld_block_identifier,
      file = extracted_file,
      cis_trans = NA
    )
    return(extraction_info)
  }) |>
    dplyr::bind_rows()

  message(glue::glue("Extracted {nrow(extracted_regions)} regions from {study$study_name}"))
  return(extracted_regions)
}

main()
