source('constants.R')
#TODO: 3 step process per 'study' (<data_source>-<specifier>-<gene>)
# 1. Find cis snps (via cis-wind)
# 2. Find all other significant hits, which will be considered trans
#  with these 2 steps, I can create extracted_snps.tsv, but with cis_trans filled out appropriately
# 3. Using extracted_snps, extract whole regions and store

parser <- argparser::arg_parser('Extract genomic regions from a BESD file')
parser <- argparser::add_argument(parser, '--extracted_study_location', help = 'Study location', type = 'character')
parser <- argparser::add_argument(parser, '--bespoke_parsing', help = 'Flag to signal bespoke parsing of besd format', type = 'character', default = 'none')
parser <- argparser::add_argument(parser, '--extracted_output_file', help = 'Extracted Output file', type = 'character')
args <- argparser::parse_args(parser)

ld_blocks <- vroom::vroom('data/ld_blocks.tsv', show_col_types = F)

main <- function() {
  study <- vroom::vroom(glue::glue('{pipeline_metadata_dir}/studies_to_process.tsv'), show_col_types = F) |>
    dplyr::filter(extracted_location == args$extracted_study_location)
  if (nrow(study) != 1) stop('Error: cant find study to process')

  if (study$reference_build != reference_builds$GRCh38) {
    stop(glue::glue('Error: Only BESD files using {reference_builds$GRCh38} is allowed right now'))
  }

  p_value_threshold <- ifelse(is.na(study$p_value_threshold), lowest_p_value_threshold, study$p_value_threshold)
  metadata <- jsonlite::fromJSON(glue::glue('{study$study_location}.json'))

  extracted_cis_snps <- extracted_trans_snps <- data.frame(chr=character(),
                                                           bp=numeric(),
                                                           log_p=numeric(),
                                                           ld_block=character(),
                                                           file=character(),
                                                           cis_trans=character()
  )
  cis_results <- trans_results <- list()
  if (metadata$cis_trans == cis_trans$cis_only || metadata$cis_trans == cis_trans$cis_trans) {
    cis_results <- extract_cis_region(study, p_value_threshold)
  }

  if (metadata$cis_trans == cis_trans$trans_only ||
     (metadata$cis_trans == cis_trans$cis_trans && !is.null(cis_results))) {
    trans_results <- extract_trans_regions(cis_results$snp_data, study, p_value_threshold)
  }

  clumped_snps <- dplyr::bind_rows(cis_results$clumped_snps, trans_results$clumped_snps)
  vroom::vroom_write(clumped_snps, glue::glue('{args$extracted_study_location}/clumped_snps.tsv'))

  extracted_snps <- dplyr::bind_rows(cis_results$snp_data, trans_results$snp_data)
  vroom::vroom_write(extracted_snps, args$extracted_output_file)
}

extract_cis_region <- function(study, p_value_threshold) {
  tmp_cis_region <- glue::glue('/tmp/{study$study_name}_top_snp')
  extract_cis_region <- glue::glue('smr --beqtl-summary {study$study_location} ',
                           '--query 1 ',
                           '--probe {study$probe} ',
                           '--cis-wind 2000 ',
                           '--out {tmp_cis_region}'
  )

  system(extract_cis_region, wait=T, ignore.stdout = T)
  if (!file.exists(glue::glue('{tmp_cis_region}.txt'))) return()
  cis_region <- vroom::vroom(glue::glue('{tmp_cis_region}.txt'), show_col_types = F) |>
    format_gwas()
  if (nrow(cis_region) < 1) return()

  top_hit <- cis_region[which.min(cis_region$P), ]
  if (top_hit$P > p_value_threshold) return()

  ld_block <- dplyr::filter(ld_blocks, chr == top_hit$CHR & start <= top_hit$BP & stop > top_hit$BP & ancestry == study$ancestry)
  ld_block_string <- ld_block_string(ld_block$ancestry, ld_block$chr, ld_block$start, ld_block$stop)

  if (nrow(ld_block) == 0) {
    missing <- data.frame(study=study$study_name, chr=top_hit$CHR, bp=top_hit$BP)
    vroom::vroom_write(missing, glue::glue('{pipeline_metadata_dir}/missing_ld_blocks.tsv'), append = T)
    message('Missing LD block for ', top_hit$SNP)
    return()
  }

  dir.create(glue::glue('{study$extracted_location}/extracted'), showWarnings = F, recursive = T)
  dir.create(glue::glue('{study$extracted_location}/standardised'), showWarnings = F, recursive = T)
  dir.create(glue::glue('{study$extracted_location}/imputed'), showWarnings = F, recursive = T)
  dir.create(glue::glue('{study$extracted_location}/finemapped'), showWarnings = F, recursive = T)

  cis_region <- dplyr::filter(cis_region, BP >= ld_block$start & BP <= ld_block$stop) 
  extracted_file <- glue::glue('{study$extracted_location}extracted/{study$ancestry}_{top_hit$CHR}_{top_hit$BP}.tsv.gz')
  vroom::vroom_write(cis_region, extracted_file)

  unlink(glue::glue('{tmp_cis_region}.txt'))

  extracted_snps <- data.frame(chr = as.character(top_hit$CHR),
                               bp = top_hit$BP,
                               log_p = -log10(top_hit$P),
                               ld_block = ld_block_string,
                               file = extracted_file,
                               cis_trans = 'cis'
  )
  message(glue::glue('found {nrow(cis_region)} cis snps for {study$study_name}'))
  return(list(snp_data=extracted_snps, clumped_snps=data.frame()))
}

#' extract_trans_regions
#' @param extracted_cis_snps: region of extracted cis snps, to be filtered out of trans results
#' @param study of interest to extract
#' @param p_value_threshold to extract
#' 1: get all hits above p-value threshold
#' 2: clump results, filter out cis region (if any)
#' 3: loop through clumped results to get all regions, then extract and store
extract_trans_regions <- function(extracted_cis_snp, study, p_value_threshold) {
  tmp_trans_snps <- glue::glue('/tmp/{study$study_name}_top_trans_snps')

  extract_top_snps <- glue::glue('smr --beqtl-summary {study$study_location} ',
                           '--query {p_value_threshold} ',
                           '--probe {study$probe} ',
                           '--out {tmp_trans_snps}'
  )
  system(extract_top_snps, wait=T, ignore.stdout = T)
  if (!file.exists(glue::glue('{tmp_trans_snps}.txt'))) return()
  probe_top_hits <- vroom::vroom(glue::glue('{tmp_trans_snps}.txt'), show_col_types = F)

  full_bfile <- glue::glue('{ld_reference_panel_dir}/{study$ancestry}/full')
  plink_command <- glue::glue('plink2 --bfile {full_bfile} ', 
                          '--clump {tmp_trans_snps}.txt ',
                          '--clump-snp-field SNP ',
                          '--clump-field p ',
                          '--clump-p1 {p_value_threshold} ',
                          '--clump-kb 1000 ',
                          '--clump-r2 0.1 ',
                          '--out {tmp_trans_snps}'
  )
  system(plink_command, wait=T, ignore.stdout = T)

  clumped_trans_snps <- data.table::fread(glue::glue('{tmp_trans_snps}.clumps')) |>
    dplyr::rename(SNP='ID', CHR='#CHROM', BP='POS') |>
    dplyr::select(SNP, CHR, BP, P) |>
    dplyr::mutate(CHR=as.numeric(CHR), BP=as.numeric(BP), P=as.numeric(P)) |>
    dplyr::arrange(P)
  if (nrow(clumped_trans_snps) < 1) return()

  ld_block_strings <- apply(clumped_trans_snps, 1, function(clump) {
    bp <- as.numeric(clump['BP'])
    ld_block <- dplyr::filter(ld_blocks, chr == as.numeric(clump['CHR']) & start <= bp & stop > bp & ancestry == study$ancestry)
    return(ld_block_string(ld_block$ancestry, ld_block$chr, ld_block$start, ld_block$stop))
  })
  clumped_trans_snps$ld_block_string <- ld_block_strings

  # removing duplicate entries per region, and the original cis region so we only grab the correct retions once.
  clumped_trans_snps <- clumped_trans_snps[!duplicated(clumped_trans_snps$ld_block_string), ]
  if (is.data.frame(extracted_cis_snp) && nrow(extracted_cis_snp) == 1) {
    clumped_trans_snps <- dplyr::filter(clumped_trans_snps, 
      !is.na(ld_block_string) & ld_block_string != extracted_cis_snp$ld_block
    )
  }

  message(glue::glue('found {nrow(clumped_trans_snps)} new regions to extract for trans results'))

  extracted_trans_snps <- apply(clumped_trans_snps, 1, function(clumped_snp) {
    trans_bp <- as.numeric(clumped_snp['BP'])
    trans_chr <- as.numeric(clumped_snp['CHR'])
    trans_p <- as.numeric(clumped_snp['P'])
    ld_block <- dplyr::filter(ld_blocks, chr == trans_chr & start <= trans_bp & stop > trans_bp & ancestry == study$ancestry)

    tmp_trans_region <- glue::glue('/tmp/{study$study_name}_trans_region')
    extract_region <- glue::glue('smr --beqtl-summary {study$study_location} ',
                            '--query 1 ',
                            '--snp {clumped_snp["SNP"]} ',
                            '--snp-wind 3000 ', # smr doesn't accept BP ranges, and errors if specific RSID isn't present
                            '--probe {study$probe} ',
                            '--out {tmp_trans_region}'
    )
    system(extract_region, wait=T, ignore.stdout = T)
    if (!file.exists(glue::glue('{tmp_trans_region}.txt'))) return()

    trans_region <- vroom::vroom(glue::glue('{tmp_trans_region}.txt'), show_col_types = F)
    trans_region <- format_gwas(trans_region) |>
      dplyr::filter(BP >= ld_block$start & BP <= ld_block$stop) 

    extracted_file <- glue::glue('{study$extracted_location}extracted/{study$ancestry}_{trans_chr}_{trans_bp}.tsv.gz')
    vroom::vroom_write(trans_region, extracted_file)

    # ld_block_string <- ld_block_string(ld_block$ancestry, ld_block$chr, ld_block$start, ld_block$stop)

    # if (nrow(ld_block) == 0) {
      # missing <- data.frame(study=study$study_name, chr=trans_chr, bp=trans_bp)
      # vroom::vroom_write(missing, glue::glue('{pipeline_metadata_dir}/missing_ld_blocks.tsv'), append = T)
      # message('Missing LD block for ', clumped_snp['SNP'])
      # return()
    # }

    unlink(glue::glue('{tmp_trans_region}.txt'))

    extracted_trans_hit <- data.frame(chr = as.character(trans_chr),
                                    bp = trans_bp,
                                    log_p = -log10(trans_p),
                                    ld_block = clumped_snp['ld_block_string'],
                                    file = extracted_file,
                                    cis_trans = 'trans'
    )
    return(extracted_trans_hit)
  }) |> dplyr::bind_rows()

  unlink(glue::glue('{tmp_trans_snps}.txt'))

  return(list(snp_data=extracted_trans_snps, clumped_snps=clumped_trans_snps))
}

format_gwas <- function(gwas) {
  gwas <- dplyr::rename(gwas, RSID='SNP', CHR='Chr', EA='A1', OA='A2', EAF='Freq', BETA='b', P='p') |>
    dplyr::select(-Probe, -Probe_Chr, -Probe_bp, -Gene, -Orientation)

  return(gwas)
}

main()
