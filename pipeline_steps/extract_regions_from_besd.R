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
  if (metadata$cis_trans == cis_trans$trans_only || metadata$cis_trans == cis_trans$cis_trans) {
    trans_results <- extract_trans_regions(extracted_cis_snps, study, p_value_threshold)
  }

  #TODO: save clumped_snps here
  # clumped_snps <- dplyr::bind_rows(cis_results$clumped_snps, trans_results$clumped_snps)
  # vroom::vroom_write(clumped_snps, glue::glue('{extracted_study_location}/clumped_snps.tsv'))

  extracted_snps <- dplyr::bind_rows(cis_results$snp_data, trans_results$snp_data)
  vroom::vroom_write(extracted_snps, args$extracted_output_file)
}


extract_cis_region <- function(study, p_value_threshold) {
  tmp_cis_snp <- glue::glue('/tmp/{study$study_name}_top_snp')
  extract_top_snp <- glue::glue('smr --beqtl-summary {study$study_location} ',
                           '--query {p_value_threshold} ',
                           '--probe {study$probe} ',
                           '--cis-wind 1 ',
                           '--out {tmp_cis_snp}'
  )

  system(extract_top_snp, wait=T, ignore.stdout = T)
  if (!file.exists(glue::glue('{tmp_cis_snp}.txt'))) return()
  top_cis_snp <- vroom::vroom(glue::glue('{tmp_cis_snp}.txt'), show_col_types = F)
  if (nrow(top_cis_snp) < 1) return()

  dir.create(glue::glue('{study$extracted_location}/extracted'), showWarnings = F, recursive = T)
  dir.create(glue::glue('{study$extracted_location}/standardised'), showWarnings = F, recursive = T)
  dir.create(glue::glue('{study$extracted_location}/imputed'), showWarnings = F, recursive = T)
  dir.create(glue::glue('{study$extracted_location}/finemapped'), showWarnings = F, recursive = T)

  top_cis_snp <- top_cis_snp[top_cis_snp$p == min(top_cis_snp$p), ][1, ]
  ld_block <- dplyr::filter(ld_blocks, chr == top_cis_snp$Chr & start <= top_cis_snp$BP & stop > top_cis_snp$BP & ancestry == study$ancestry)
  ld_block_string <- ld_block_string(ld_block$ancestry, ld_block$chr, ld_block$start, ld_block$stop)

  if (nrow(ld_block) == 0) {
    missing <- data.frame(study=study$study_name, chr=top_cis_snp$Chr, bp=top_cis_snp$BP)
    vroom::vroom_write(missing, glue::glue('{pipeline_metadata_dir}/missing_ld_blocks.tsv'), append = T)
    message('Missing LD block for ', top_cis_snp$SNP)
  }

  tmp_cis_region <- glue::glue('/tmp/{study$study_name}_cis_region')
  extract_region <- paste('smr --beqtl-summary', study$study_location,
                              '--query 1',
                              '--snp', top_cis_snp$SNP,
                              '--snp-wind 1000', # smr doesn't accept BP ranges, and errors if specific RSID isn't present
                              '--probe ', study$probe,
                              '--out ', tmp_cis_region
  )
  system(extract_region, wait=T, ignore.stdout = T)

  cis_region <- vroom::vroom(glue::glue('{tmp_cis_region}.txt'), show_col_types = F)
  cis_region <- format_gwas(cis_region) |>
    dplyr::filter(BP >= ld_block$start & BP <= ld_block$stop) |>
    gwas_health_check() |>
    filter_gwas() 

  extracted_file <- glue::glue('{study$extracted_location}extracted/{study$ancestry}_{top_cis_snp$Chr}_{top_cis_snp$BP}.tsv.gz')
  vroom::vroom_write(cis_region, extracted_file)

  extracted_snps <- data.frame(chr = as.character(top_cis_snp$Chr),
                               bp = top_cis_snp$BP,
                               log_p = -log10(top_cis_snp$p),
                               ld_block = ld_block_string,
                               file = extracted_file,
                               cis_trans = 'cis'
  )
  return(list(snp_data=extracted_snps, clumped_snps=NULL))
}

# TODO: WARNING - UNTESTED 
#' extract_trans_regions
#' @param extracted_cis_snps: region of extracted cis snps, to be filtered out of trans results
#' 1: get all hits above p-value threshold
#' 2: clump results, somehow filter out cis region
#' 3: loop through clumped results to get all regions
extract_trans_regions <- function(extracted_cis_snps, study, p_value_threshold) {
  tmp_trans_snps <- glue::glue('/tmp/{study$study_name}_top_trans_snps')
  extract_top_snps <- glue::glue('smr --beqtl-summary {study$study_location} ',
                           '--query {p_value_threshold} ',
                           '--probe {study$probe} ',
                           '--cis-wind 1 ',
                           '--out {tmp_trans_snps}'
  )
  system(extract_top_snps, wait=T, ignore.stdout = T)
  if (!file.exists(glue::glue('{tmp_trans_snps}.txt'))) return()
  probe_top_hits <- vroom::vroom(glue::glue('{tmp_smr_result}.txt'), show_col_types = F)

  #filter out all cis snps (if there are any)
  if (nrow(extracted_cis_snps) > 0) {
    min_bp <- min(extracted_cis_snps$bp)
    max_bp <- max(extracted_cis_snps$bp)

    probe_top_hits <- dplyr::filter(probe_top_hits, !(Chr == extracted_cis_snps$Chr[1] & BP > min_bp & BP < max_bp))
  }

  full_bfile <- glue::glue('{ld_reference_panel_dir}/{study$ancestry}/full')
  plink_command <- glue::glue('plink1.9 --bfile {full_bfile} ', 
                          '--clump {probe_top_hits}.txt ',
                          '--clump-snp-field SNP ',
                          '--clump-p1 {p_value_threshold} ',
                          '--clump-kb 1000 ',
                          '--clump-r2 0.1 ',
                          '--out {tmp_trans_snps}'
  )

  system(plink_command, wait=T, ignore.stdout = T)

  clumped_trans_snps <- data.table::fread(glue::glue('{tmp_trans_snps}.clumped'))
  if (nrow(clumped_trans_snps) < 1) return()

  extracted_trans_snps <- apply(clumped_trans_snps, 1, function(clumped_snp) {
    tmp_trans_region <- glue::glue('/tmp/{study$study_name}_cis_region')
    extract_region <- glue::glue('smr --beqtl-summary {study$study_location} ',
                            '--query 1 ',
                            '--snp {clumped_snp["SNP"]} ',
                            '--snp-wind 1000 ', # smr doesn't accept BP ranges, and errors if specific RSID isn't present
                            '--probe {study$probe} ',
                            '--out {tmp_trans_region}'
    )
    system(extract_region, wait=T, ignore.stdout = T)
    trans_region <- vroom::vroom(glue::glue('{tmp_trans_region}.txt'), show_col_types = F)
    trans_region <- format_gwas(trans_region) |>
      gwas_health_check() |>
      filter_gwas()

    extracted_file <- glue::glue('{study$extracted_location}extracted/{study$ancestry}_{clumped_snp["CHR"]}_{clumped_snp["BP"]}.tsv.gz')
    vroom::vroom_write(trans_region, extracted_file)
    ld_block <- dplyr::filter(ld_blocks, chr == clumped_snp['CHR'] & start <= clumped_snp['BP'] & stop > clumped_snp['BP'] & ancestry == study$ancestry)

    ld_block_string <- ld_block_string(ld_block$ancestry, ld_block$chr, ld_block$start, ld_block$stop)

    extracted_trans_hit <- data.frame(chr = as.character(clumped_snp['CHR']),
                                    bp = clumped_snp['BP'],
                                    log_p = -log10(clumped_snp['P']),
                                    ancestry = study$ancestry,
                                    ld_block = ld_block_string,
                                    file = extracted_file,
                                    cis_trans = 'trans',
                                    reference_build=study$reference_build
    )
    return(extracted_trans_hit)
  }) |> dplyr::bind_rows()

  return(list(snp_data=extracted_trans_snps, clumped_snps=NULL))
}

format_gwas <- function(gwas) {
  gwas <- dplyr::rename(gwas, RSID='SNP', CHR='Chr', EA='A1', OA='A2', EAF='Freq', BETA='b', P='p') |>
    dplyr::select(-Probe, -Probe_Chr, -Probe_bp, -Gene, -Orientation)

  return(gwas)
}

main()
