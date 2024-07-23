source('constants.R')
#TODO: 3 step process per 'study' (<data_source>-<specifier>-<gene>)
# 1. Find cis snps (via cis-wind)
# 2. Find all other significant hits, which will be considered trans
#  with these 2 steps, I can create extracted_snps.tsv, but with cis_trans filled out appropriately
# 3. Using extracted_snps, extract whole regions

parser <- argparser::arg_parser('Extract genomic regions from a BESD file')
parser <- argparser::add_argument(parser, '--extracted_study_location', help = 'Study location', type = 'character')
parser <- argparser::add_argument(parser, '--extracted_output_file', help = 'Extracted Output file', type = 'character')
args <- argparser::parse_args(parser)

ld_regions <- vroom::vroom('data/ld_regions.tsv', show_col_types = F)

main <- function(args) {
  study <- vroom::vroom(paste0(pipeline_metadata_dir, '/studies_to_process.tsv'), show_col_types = F) |>
    dplyr::filter(extracted_location == args$extracted_study_location)
  if (nrow(study) != 1) stop('Error: cant find study to process')

  p_value_threshold <- ifelse(is.na(study$p_value_threshold), DEFAULT_P_VALUE_THRESHOLD, study$p_value_threshold)
  metadata <- jsonlite::fromJSON(paste0(study$study_location, '.json'))

  extracted_cis_snps <- extracted_trans_snps <- data.frame(chr=character(),
                                                           bp=numeric(),
                                                           log_p=numeric(),
                                                           ancestry=character(),
                                                           ld_region=character(),
                                                           file=character(),
                                                           cis_trans=character()
  )
  if (metadata$cis_trans == cis_trans$cis_only || metadata$cis_trans == cis_trans$cis_trans) {
    extracted_cis_snps <- extract_cis_region(study, p_value_threshold)
  }
  if (metadata$cis_trans == cis_trans$trans_only || metadata$cis_trans == cis_trans$cis_trans) {
    extracted_trans_snps <- extract_trans_regions(study, p_value_threshold)
  }

  extracted_snps <- dplyr::bind_rows(extracted_cis_snps, extracted_trans_snps)
  vroom::vroom_write(extracted_snps, args$extracted_output_file)
}


extract_cis_region <- function(study, p_value_threshold) {
  tmp_cis_snp <- paste0('/tmp/', study$study_name, '_top_snp')
  extract_top_snp <- paste('smr --beqtl-summary', study$study_location,
                           '--query', p_value_threshold,
                           '--probe ', study$gene,
                           '--cis-wind 1',
                           '--out ', tmp_cis_snp
  )

  system(extract_top_snp, wait=T)
  if (!file.exists(paste0(tmp_cis_snp, '.txt'))) return()
  top_cis_snp <- vroom::vroom(paste0(tmp_cis_snp, '.txt'), show_col_types = F)
  if (nrow(top_cis_snp) < 1) return()
  top_cis_snp <- top_cis_snp[top_cis_snp$p == min(top_cis_snp$p), ][1, ]

  ld_block <- dplyr::filter(ld_regions, chr == top_cis_snp$Chr & start < top_cis_snp$BP & stop > top_cis_snp$BP & ancestry == study$ancestry)
  ld_block_string <- ld_block_string(ld_block$ancestry, ld_block$chr, ld_block$start, ld_block$stop)
  ld_info <- ld_block_dirs(ld_block_string)

  tmp_cis_region <- paste0('/tmp/', study$study_name, '_cis_region')
  extract_region <- paste('smr --beqtl-summary', study$study_location,
                              '--query 1',
                              '--snp', top_cis_snp$SNP,
                              '--snp-wind 4000', # smr doesn't accept BP ranges, and errors if specific RSID isn't present
                              '--probe ', study$gene,
                              '--out ', tmp_cis_region
  )
  system(extract_region, wait=T)
  cis_region <- vroom::vroom(paste0(tmp_cis_region, '.txt'), show_col_types = F)
  cis_region <- format_gwas(cis_region)

  extracted_file <-paste0(study$extracted_location, '/original/', study$ancestry, '_', top_cis_snp$Chr, '_', top_cis_snp$BP, '.tsv.gz')
  vroom::vroom_write(cis_region, extracted_file)

  top_hits_per_probe <- data.frame(chr = as.character(top_cis_snp$Chr),
                                   bp = top_cis_snp$BP,
                                   log_p = -log10(top_cis_snp$p),
                                   ancestry = study$ancestry,
                                   ld_region = ld_block_string,
                                   file = extracted_file,
                                   cis_trans = 'cis'
  )
  return(top_hits_per_probe)
}

#TODO: not complete
extract_trans_regions <- function() {
  tmp_smr_result <- paste0('/tmp/', study$gene)
  extract_to_cis_hits <- paste('smr --beqtl-summary', study$study_location,
                               '--query', p_value,
                               '--probe ', study$gene,
                               '--cis-wind 1000',
                               '--out ', tmp_smr_result
  )
  system(extract_to_cis_hits, wait=T)
  probe_top_hits <- vroom::vroom(paste0(tmp_smr_result, '.txt'), show_col_types = F)

  if (nrow(probe_top_hits) < 1) return()
  top_hits_per_probe <- data.frame(file_prefix = study$study_location, SNP = probe_top_hits$SNP, probe = study$gene, cis_trans = 'cis')

  vroom::vroom_write(top_hits_per_probe, 'all_cis_top_hits.tsv')
  print('starting trans')

  all_results <- apply(top_hits_per_probe, 1, function(result) {
    tmp_smr_result <- paste0('/tmp/', result['probe'], '_', result['SNP'])
    extract_trans_snps <- paste('smr --beqtl-summary', result['file_prefix'],
                                '--query', p_value,
                                '--snp ', result['SNP'],
                                '--out ', tmp_smr_result
    )
    print('extract_trans_snps')
    print(extract_trans_snps)
    system(extract_trans_snps, wait=T)
    probe_top_hits <- vroom::vroom(paste0(tmp_smr_result, '.txt'), show_col_types = F)
    print(probe_top_hits)
    if (nrow(probe_top_hits) <= 1) return()
    return(probe_top_hits)
  }) |> dplyr::bind_rows()



  vroom::vroom_write(all_results, paste0(basename(study$study_location), '_top_hits.tsv'))
  all_results <- dplyr::left_join(all_results, all_cis_top_hits, by=dplyr::join_by(SNP==SNP, Probe==probe)) |>
    dplyr::select(-file_prefix)
  all_results$cis_trans[is.na(all_results$cis_trans)] <- 'trans'

  apply(all_results, 1, function(result) {
    bp <- as.numeric(result['BP'])
    extracted_chr <- result['Chr']
    ld_block <- dplyr::filter(ld_regions, chr == extracted_chr & start < bp & stop > bp & ancestry == study['ancestry'])
    if (nrow(ld_block) > 1) stop(paste('Error: More than 1 LD Block associated with', extracted_chr, bp))

    tmp_smr_result <- paste0('/tmp/', result['probe'], '_', result['SNP'], '_region')
    extract_region <- paste('smr --beqtl-summary',result['file_prefix'],
                            '--query 1',
                            '--snp-chr',extracted_chr,
                            '--from-snp',ld_block$start,
                            '--to-snp',ld_block$stop,
                            '--out', tmp_smr_result
    )
    system(extract_region, wait = T)
  })


  vroom::vroom_write(all_results, paste0(basename(file_prefix), 'merged_top_hits.tsv'))
  print(all_results)

  vroom::vroom_write(extracted_snps, args$extracted_output_file)

}

format_gwas <- function(gwas) {
  gwas <- dplyr::rename(gwas, RSID='SNP', CHR='Chr', EA='A1', OA='A2', EAF='Freq', BETA='b', P='p') |>
    dplyr::select(-Probe, -Probe_Chr, -Probe_bp, -Gene, -Orientation)
  return(gwas)
}


main(args)
