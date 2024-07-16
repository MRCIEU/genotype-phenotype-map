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
study <- vroom::vroom(paste0(pipeline_metadata_dir, '/studies_to_process.tsv'), show_col_types = F) |>
  dplyr::filter(extracted_location == args$extracted_study_location)

if (nrow(study) != 1) quit('Error: cant find study to process', status = 1)

p_value <- ifelse(is.na(study$p_value_threshold), DEFAULT_P_VALUE_THRESHOLD, study$p_value_threshold)

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
  print(result)
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
