source('constants.R')
# 1. Clump GWAS at requested --value threshold to find all regions to extract
# 2. Deduplicate list to ensure we are only extracting 1 SNP per region (and using the top hit)
#  with these 2 steps, I can create extracted_snps.tsv, but with cis_trans filled out appropriately
# 3. Using extracted_snps, extract whole regions and store

parser <- argparser::arg_parser('Extract genomic regions from a BESD file')
parser <- argparser::add_argument(parser, '--extracted_study_location', help = 'Study location', type = 'character')
parser <- argparser::add_argument(parser, '--extracted_output_file', help = 'Extracted Output file', type = 'character')
args <- argparser::parse_args(parser)

ld_blocks <- vroom::vroom('data/ld_blocks.tsv', show_col_types = F)

main <- function() {
  study <- vroom::vroom(glue::glue('{pipeline_metadata_dir}/studies_to_process.tsv'), show_col_types = F) |>
    dplyr::filter(extracted_location == args$extracted_study_location)
  if (nrow(study) != 1) stop('Error: cant find study to process')

  p_value_threshold <- ifelse(is.na(study$p_value_threshold), lowest_p_value_threshold, study$p_value_threshold)

  metadata <- jsonlite::fromJSON(glue::glue('{study$study_location}/{study$study_name}.json'))
  vcf_file <- glue::glue('{study$study_location}/{study$study_name}.vcf.gz')
  clumped_hits_file <- glue::glue('{args$extracted_study_location}/clumped_snps.tsv')

  if (study$reference_build == reference_builds$GRCh37) {
    vcf_file <- convert_reference_build(study, vcf_file)
  }

  clumped_snps <- find_clumped_hits(study, vcf_file, p_value_threshold)
  vroom::vroom_write(clumped_snps, clumped_hits_file)

  extracted_snps <- extract_clumped_regions(study, vcf_file, clumped_snps)
  create_svgs_from_gwas(study)

  vroom::vroom_write(extracted_snps, args$extracted_output_file)
}

convert_reference_build <- function(study,
                                    vcf_file,
                                    input_reference_build=reference_builds$GRCh37,
                                    output_reference_build=reference_builds$GRCh38) {

  if (input_reference_build == output_reference_build) return(vcf_file)

  dir.create(glue::glue('{study$extracted_location}/vcf'), showWarnings = F, recursive = T)
  liftover_conversion <- available_liftover_conversions[[glue::glue('{input_reference_build}{output_reference_build}')]]
  if (is.null(liftover_conversion)) {
    stop(paste(c("Error: liftOver combination of", input_build, output_build, "not recocognised.",
                 "Reference builds must be one of:", reference_builds), collapse = " "))
  }
  output_file <- glue::glue('{study$extracted_location}vcf/hg38.vcf.gz')
  rejected_file <- glue::glue('{study$extracted_location}vcf/hg38_rejected.vcf')
  fasta_file <- glue::glue('{liftover_dir}/hg38.fa')

  if (file.exists(output_file)) {
    message(glue::glue('{vcf_file} already converted to hg38.'))
    return(output_file)
  }

  bcf_liftover_command <- glue::glue(
    '/home/bcftools/bcftools annotate --rename-chrs {liftover_dir}/num_to_chr.txt {vcf_file} | ',
      '/home/bcftools/bcftools +liftover --no-version -Ou -- ',
      '-s {liftover_dir}/hg19.fa ',
      '-f {liftover_dir}/hg38.fa ',
      '-c {liftover_conversion} ',
      '--reject {rejected_file} | ',
        '/home/bcftools/bcftools annotate --rename-chrs {liftover_dir}/chr_to_num.txt | ',
          '/home/bcftools/bcftools sort -Oz -o {output_file} -W=tbi'
  )
  system(bcf_liftover_command, wait = T, ignore.stdout = T)

  return(output_file)
}

find_clumped_hits <- function(study, vcf_file, p_value_threshold) {
  gwasvcf::set_bcftools('/home/bcftools/bcftools')

  significant_hits <- gwasvcf::query_gwas(vcf_file, pval=p_value_threshold) |>
    gwasvcf::vcf_to_tibble() |>
    dplyr::mutate(pval=10^{-LP}, rsid=ID) |>
    dplyr::select(rsid, pval)

  significant_hits_file <- withr::local_tempfile()
  clumped_hits_file <- withr::local_tempfile()
  vroom::vroom_write(significant_hits, significant_hits_file, delim=' ')

  bfile <- glue::glue('{ld_reference_panel_dir}{study$ancestry}/full_rsid') 
  plink_command <- glue::glue('plink2 --bfile {bfile} ',
    '--clump {significant_hits_file} ',
    '--clump-p1 {p_value_threshold} ',
    '--clump-r2 0.1 ',
    '--clump-kb 1000 ',
    '--clump-snp-field rsid ',
    '--clump-field pval ',
    '--out {clumped_hits_file}'
  )
  system(plink_command, wait = T, ignore.stdout = T)

  clumped_snps <- data.table::fread(glue::glue('{clumped_hits_file}.clumps')) |>
    dplyr::rename(RSID='ID', CHR='#CHROM', BP='POS') |>
    dplyr::select(RSID, CHR, BP, P) |>
    dplyr::mutate(CHR=as.numeric(CHR), BP=as.numeric(BP), P=as.numeric(P)) |>
    dplyr::arrange(P)

  return(clumped_snps)
}

extract_clumped_regions <- function(study, vcf_file, clumped_snps) {
  extracted_snps <- data.frame(
    chr=character(),
    bp=numeric(),
    log_p=numeric(),
    ld_block=character(),
    file=character(),
    cis_trans=character()
  )

  if (nrow(clumped_snps) == 0) {
    message(glue::glue('{study["study_location"]}: No clumped results'))
    return(extracted_snps)
  }

  dir.create(glue::glue('{study$extracted_location}/svgs'), showWarnings = F, recursive = T)
  dir.create(glue::glue('{study$extracted_location}/extracted'), showWarnings = F, recursive = T)
  dir.create(glue::glue('{study$extracted_location}/standardised'), showWarnings = F, recursive = T)
  dir.create(glue::glue('{study$extracted_location}/imputed'), showWarnings = F, recursive = T)
  dir.create(glue::glue('{study$extracted_location}/finemapped'), showWarnings = F, recursive = T)

  ld_blocks$bcf_region <- glue::glue('{ld_blocks$chr}:{ld_blocks$start}-{ld_blocks$stop}') 
  ld_blocks$string_region <- ld_block_string(ld_blocks$ancestry, ld_blocks$chr, ld_blocks$start, ld_blocks$stop)

  clumped_snps <- apply(clumped_snps, 1, function(clump) {
    clump_chr <- as.numeric(clump['CHR']) 
    clump_bp <- as.numeric(clump['BP']) 

    region <- dplyr::filter(ld_blocks, chr == clump_chr & start <= clump_bp & stop > clump_bp & ancestry == study$ancestry)
    if (nrow(region) == 0) {
      # message(glue::glue('no region for {clump_chr}:{clump_bp}:{study$ancestry}'))
      return()
    }

    return(data.frame(CHR = clump_chr,
      BP = clump_bp,
      P = as.numeric(clump['P']),
      bcf_region = region$bcf_region,
      region_start = region$start,
      region_stop = region$stop,
      string_region = region$string_region)
    )
  }) |> dplyr::bind_rows()

  #removing duplicate entries per region, so we only grab the region once.
  clumped_snps <- clumped_snps[!duplicated(clumped_snps$bcf_region), ]
  message(glue::glue('num to extract: {nrow(clumped_snps)}'))

  if (nrow(clumped_snps) == 0) {
    return(data.frame())
  }

  regions_file <- withr::local_tempfile()
  vroom::vroom_write(dplyr::select(clumped_snps, CHR, region_start, region_stop), regions_file, delim = '\t', col_names = F)
  
  bcf_query <- glue::glue('/home/bcftools/bcftools query ',
    '--regions-file {regions_file} ',
    '--format "[%ID]\t[%CHROM]\t[%POS]\t[%ALT]\t[%REF]\t[%AF]\t[%ES]\t[%SE]\t[%LP]" ',
    '{vcf_file}'
  )
  extracted_regions <- system(bcf_query, wait = T, intern = T)
  extracted_regions <- data.table::fread(text = extracted_regions)
  colnames(extracted_regions) <- c('RSID', 'CHR', 'BP', 'EA', 'OA', 'EAF', 'BETA', 'SE', 'LP')

  extracted_snp_info <- apply(clumped_snps, 1, function(clump) {
    extracted_region <- dplyr::filter(extracted_regions, CHR == as.numeric(clump['CHR']) & BP >= as.numeric(clump['region_start']) & BP <= as.numeric(clump['region_stop']))

    extracted_file <- glue::glue('{study$extracted_location}extracted/{study$ancestry}_{clump["CHR"]}_{clump["BP"]}.tsv.gz')
    extracted_file <- gsub(' ', '', extracted_file)
    vroom::vroom_write(extracted_region, extracted_file)

    extraction_info <- data.frame(chr = as.character(clump[['CHR']]),
                                  bp = as.numeric(clump[['BP']]),
                                  log_p = -log10(as.numeric(clump[['P']])),
                                  ld_block = clump[['string_region']],
                                  file = extracted_file,
                                  cis_trans = NA
    )
    return(extraction_info)
  }) |> dplyr::bind_rows()

  message(glue::glue('{study["extracted_location"]}: Extracted {nrow(extracted_snp_info)} regions'))
  return(extracted_snp_info)
}

main()
