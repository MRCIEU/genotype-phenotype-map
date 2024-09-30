source('constants.R')
# 1. Clump GWAS at requested --value threshold to find all regions to extract
# 2. Deduplicate list to ensure we are only extracting 1 SNP per region (and using the top hit)
#  with these 2 steps, I can create extracted_snps.tsv, but with cis_trans filled out appropriately
# 3. Using extracted_snps, extract whole regions and store

parser <- argparser::arg_parser('Extract genomic regions from a BESD file')
parser <- argparser::add_argument(parser, '--extracted_study_location', help = 'Study location', type = 'character')
parser <- argparser::add_argument(parser, '--extracted_output_file', help = 'Extracted Output file', type = 'character')
args <- argparser::parse_args(parser)

ld_regions <- vroom::vroom('data/ld_regions.tsv', show_col_types = F)

main <- function() {
  study <- vroom::vroom(paste0(pipeline_metadata_dir, '/studies_to_process.tsv'), show_col_types = F) |>
    dplyr::filter(extracted_location == args$extracted_study_location)
  if (nrow(study) != 1) stop('Error: cant find study to process')

  p_value_threshold <- ifelse(is.na(study$p_value_threshold), genome_wide_p_value_threshold, study$p_value_threshold)
  metadata <- jsonlite::fromJSON(glue::glue('{study$study_location}/{study$study_name}.json'))
  vcf_file <- glue::glue('{study$study_location}/{study$study_name}.vcf.gz')

  if (study$reference_build == reference_builds$GRCh37) {
    vcf_file <- convert_reference_build_using_picard(study, vcf_file)
  }

  clumped_snps <- find_clumped_hits(study, vcf_file, p_value_threshold)
  vroom::vroom_write(clumped_snps, glue::glue('{args$extracted_study_location}/clumped_snps.tsv'))
  extracted_snps <- extract_clumped_regions(study, vcf_file, clumped_snps)

  vroom::vroom_write(extracted_snps, args$extracted_output_file)
}

convert_reference_build_using_picard <- function(study,
                                                 vcf_file,
                                                 input_reference_build=reference_builds$GRCh37,
                                                 output_reference_build=reference_builds$GRCh38) {

  if (input_reference_build == output_reference_build) return(vcf_file)

  dir.create(glue::glue('{study$extracted_location}/vcf'), showWarnings = F, recursive = T)
  liftover_conversion <- available_liftover_conversions[[paste0(input_reference_build, output_reference_build)]]
  if (is.null(liftover_conversion)) {
    stop(paste(c("Error: liftOver combination of", input_build, output_build, "not recocognised.",
                 "Reference builds must be one of:", reference_builds), collapse = " "))
  }
  output_file <- glue::glue('{study$extracted_location}vcf/hg38.vcf.gz')
  rejected_file <- glue::glue('{study$extracted_location}vcf/hg38_rejected.vcf.gz')
  fasta_file <- glue::glue('{liftover_dir}/hg38.fa')

  # https://broadinstitute.github.io/picard/command-line-overview.html
  # picard_command <- glue::glue(
    # 'java -jar /usr/bin/picard.jar LiftoverVcf ',
    # 'I={vcf_file} O={output_file} ',
    # 'CHAIN={liftover_conversion} REJECT={rejected_file} R={fasta_file}') #dont know what R does

  bcf_liftover_command <- glue::glue(
    '/home/bcftools/bcftools annotate --rename-chrs {liftover_dir}/chr_conversion.txt {vcf_file} | ',
      '/home/bcftools/bcftools +liftover --no-version -Ou -- ',
      '-s {liftover_dir}/hg19.fa ',
      '-f {liftover_dir}/hg38.fa ',
      '-c {liftover_conversion} ',
      '--reject {rejected_file} | ',
      # '--reject-type z'
      # '--reject-type z | ',
        '/home/bcftools/bcftools sort -Oz -o {output_file} -W=tbi'
  )
  message(bcf_liftover_command)
  system(bcf_liftover_command, wait = T, ignore.stdout = T)

  return(output_file)
}

find_clumped_hits <- function(study, vcf_file, p_value_threshold) {
  gwasvcf::set_bcftools('/home/bcftools/bcftools')

  significant_hits <- gwasvcf::query_gwas(vcf_file, pval=p_value_threshold) |>
		gwasvcf::vcf_to_tibble() |>
		dplyr::mutate(pval=10^{-LP}) |>
		dplyr::select(rsid, pval)

  significant_hits_file <- tempfile()
  clumped_hits_file <- tempfile()
  vroom::vroom_write(significant_hits, significant_hits_file, delim=' ')

  bfile <- glue::glue('{ld_reference_panel_dir}{study$ancestry}/full_rsid') 
  plink_command <- glue::glue('plink2 --bfile {bfile} ',
    '--clump {significant_hits_file} ',
    '--clump-p1 {p_value_threshold} ',
    '--clump-r2 0.1 ',
    '--clump-kb 1000 ',
    '--clump-snp-field rsid ',
    '--clump-field pval ',
    '--out {clumped_hits_file}')
  system(plink_command, wait = T, ignore.stdout = T)

  clumped_snps <- data.table::fread(glue::glue('{clumped_hits_file}.clumps'))

  clumped_snps <- data.table::fread(glue::glue('{clumped_hits_file}.clumps')) |>
    dplyr::rename(RSID='ID', CHR='#CHROM', BP='POS') |>
    dplyr::select(RSID, CHR, BP, P) |>
    dplyr::arrange(P)

  return(clumped_snps)
}

extract_clumped_regions <- function(study, vcf_file, clumped_snps) {
  extracted_snps <- data.frame(
    chr=character(),
    bp=numeric(),
    log_p=numeric(),
    ld_region=character(),
    file=character(),
    cis_trans=character()
  )

  if (nrow(clumped_snps) == 0) {
    message(glue::glue('{study["study_location"]}: No clumped results'))
    return(extracted_snps)
  }

  dir.create(glue::glue('{study$extracted_location}/original'), showWarnings = F, recursive = T)
  dir.create(glue::glue('{study$extracted_location}/standardised'), showWarnings = F, recursive = T)
  dir.create(glue::glue('{study$extracted_location}/imputed'), showWarnings = F, recursive = T)
  dir.create(glue::glue('{study$extracted_location}/finemapped'), showWarnings = F, recursive = T)

  ld_regions$bcf_region <- glue::glue('{ld_regions$chr}:{ld_regions$start}-{ld_regions$stop}') 
  ld_regions$string_region <- ld_block_string(ld_regions$ancestry, ld_regions$chr, ld_regions$start, ld_regions$stop)

  clumped_snps <- apply(clumped_snps, 1, function(clump) {
    region <- dplyr::filter(ld_regions, chr == clump['CHR'] & start < as.numeric(clump['BP']) & stop > as.numeric(clump['BP']) & ancestry == study$ancestry)
    if (nrow(region) == 0) return()

    return(data.frame(CHR = clump['CHR'], BP = as.numeric(clump['BP']), P=clump['P'], bcf_region = region$bcf_region, string_region = region$string_region))
  }) |> dplyr::bind_rows()
  
  #removing duplicate entries per region, so we only grab the region once.
  clumped_snps <- clumped_snps[!duplicated(clumped_snps$bcf_region), ]
  print(glue::glue('num to extract: {nrow(clumped_snps)}'))

  extracted_snps <- apply(clumped_snps, 1, function(clump) {
    bcf_query <- glue::glue('/home/bcftools/bcftools query ',
      '--regions {clump[["bcf_region"]]} ',
      '--format "[%ID]\t[%CHROM]\t[%POS]\t[%REF]\t[%ALT]\t[%AF]\t[%ES]\t[%SE]\t[%LP]" ',
      '{vcf_file}'
    )

    extracted_region <- system(bcf_query, wait = T, intern = T)
    extracted_region <- data.table::fread(text = extracted_region)
    colnames(extracted_region) <- c('RSID', 'CHR', 'BP', 'EA', 'OA', 'EAF', 'BETA', 'SE', 'LP')

    extracted_file <- glue::glue('{study$extracted_location}original/{study$ancestry}_{clump["CHR"]}_{clump["BP"]}.tsv.gz')
    extracted_file <- gsub(' ', '', extracted_file)
    vroom::vroom_write(extracted_region, extracted_file)

    extraction_info <- data.frame(chr = as.character(clump[['CHR']]),
                                  bp = as.numeric(clump[['BP']]),
                                  log_p = -log10(as.numeric(clump[['P']])),
                                  ld_region = clump[['string_region']],
                                  file = extracted_file,
                                  cis_trans = NA
    )
    return(extraction_info)
  }) |> dplyr::bind_rows()

  message(glue::glue('{study["extracted_location"]}: Extracted {nrow(extracted_snps)} regions'))
  return(extracted_snps)
}

main()
