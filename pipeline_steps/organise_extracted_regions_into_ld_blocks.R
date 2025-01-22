source('constants.R')

parser <- argparser::arg_parser('Organise Extracted Regions into LD regions')
parser <- argparser::add_argument(parser, '--output_file', help = 'Output file', type = 'character')
args <- argparser::parse_args(parser)

ld_blocks <- vroom::vroom('data/ld_blocks.tsv', show_col_types = F)
studies_to_process <- vroom::vroom(glue::glue('{data_dir}/pipeline_metadata/studies_to_process.tsv'), show_col_types = F)

all_extracted_snp_files <- glue::glue('{studies_to_process$extracted_location}extracted_snps.tsv')
#filtering out results without any significant SNPs, so we don't hit ulimits on the box
all_extracted_snp_files  <- Filter(function(file) file.info(file)$size > 70, all_extracted_snp_files)
all_extracted_snps <- vroom::vroom(all_extracted_snp_files, show_col_types = F)

all_extracted_snps$study_name <- stringr::str_extract(all_extracted_snps$file, '(?<=study/)[\\w-_:]+')
extracted_snps_by_region <- split(all_extracted_snps, all_extracted_snps$ld_block)

results <- lapply(extracted_snps_by_region, function(extracted_snps) {
  ld_block <- unique(extracted_snps$ld_block)
  ld_info <- ld_block_dirs(ld_block)
  merged_data <- merge(extracted_snps, studies_to_process, by='study_name')
  if (nrow(extracted_snps) == 0) return()

  extracted_studies <- data.frame(study = merged_data$study_name,
                                  file = merged_data$file,
                                  ancestry = merged_data$ancestry,
                                  chr = merged_data$chr,
                                  bp = merged_data$bp,
                                  p_value_threshold = merged_data$p_value_threshold,
                                  category = merged_data$category,
                                  sample_size = merged_data$sample_size,
                                  cis_trans = merged_data$cis_trans,
                                  reference_build = reference_builds$GRCh38,
                                  ld_block = ld_info$block,
                                  variant_type = merged_data$variant_type
  )

  if (!dir.exists(ld_info$ld_block_data)) dir.create(ld_info$ld_block_data, recursive = T)
  extracted_studies_file <- glue::glue('{ld_info$ld_block_data}/extracted_studies.tsv')

  if (file.exists(extracted_studies_file)) {
    existing_extracted_studies <- vroom::vroom(extracted_studies_file, show_col_types = F)
    extracted_studies <- rbind(existing_extracted_studies, extracted_studies)
    extracted_studies <- extracted_studies[!duplicated(extracted_studies), ]
  }

  vroom::vroom_write(extracted_studies, extracted_studies_file)
})

ld_info <- construct_ld_block(ld_blocks$ancestry, ld_blocks$chr, ld_blocks$start, ld_blocks$stop)
ld_blocks$ld_block <- ld_info$block 
ld_blocks$data_dir <- ld_info$ld_block_data
all_updated_ld_blocks <- dplyr::filter(ld_blocks, dir.exists(data_dir)) |> dplyr::arrange(chr)

vroom::vroom_write(all_updated_ld_blocks, args$output_file)
