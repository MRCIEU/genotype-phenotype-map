source('constants.R')
library(argparser, quietly = TRUE)

parser <- arg_parser("Sorting into LD blocks and colocalising results")
parser <- add_argument(parser,
                       "--chr",
                       help = "Chromosome to run script on",
                       type = "numeric"
)

load("data/ldetect.rdata")
current_state <- vroom::vroom(paste0(data_dir, "/pipeline_metadata/current_state.tsv"), show_col_types = F) |>
  dplyr::filter(extracted == F)

updated_ld_blocks <- apply(current_state, 1, function(study) {
  study_dir <- study[['extracted_location']]
  extracted_snps <- vroom::vroom(paste0(study_dir, "/extracted_snps.tsv"), show_col_types = F) |>
    dplyr::filter(CHR == args$chr)

  updated_ld_blocks <- apply(extracted_snps, 1, function(extracted) {
    bp <- as.numeric(extracted[['BP']])
    extracted_chr <- as.numeric(extracted[['CHR']])
    ancestry <- extracted[['ANCESTRY']]
    ld_block <- dplyr::filter(ldetect, chr == paste0("chr", extracted_chr) & start < bp & stop > bp & pop == ancestry)
    print(ld_block)
    print(nrow(ld_block))
    if (nrow(ld_block) > 1) stop(paste("Error: More than 1 LD Block associated with", extracted_chr, bp))

    ld_block_dir <- paste0(ld_block_dir, ancestry, "/", extracted_chr, "/", ld_block$start, "_", ld_block$stop)
    print(ld_block_dir)
    if(!dir.exists(ld_block_dir)) dir.create(ld_block_dir, recursive = T)

    #TODO: change to tsv.gz
    study_file <- paste0(study_dir, ancestry, "_", extracted_chr, "_", bp, ".z")

    file.symlink(study_file, ld_block_dir)
    return(ld_block)
  }) |> dplyr::bind_rows()
  return(updated_ld_blocks)
})

updated_ld_blocks <- dplyr::bind_rows(updated_ld_blocks) |> dplyr::distinct()


#TODO: actually run hyprcoloc below
hyprcoloc::hyprcoloc()