source("constants.R")
ld_regions <- vroom::vroom("data/ld_regions.tsv", show_col_types = F)

for (region in seq_len(nrow(ld_regions))) {
  ld_region <- ld_regions[region,]
  ld_region_string <- paste0(ld_block_matrices_dir, ld_region$pop, "/", ld_region$chr, "_", ld_region$start, "_", ld_region$stop)
  print(ld_region_string)

  snp_list <- vroom::vroom(paste0(ld_region_string, ".snplist"), delim = " ", col_names = F, show_col_types = F)
  ld_matrix <- vroom::vroom(paste0(ld_region_string, ".ld"), delim = " ", col_names = F, show_col_types = F)
  ld_matrix <- dplyr::bind_cols(snp_list, ld_matrix)

  vroom::vroom_write(ld_matrix, paste0(ld_region_string, ".tsv.gz"), col_names = F)
}
