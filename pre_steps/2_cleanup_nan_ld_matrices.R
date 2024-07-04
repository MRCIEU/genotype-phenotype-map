source("constants.R")

ld_regions <- vroom::vroom("data/ld_regions.tsv", show_col_types = F)

for (region in seq_len(nrow(ld_regions))) {
  ld_region <- ld_regions[region,]
  ld_region_string <- paste0(ld_block_matrices_dir, ld_region$pop, "/", ld_region$chr, "_", ld_region$start, "_", ld_region$stop, ".ld")

  ld_matrix <- vroom::vroom(ld_region_string, col_names = F)
  keep <- !is.na(ld_matrix$X1)
  ld_matrix <- ld_matrix[keep, keep]

  vroom::vroom_write(gwas, ld_region_string, col_names = F)
}
