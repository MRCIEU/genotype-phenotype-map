source("constants.R")

ld_regions <- vroom::vroom("data/ld_regions.tsv", show_col_types = F)

for (region in seq_len(nrow(ld_regions))) {
  ld_region <- ld_regions[region,]
  ld_region_string <- paste0(ld_block_matrices_dir, ld_region$pop, "/", ld_region$chr, "_", ld_region$start, "_", ld_region$stop, ".ld")

  ld_matrix <- vroom::vroom(ld_region_string, col_names = F, show_col_types=F)

  if (nrow(ld_matrix) == ncol(ld_matrix)) {
    next
  }

  keep <- !is.na(ld_matrix$X1)

  print(paste(ld_region_string))
  print(paste('deleting', sum(!keep)))

  ld_matrix <- ld_matrix[keep, c(keep, F)]
  vroom::vroom_write(ld_matrix, ld_region_string, col_names = F, delim=' ')

  if (sum(!keep) > 0) {
    ld_data_string <- paste0(ld_block_matrices_dir, ld_region$pop, "/", ld_region$chr, "_", ld_region$start, "_", ld_region$stop, ".tsv")
    ld_data <- vroom::vroom(ld_data_string, show_col_types=F)
    ld_data <- ld_data[keep, ]
    vroom::vroom_write(ld_data, ld_data_string)
  }
}

