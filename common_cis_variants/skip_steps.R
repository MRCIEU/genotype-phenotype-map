source("constants.R")
ld_blocks <- vroom::vroom('data/ld_regions.tsv', show_col_types=F)

for (i in 1:nrow(ld_blocks)) {
  ld_block <- ld_blocks[i, ]
  ld_block_dir <- paste0(ld_block_data_dir, ld_block$pop, "/", ld_block$chr, "/", ld_block$start, "_", ld_block$stop)
  vroom::vroom_write(data.frame(), paste0(ld_block_dir, '/imputation_complete'))
  #vroom::vroom_write(data.frame(), paste0(ld_block_dir, '/finemapping_complete'))
}

