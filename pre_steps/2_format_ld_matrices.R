source("constants.R")
ld_regions <- vroom::vroom("data/ld_regions.tsv", show_col_types = F)

for (region in seq_len(nrow(ld_regions))) {
  ld_region <- ld_regions[region,]
  ld_region_string <- paste0(ld_block_matrices_dir, ld_region$ancestry, "/", ld_region$chr, "_", ld_region$start, "_", ld_region$stop)

  bim <- vroom::vroom(paste0(ld_region_string, ".bim"), delim='\t', col_names = F, show_col_types = F) |>
    dplyr::select(X2, X4)
  gwas <- data.table::fread(paste0(ld_region_string, ".frq")) |>
    dplyr::rename(RSID='SNP', EA='A1', OA='A2', EAF='MAF')
  gwas$BP <- bim$X4[match(gwas$RSID, bim$X2)]
  gwas <- dplyr::select(gwas, RSID, CHR, BP, EA, OA, EAF, -NCHROBS)
  gwas <- standardise_alleles(gwas)

  vroom::vroom_write(gwas, paste0(ld_region_string, ".tsv"))
}
