source("../../pipeline_steps/constants.R")
ld_regions <- vroom::vroom("../../pipeline_steps/data/ld_regions.tsv", show_col_types = F)
ld_info <- construct_ld_block(ld_regions$ancestry, ld_regions$chr, ld_regions$start, ld_regions$stop)

NUM_PARALLEL_JOBS = 150
results <- parallel::mclapply(X=ld_info$ld_reference_panel_prefix, mc.cores=NUM_PARALLEL_JOBS, FUN=function(ld_region_string) {
  print(ld_region_string)

  bim <- vroom::vroom(paste0(ld_region_string, ".bim"), delim='\t', col_names = F, show_col_types = F) |>
    dplyr::select(X2, X4)
  gwas <- data.table::fread(paste0(ld_region_string, ".frq")) |>
    dplyr::rename(RSID='SNP', EA='A1', OA='A2', EAF='MAF')

  gwas$BP <- bim$X4[match(gwas$RSID, bim$X2)]
  gwas <- dplyr::select(gwas, RSID, CHR, BP, EA, OA, EAF, -NCHROBS)
  gwas <- standardise_alleles(gwas)

  ld_matrix <- vroom::vroom(paste0(ld_region_string, '.ld'), col_names = F, show_col_types=F)

  if (nrow(ld_matrix) == ncol(ld_matrix)) {
    print('already processed, skipping')
    return()
  }

  keep <- !is.na(ld_matrix$X1)
  print(paste('deleting', sum(!keep)))
  ld_matrix <- ld_matrix[keep, c(keep, F)]

  if (sum(!keep) > 0) {
    gwas <- gwas[keep, ]
  }

  vroom::vroom_write(ld_matrix, paste0(ld_region_string, '.ld'), col_names = F, delim=' ')
  vroom::vroom_write(gwas, paste0(ld_region_string, ".tsv"))
})
