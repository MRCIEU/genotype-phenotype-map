source('../../pipeline_steps/constants.R')
ld_blocks <- vroom::vroom('../../pipeline_steps/data/ld_blocks.tsv', show_col_types = F) |>
  dplyr::filter(ancestry == 'EUR')
ld_info <- construct_ld_block(ld_blocks$ancestry, ld_blocks$chr, ld_blocks$start, ld_blocks$stop)

merge_all_bfiles <- function() {
  mergefile <- tempfile()
  write.table(ld_info$ld_reference_panel_prefix, file=mergefile, row=F, col=F, qu=F)
  output <- glue::glue("{ld_reference_panel_dir}EUR/full")
  glue::glue("plink1.9 --merge-list {mergefile} --make-bed --out {output} --keep-allele-order") |> system()
}


create_tsv_from_afreq <- function() {
  lapply(ld_info$ld_reference_panel_prefix, function(prefix) {
    freq <- vroom::vroom(glue::glue('{prefix}.afreq'), show_col_types = F) |>
      dplyr::rename(CHR='#CHROM', SNP='ID', EA='ALT', OA='REF', EAF='ALT_FREQS') |>
      dplyr::mutate(BP=as.numeric(sub('.*:(\\d+)_.*', '\\1', SNP))) |>
      dplyr::select(-`PROVISIONAL_REF?`, -OBS_CT)
    print(head(freq))

    vroom::vroom_write(freq, glue::glue('{prefix}.tsv'))
  })
}

populate_rsid_from_snp <- function() {
  bim <- vroom::vroom('full.bim', col_names=F)
  full_rsid <- vroom::vroom('full_std.tsv.gz')
  matching <- match(bim$X2, full_rsid$SNP)
  bim$RSID <- full_rsid$RSID[matching]
  bim$X2[!is.na(bim$RSID)] <- bim$RSID[!is.na(bim$RSID)]
  bim <- dplyr::select(bim, -RSID)
  vroom::vroom_write(bim, 'full_rsid.bim', col_names=F)
}