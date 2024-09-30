source('../../pipeline_steps/constants.R')
ld_regions <- vroom::vroom('../../pipeline_steps/data/ld_regions.tsv', show_col_types = F) |>
  dplyr::filter(ancestry == 'EUR')
ld_info <- construct_ld_block(ld_regions$ancestry, ld_regions$chr, ld_regions$start, ld_regions$stop)

mergefile <- tempfile()
write.table(ld_info$ld_reference_panel_prefix, file=mergefile, row=F, col=F, qu=F)
output <- glue::glue("{ld_reference_panel_dir}EUR/full")
glue::glue("plink1.9 --merge-list {mergefile} --make-bed --out {output} --keep-allele-order") |> system()


#then create a map of SNP to RSID