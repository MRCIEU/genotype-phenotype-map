source("../../pipeline_steps/constants.R")

merge_all_bfiles <- function(ld_info) {
  mergefile <- withr::local_tempfile()
  write.table(ld_info$ld_reference_panel_prefix, file = mergefile, row = F, col = F, qu = F)
  output <- glue::glue("{ld_reference_panel_dir}EUR/new_full")
  glue::glue("plink1.9 --merge-list {mergefile} --make-bed --out {output} --keep-allele-order") |> system()
  return()
}

create_tsv_from_afreq <- function(prefix) {
  freq <- vroom::vroom(glue::glue("{prefix}.afreq"), show_col_types = F) |>
    dplyr::rename(CHR = "#CHROM", SNP = "ID", EA = "ALT", OA = "REF", EAF = "ALT_FREQS") |>
    dplyr::mutate(BP = as.numeric(sub(".*:(\\d+)_.*", "\\1", SNP))) |>
    dplyr::select(-`PROVISIONAL_REF?`, -OBS_CT)

  vroom::vroom_write(freq, glue::glue("{prefix}.tsv"))
  return()
}

populate_rsid_from_snp <- function() {
  bim <- vroom::vroom(glue::glue("{ld_reference_panel_dir}/EUR/new_full.bim"), col_names = F)
  full_rsid <- vroom::vroom(glue::glue("{ld_reference_panel_dir}/EUR/chrbp_rsid_map.tsv.gz"))
  matching <- match(bim$X2, full_rsid$SNP)
  bim$RSID <- full_rsid$RSID[matching]
  bim$X2[!is.na(bim$RSID)] <- bim$RSID[!is.na(bim$RSID)]
  bim <- dplyr::select(bim, -RSID)
  vroom::vroom_write(bim, glue::glue("{ld_reference_panel_dir}/EUR/new_full_rsid.bim"), col_names = F)
  file.copy(
    glue::glue("{ld_reference_panel_dir}/EUR/new_full.bed"),
    glue::glue("{ld_reference_panel_dir}/EUR/new_full_rsid.bed")
  )
  file.copy(
    glue::glue("{ld_reference_panel_dir}/EUR/new_full.fam"),
    glue::glue("{ld_reference_panel_dir}/EUR/new_full_rsid.fam")
  )
  return()
}

ld_blocks <- vroom::vroom("../../pipeline_steps/data/ld_blocks.tsv", show_col_types = F) |>
  dplyr::filter(ancestry == "EUR")
ld_info <- construct_ld_block(ld_blocks$ancestry, ld_blocks$chr, ld_blocks$start, ld_blocks$stop)

message("creating tsv files")
result <- lapply(ld_info$ld_reference_panel_prefix, function(ld_reference_panel_prefix) {
  if (file.exists(glue::glue("{ld_reference_panel_prefix}.tsv"))) {
    return()
  }

  print(ld_reference_panel_prefix)
  create_tsv_from_afreq(ld_reference_panel_prefix)
  return()
})

message("merging bfiles")
merge_all_bfiles(ld_info)

message("populating rsid bfiles")
for_snp <- vroom::vroom(glue::glue("{ld_reference_panel_dir}EUR/new_full.bim"), show_col_types = F, col_names = F) |>
  dplyr::mutate(BETA = 0, SE = 0.1, P = 1, EAF = 0.1) |>
  dplyr::rename(CHR = X1, SNP = X2, BP = X4, EA = X5, OA = X6)
vroom::vroom_write(for_snp, glue::glue("{ld_reference_panel_dir}EUR/new_full.tsv.gz"))

# note: must have chrbp_rsid_map.tsv.gz in place to do this step... Use GeneHackman to populate RSIDs if you don't have them handy # nolint: line_length_linter.
message("populating rsid bfiles")
populate_rsid_from_snp()
