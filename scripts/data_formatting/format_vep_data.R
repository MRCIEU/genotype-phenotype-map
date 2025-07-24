source('../../pipeline_steps/constants.R')
vep <- vroom::vroom(glue::glue('{variant_annotation_dir}/vep_variantannotations_hg38_altered.txt'))
rare_vep <- vroom::vroom(glue::glue('{variant_annotation_dir}/vep_rare_variantannotations_hg38_altered.txt'))

extract_columns <- function(vep) {
  vep |>
    tidyr::separate(col = "Uploaded_variation", into = c("CHR", "BP", "EA", "OA"), sep = "[:_]", remove = F) |>
    dplyr::mutate(impact = stringr::str_extract(vep$Extra, "(?<=IMPACT=)[^;]+"),
      symbol = stringr::str_extract(vep$Extra, "(?<=SYMBOL=)[^;]+"),
      biotype = stringr::str_extract(vep$Extra, "(?<=BIOTYPE=)[^;]+"),
      strand = stringr::str_extract(vep$Extra, "(?<=STRAND=)[^;]+"),
      canonical = stringr::str_extract(vep$Extra, "(?<=CANONICAL=)[^;]+"),
      ref_allele = stringr::str_extract(vep$Extra, "(?<=REF_ALLELE=)[^;]+"),
      all_af = stringr::str_extract(vep$Extra, "(?<=AF=)[^;]+"),
      eur_af = stringr::str_extract(vep$Extra, "(?<=EUR_AF=)[^;]+"),
      eas_af = stringr::str_extract(vep$Extra, "(?<=EAS_AF=)[^;]+"),
      amr_af = stringr::str_extract(vep$Extra, "(?<=AMR_AF=)[^;]+"),
      afr_af = stringr::str_extract(vep$Extra, "(?<=AFR_AF=)[^;]+"),
      sas_af = stringr::str_extract(vep$Extra, "(?<=SAS_AF=)[^;]+"),
    ) |>
    dplyr::rename(SNP = Uploaded_variation, RSID = Existing_variation) |>
    dplyr::mutate(flipped = ifelse(EA == ref_allele, F, T)) |>
    dplyr::mutate(display_snp = ifelse(flipped,
      paste0(CHR, ':', BP, ' ', OA, '/', EA),
      paste0(CHR, ':', BP, ' ', EA, '/', OA))
    ) |>
    dplyr::select(-Location, -Allele, -Feature, -Extra) |>
    dplyr::rename_with(tolower)
}

vep <- extract_columns(vep)
rare_vep <- extract_columns(rare_vep)

all_vep <- rbind(vep, rare_vep)

vroom::vroom_write(all_vep, glue::glue('{variant_annotation_dir}/vep_annotations_hg38.tsv.gz'))