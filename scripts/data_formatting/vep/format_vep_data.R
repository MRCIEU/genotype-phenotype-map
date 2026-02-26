source("../../../pipeline_steps/constants.R")
source("../../../pipeline_steps/common_extraction_functions.R")

main <- function() {
  vep_annotations_file <- glue::glue("{variant_annotation_dir}/vep_annotations_hg38.tsv.gz")
  if (file.exists(vep_annotations_file)) {
    existing_vep <- vroom::vroom(vep_annotations_file)
  } else {
    existing_vep <- data.frame()
  }
  genebass_vep <- vroom::vroom(
    glue::glue("{variant_annotation_dir}/genebass_vep_rare_variantannotations_hg38_altered.txt")
  ) |>
    dplyr::filter(!Uploaded_variation %in% existing_vep$snp)
  backman_vep <- vroom::vroom(
    glue::glue("{variant_annotation_dir}/backman_vep_rare_variantannotations_hg38_altered.txt")
  ) |>
    dplyr::filter(!Uploaded_variation %in% existing_vep$snp)
  azphewas_vep <- vroom::vroom(
    glue::glue("{variant_annotation_dir}/azphewas_vep_rare_variantannotations_hg38_altered.txt")
  ) |>
    dplyr::filter(!Uploaded_variation %in% existing_vep$snp)
  common_vep <- vroom::vroom(glue::glue("{variant_annotation_dir}/vep_variantannotations_hg38_altered.txt")) |>
    dplyr::filter(!Uploaded_variation %in% existing_vep$snp)

  extract_columns <- function(vep) {
    return(vep |>
      tidyr::separate(col = "Uploaded_variation", into = c("CHR", "BP", "EA", "OA"), sep = "[:_]", remove = F) |>
      dplyr::mutate(CHR = as.numeric(CHR)) |>
      dplyr::filter(!is.na(CHR) & CHR != 23) |>
      dplyr::mutate(
        impact = stringr::str_extract(Extra, "(?<=IMPACT=)[^;]+"),
        symbol = stringr::str_extract(Extra, "(?<=SYMBOL=)[^;]+"),
        biotype = stringr::str_extract(Extra, "(?<=BIOTYPE=)[^;]+"),
        strand = stringr::str_extract(Extra, "(?<=STRAND=)[^;]+"),
        canonical = stringr::str_extract(Extra, "(?<=CANONICAL=)[^;]+"),
        ref_allele = stringr::str_extract(Extra, "(?<=REF_ALLELE=)[^;]+"),
        all_af = stringr::str_extract(Extra, "(?<=AF=)[^;]+"),
        eur_af = stringr::str_extract(Extra, "(?<=EUR_AF=)[^;]+"),
        eas_af = stringr::str_extract(Extra, "(?<=EAS_AF=)[^;]+"),
        amr_af = stringr::str_extract(Extra, "(?<=AMR_AF=)[^;]+"),
        afr_af = stringr::str_extract(Extra, "(?<=AFR_AF=)[^;]+"),
        sas_af = stringr::str_extract(Extra, "(?<=SAS_AF=)[^;]+"),
        flipped = EA == ref_allele,
        EA_new = ifelse(flipped, OA, EA),
        OA_new = ifelse(flipped, EA, OA),
      ) |>
      dplyr::rename(rsid = Existing_variation, snp = Uploaded_variation) |>
      dplyr::mutate(display_snp = paste0(CHR, ":", BP, " ", OA_new, "/", EA_new)) |>
      dplyr::select(-Location, -Allele, -Feature, -Extra, -EA, -OA) |>
      dplyr::rename(EA = EA_new, OA = OA_new) |>
      dplyr::rename_with(tolower))
  }

  all_new_vep <- rbind(common_vep, genebass_vep, backman_vep, azphewas_vep)

  all_new_vep <- all_new_vep |>
    dplyr::distinct(Uploaded_variation, .keep_all = T) |>
    extract_columns()

  duplicates <- all_new_vep[duplicated(all_new_vep$snp) | duplicated(all_new_vep$display_snp), ]
  if (nrow(duplicates) > 0) {
    message(dplyr::select(duplicates, snp, display_snp))
    stop("ERROR: Found ", nrow(duplicates), " rows with duplicate snp or display_snp values.")
  }

  message(glue::glue("Formatting {nrow(all_new_vep)} new VEP annotations"))
  all_vep <- rbind(existing_vep, all_new_vep)
  vroom::vroom_write(all_vep, vep_annotations_file)
  return()
}

main()
