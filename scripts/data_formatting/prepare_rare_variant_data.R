backman <- function() {
  backman_metadata_file <- "/local-scratch/data/ukb-seq/downloads/backmanexwas/ukb-wes-bm-metadata.tsv"
  new_backman_dir <- "/local-scratch/data/hg38/rare/backman/"

  backman_metadata <- vroom::vroom(backman_metadata_file, show_col_types = F) |>
    dplyr::select(-extracted_location)

  new_metadata <- apply(backman_metadata, 1, function(study) {
    rv <- vroom::vroom(study[["study_location"]], show_col_types = F)
    new_basename <- basename(sub("\\..*$", "", study[["study_location"]]))
    new_file_name <- paste0(new_backman_dir, new_basename, ".tsv.gz")
    print(new_file_name)
    study[["study_location"]] <- new_file_name
    if (file.exists(new_file_name)) {
      return(study)
    }

    if ("beta" %in% colnames(rv)) {
      rv <- dplyr::rename(
        rv,
        CHR = chromosome,
        BP = base_pair_location,
        P = p_value,
        EA = effect_allele,
        OA = other_allele,
        trait = Trait,
        BETA = beta,
        EAF = effect_allele_frequency,
        SE = standard_error
      )
    } else if ("odds_ratio" %in% colnames(rv)) {
      rv <- dplyr::rename(
        rv,
        CHR = chromosome,
        BP = base_pair_location,
        P = p_value,
        EA = effect_allele,
        OA = other_allele,
        trait = Trait,
        OR = odds_ratio,
        CI_LOWER = ci_lower,
        CI_UPPER = ci_upper,
        EAF = effect_allele_frequency
      )
    }
    vroom::vroom_write(rv, new_file_name)
    return(study)
  })
  new_metadata <- dplyr::bind_rows(as.data.frame(t(new_metadata)))
  vroom::vroom_write(new_metadata, glue::glue("{new_backman_dir}/metadata.tsv"))
  return()
}

genebass <- function() {
  genebass_metadata_file <- "/local-scratch/data/ukb-seq/downloads/genebass/rare/ukb-wes-gb-metadata.tsv"
  new_genebass_dir <- "/local-scratch/data/hg38/rare/genebass/"

  genebass_metadata <- vroom::vroom(genebass_metadata_file, show_col_types = F) |>
    dplyr::select(-extracted_location)

  new_metadata <- apply(genebass_metadata, 1, function(study) {
    original_study_location <- study[["study_location"]]
    new_basename <- basename(sub("\\..*$", "", study[["study_location"]]))
    new_file_name <- paste0(new_genebass_dir, new_basename, ".tsv.gz")
    study[["study_location"]] <- new_file_name
    if (file.exists(new_file_name)) {
      return(study)
    }

    rv <- vroom::vroom(original_study_location, show_col_types = F)
    print(new_file_name)
    rv <- tidyr::separate(rv, col = markerID, into = c("CHR", "BP", "OA", "EA"), sep = "[:_/]", remove = F) |>
      dplyr::mutate(CHR = sub("chr", "", CHR)) |>
      dplyr::rename(P = Pvalue, EAF = AF, GENE = "gene", ANNOTATION = "annotation")
    vroom::vroom_write(rv, new_file_name)
    return(study)
  })
  new_metadata <- dplyr::bind_rows(as.data.frame(t(new_metadata)))
  vroom::vroom_write(new_metadata, glue::glue("{new_genebass_dir}/metadata.tsv"))
  return()
}

az <- function() {
  az_metadata_file <- "/local-scratch/data/ukb-seq/downloads/azexwas/ukb-wes-az-metadata.tsv"
  new_az_dir <- "/local-scratch/data/hg38/rare/azexwas/"

  az_metadata <- vroom::vroom(az_metadata_file, show_col_types = F) |>
    dplyr::select(-extracted_location)

  new_metadata <- apply(az_metadata, 1, function(study) {
    rv <- vroom::vroom(study[["study_location"]], show_col_types = F)
    new_basename <- basename(sub("\\..*$", "", study[["study_location"]]))
    new_file_name <- paste0(new_az_dir, new_basename, ".tsv.gz")
    study[["study_location"]] <- new_file_name
    if (file.exists(new_file_name)) {
      return(study)
    }

    if ("Effect size" %in% colnames(rv)) {
      rv <- tidyr::separate(rv, col = Variant, into = c("CHR", "BP", "OA", "EA"), sep = "-", remove = F) |>
        dplyr::rename(
          P = "p-value",
          EAF = "AAF",
          BETA = "Effect size",
          SE = "Effect size standard error",
          GENE = "Gene",
          ANNOTATION = "Consequence type"
        ) |>
        dplyr::mutate(Phenotype = sub(".*#", "", Phenotype))
    } else if ("Odds ratio" %in% colnames(rv)) {
      rv <- tidyr::separate(rv, col = Variant, into = c("CHR", "BP", "OA", "EA"), sep = "-", remove = F) |>
        dplyr::rename(
          P = "p-value",
          EAF = "Control AAF",
          OR = "Odds ratio",
          CI_LOWER = "Odds ratio LCI",
          CI_UPPER = "Odds ratio UCI",
          GENE = "Gene",
          ANNOTATION = "Consequence type"
        ) |>
        dplyr::mutate(Phenotype = sub(".*#", "", Phenotype))
    }
    vroom::vroom_write(rv, new_file_name)
    return(study)
  })
  new_metadata <- dplyr::bind_rows(as.data.frame(t(new_metadata)))
  vroom::vroom_write(new_metadata, glue::glue("{new_az_dir}/metadata.tsv"))
  return()
}

backman()
genebass()
az()
