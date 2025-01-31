vep <- vep |> 
  tidyr::separate(col = "Uploaded_variation", into = c("CHR", "BP", "EA", "OA"), sep = "[:_]", remove = F) |>
  dplyr::mutate(impact = stringr::str_extract(vep$Extra, "(?<=IMPACT=)[^;]+"),
    symbol = stringr::str_extract(vep$Extra, "(?<=SYMBOL=)[^;]+"),
    biotype = stringr::str_extract(vep$Extra, "(?<=BIOTYPE=)[^;]+"),
    strand = stringr::str_extract(vep$Extra, "(?<=STRAND=)[^;]+"),
    canonical = stringr::str_extract(vep$Extra, "(?<=CANONICAL=)[^;]+")
    canonical = stringr::str_extract(vep$Extra, "(?<=CANONICAL=)[^;]+")
    ALL_AF = stringr::str_extract(vep$Extra, "(?<=AF=)[^;]+")
    EUR_AF = stringr::str_extract(vep$Extra, "(?<=EUR_AF=)[^;]+")
    EAS_AF = stringr::str_extract(vep$Extra, "(?<=EAS_AF=)[^;]+")
    AMR_AF = stringr::str_extract(vep$Extra, "(?<=AMR_AF=)[^;]+")
    AFR_AF = stringr::str_extract(vep$Extra, "(?<=AFR_AF=)[^;]+")
    SAS_AF = stringr::str_extract(vep$Extra, "(?<=SAS_AF=)[^;]+")
  ) |>
  dplyr::rename(SNP = Uploaded_variation, RSID = Existing_variation) |>
  dplyr::select(-Location, -Allele, -Feature, -Extra)

