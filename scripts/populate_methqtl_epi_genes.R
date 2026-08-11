source("../pipeline_steps/constants.R")

manifest_file <- file.path(variant_annotation_dir, "MethylationEPIC_v-1-0_B4.csv.gz")
epi_dir <- "/local-scratch/data/hg38/GTEx-methQTL-v10"
gene_column <- "GencodeBasicV12_NAME"

main <- function() {
  probe_gene_map <- load_probe_gene_map()

  epi_files <- list.files(epi_dir, pattern = "\\.epi$", recursive = TRUE, full.names = TRUE)
  if (length(epi_files) == 0) {
    stop("No .epi files found in: ", epi_dir)
  }

  message("Updating ", length(epi_files), " .epi files in ", epi_dir)

  all_stats <- lapply(epi_files, update_epi_file, probe_gene_map = probe_gene_map)
  stats_df <- dplyr::bind_rows(all_stats)

  message(
    "Done. Updated ", nrow(stats_df), " files (",
    sum(stats_df$rows), " probes: ",
    sum(stats_df$annotated), " with gene, ",
    sum(stats_df$na), " NA)"
  )
  return(invisible(NULL))
}

first_manifest_gene <- function(gene_field) {
  if (length(gene_field) != 1 || is.na(gene_field) || gene_field == "" || gene_field == "NA") {
    return(NA_character_)
  }
  genes <- strsplit(gene_field, ";", fixed = TRUE)[[1]]
  genes <- trimws(genes)
  genes <- genes[genes != ""]
  if (length(genes) == 0) {
    return(NA_character_)
  }
  return(genes[1])
}

format_epi_gene <- function(gene) {
  return(ifelse(is.na(gene) | gene == "", "NA", gene))
}

load_probe_gene_map <- function() {
  if (!file.exists(manifest_file)) {
    stop("Manifest file not found: ", manifest_file)
  }

  manifest <- data.table::fread(
    manifest_file,
    select = c("IlmnID", gene_column, "UCSC_RefGene_Name"),
    showProgress = FALSE
  )

  manifest <- as.data.frame(manifest)
  manifest$gene <- vapply(
    seq_len(nrow(manifest)),
    function(i) {
      primary <- first_manifest_gene(manifest[[gene_column]][i])
      if (!is.na(primary)) {
        return(primary)
      }
      return(first_manifest_gene(manifest$UCSC_RefGene_Name[i]))
    },
    character(1)
  )

  probe_gene_map <- manifest |>
    dplyr::transmute(probe = IlmnID, gene = gene) |>
    dplyr::filter(!is.na(probe) & probe != "") |>
    dplyr::distinct(probe, .keep_all = TRUE)

  message(
    "Loaded probe-gene map: ",
    nrow(probe_gene_map),
    " probes (",
    sum(!is.na(probe_gene_map$gene)),
    " annotated, ",
    sum(is.na(probe_gene_map$gene)),
    " without gene)"
  )

  return(probe_gene_map)
}

update_epi_file <- function(epi_file, probe_gene_map) {
  epi <- vroom::vroom(epi_file, col_names = FALSE, show_col_types = FALSE)
  if (ncol(epi) < 5) {
    stop("Unexpected .epi format (expected at least 5 columns): ", epi_file)
  }

  mapped_genes <- probe_gene_map$gene[match(epi$X2, probe_gene_map$probe)]
  epi$X5 <- format_epi_gene(mapped_genes)
  vroom::vroom_write(epi, epi_file, col_names = FALSE)

  return(list(
    file = epi_file,
    rows = nrow(epi),
    annotated = sum(epi$X5 != "NA", na.rm = TRUE),
    na = sum(epi$X5 == "NA", na.rm = TRUE)
  ))
}

invisible(main())
