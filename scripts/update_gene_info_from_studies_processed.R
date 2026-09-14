source("../pipeline_steps/constants.R")
suppressMessages(require(biomaRt))

# select which mart to use, in this case ensembl
mart <- useMart("ensembl")
mart <- useDataset("hsapiens_gene_ensembl", mart = mart)

# Strip Ensembl gene ID version suffixes (e.g. ENSG00000149531.15 -> ENSG00000149531)
strip_ensembl_version <- function(x) {
  return(sub("\\.\\d+$", "", as.character(x)))
}

# ---------------------------------------------------------------------------
# Canonical gene-name <-> ENSG map from the HGNC complete set
# https://www.genenames.org/download/archive/
# ---------------------------------------------------------------------------
hgnc_complete <- vroom::vroom(
  file.path(variant_annotation_dir, "hgnc_complete_set.txt"),
  show_col_types = FALSE,
  col_types = vroom::cols(ensembl_gene_id = vroom::col_character())
) |>
  dplyr::mutate(ensembl_gene_id = strip_ensembl_version(ensembl_gene_id)) |>
  dplyr::filter(!is.na(ensembl_gene_id), ensembl_gene_id != "")

# Display symbol per gene: current approved HGNC symbol (falls back to any
# symbol for genes whose HGNC entry has been withdrawn)
hgnc_symbols <- hgnc_complete |>
  dplyr::filter(!is.na(symbol), symbol != "") |>
  dplyr::group_by(ensembl_gene_id) |>
  dplyr::arrange(dplyr::desc(status == "Approved"), symbol) |>
  dplyr::slice(1) |>
  dplyr::ungroup() |>
  dplyr::select(ensembl_gene_id, gene_name_hgnc = symbol, entrez_id) |>
  dplyr::mutate(entrez_id = as.character(entrez_id))

# Comma-separated alias string per gene from HGNC previous + alias symbols,
# so old/alternative names keep resolving in searches
format_hgnc_aliases <- function(alias_symbol, prev_symbol, symbol, ensg) {
  symbol <- if (length(symbol) == 0 || is.na(symbol)) "" else symbol
  ensg <- if (length(ensg) == 0 || is.na(ensg)) "" else ensg
  parts <- trimws(unlist(strsplit(
    paste(
      ifelse(is.na(prev_symbol), "", prev_symbol),
      ifelse(is.na(alias_symbol), "", alias_symbol),
      sep = "|"
    ),
    "[|,]"
  ), use.names = FALSE))
  parts <- parts[parts != "" & !is.na(parts) & parts != symbol & parts != ensg]
  parts <- unique(parts)
  if (length(parts) == 0) {
    return(NA_character_)
  }
  return(paste(parts, collapse = ", "))
}

hgnc_aliases <- hgnc_complete |>
  dplyr::filter(
    (!is.na(alias_symbol) & alias_symbol != "") |
      (!is.na(prev_symbol) & prev_symbol != "")
  ) |>
  dplyr::group_by(ensembl_gene_id) |>
  dplyr::summarise(
    gene_alias_hgnc = format_hgnc_aliases(alias_symbol, prev_symbol, symbol[1], ensembl_gene_id[1]),
    .groups = "drop"
  )

hgnc_map <- hgnc_symbols |>
  dplyr::left_join(hgnc_aliases, by = "ensembl_gene_id")

# Resolve display name: HGNC approved symbol > biomaRt/stored symbol > ENSG id
resolve_gene_name <- function(ensembl_id, external_gene_name, gene_name_hgnc) {
  return(dplyr::case_when(
    !is.na(gene_name_hgnc) & gene_name_hgnc != "" ~ gene_name_hgnc,
    !is.na(external_gene_name) & external_gene_name != "" & external_gene_name != ensembl_id ~ external_gene_name,
    TRUE ~ ensembl_id
  ))
}

# ---------------------------------------------------------------------------
# Existing map + newly-observed genes from processed studies
# ---------------------------------------------------------------------------
gene_info <- vroom::vroom(glue::glue("{variant_annotation_dir}/gene_info.tsv"), show_col_types = F)
studies_processed <- vroom::vroom(glue::glue("{latest_results_dir}/studies_processed.tsv.gz"), show_col_types = F)

# Compare ENSG ids without version suffixes (change #3)
studies_ensg <- strip_ensembl_version(studies_processed$ensg)
gene_info <- gene_info |>
  dplyr::mutate(ensembl_id = strip_ensembl_version(ensembl_id))

new_ensg_ids <- unique(studies_ensg[!studies_ensg %in% gene_info$ensembl_id]) |>
  na.omit() |>
  as.character()
new_ensg_ids <- new_ensg_ids[grepl("^ENSG", new_ensg_ids)]

new_genes <- getBM(
  filters = "ensembl_gene_id",
  attributes = c(
    "ensembl_gene_id",
    "external_gene_name",
    "description",
    "gene_biotype",
    "chromosome_name",
    "start_position",
    "end_position",
    "strand"
  ),
  values = new_ensg_ids, mart = mart
)

# Fetch canonical ENSP (protein) IDs for the new genes
new_ensp_mapping <- getBM(
  filters = c("ensembl_gene_id", "transcript_is_canonical"),
  attributes = c("ensembl_gene_id", "ensembl_peptide_id"),
  values = list(new_ensg_ids, TRUE),
  mart = mart
) |>
  dplyr::filter(ensembl_peptide_id != "") |>
  dplyr::distinct(ensembl_gene_id, .keep_all = TRUE)

new_genes <- new_genes |>
  dplyr::rename(
    chr = "chromosome_name",
    start = "start_position",
    stop = "end_position",
    ensembl_id = "ensembl_gene_id"
  ) |>
  dplyr::mutate(source = sub(".*\\[(.*)\\]", "\\1", description)) |>
  dplyr::mutate(description = sub(" \\[.*\\]", "", description)) |>
  dplyr::distinct(ensembl_id, .keep_all = TRUE) |>
  dplyr::left_join(hgnc_map, by = c("ensembl_id" = "ensembl_gene_id")) |>
  dplyr::mutate(
    gene_name = resolve_gene_name(ensembl_id, external_gene_name, gene_name_hgnc),
    gene_alias = dplyr::coalesce(gene_alias_hgnc, "")
  ) |>
  dplyr::select(-gene_name_hgnc, -gene_alias_hgnc, -entrez_id, -external_gene_name) |>
  dplyr::rename(gene = gene_name)

# Re-sync existing rows to current HGNC symbols (and aliases), fixing rows
# whose display name fell back to the ENSG id (change #4)
gene_info <- gene_info |>
  dplyr::left_join(hgnc_map, by = c("ensembl_id" = "ensembl_gene_id")) |>
  dplyr::mutate(
    gene_name = resolve_gene_name(ensembl_id, gene, gene_name_hgnc),
    gene_alias = dplyr::coalesce(gene_alias_hgnc, "")
  ) |>
  dplyr::select(-gene_name_hgnc, -gene_alias_hgnc, -entrez_id) |>
  dplyr::rename(gene = gene_name)

gene_info <- dplyr::bind_rows(gene_info, new_genes) |>
  dplyr::distinct(ensembl_id, .keep_all = TRUE) |>
  dplyr::select(ensembl_id, gene, gene_alias, description, gene_biotype, chr, start, stop, strand, source)

# Backfill canonical ENSP for any gene missing it (used by the KEGG output)
ensp_map <- new_ensp_mapping
missing_ensp <- setdiff(gene_info$ensembl_id, ensp_map$ensembl_gene_id)
if (length(missing_ensp) > 0) {
  backfill_ensp <- getBM(
    filters = c("ensembl_gene_id", "transcript_is_canonical"),
    attributes = c("ensembl_gene_id", "ensembl_peptide_id"),
    values = list(missing_ensp, TRUE),
    mart = mart
  ) |>
    dplyr::filter(ensembl_peptide_id != "") |>
    dplyr::distinct(ensembl_gene_id, .keep_all = TRUE)

  ensp_map <- dplyr::bind_rows(ensp_map, backfill_ensp) |>
    dplyr::distinct(ensembl_gene_id, .keep_all = TRUE)
}

vroom::vroom_write(gene_info, glue::glue("{variant_annotation_dir}/gene_info.tsv"))

# --- KEGG pathway mapping ---
# Entrez IDs from HGNC where available, falling back to biomaRt for genes
# that are absent from HGNC
all_ensg <- unique(gene_info$ensembl_id[!is.na(gene_info$ensembl_id)])
entrez_mapping <- hgnc_map |>
  dplyr::select(ensembl_gene_id, entrezgene_id = entrez_id) |>
  dplyr::filter(!is.na(entrezgene_id), entrezgene_id != "", entrezgene_id != "NA") |>
  dplyr::mutate(entrezgene_id = as.character(entrezgene_id)) |>
  dplyr::distinct(ensembl_gene_id, entrezgene_id)

missing_entrez <- setdiff(all_ensg, entrez_mapping$ensembl_gene_id)
if (length(missing_entrez) > 0) {
  biomaRt_entrez <- getBM(
    filters = "ensembl_gene_id",
    attributes = c("ensembl_gene_id", "entrezgene_id"),
    values = missing_entrez,
    mart = mart
  ) |>
    dplyr::filter(!is.na(entrezgene_id)) |>
    dplyr::mutate(entrezgene_id = as.character(entrezgene_id)) |>
    dplyr::distinct(ensembl_gene_id, entrezgene_id)

  entrez_mapping <- dplyr::bind_rows(entrez_mapping, biomaRt_entrez)
}

# Query KEGG REST API: gene → pathway mapping
kegg_gene_to_pathway <- data.table::fread(
  "https://rest.kegg.jp/link/hsa/pathway",
  header = FALSE, col.names = c("pathway_id", "kegg_gene_id")
) |>
  dplyr::mutate(
    pathway_id = sub("^path:", "", pathway_id),
    entrezgene_id = sub("^hsa:", "", kegg_gene_id)
  ) |>
  dplyr::select(pathway_id, entrezgene_id)

# Query KEGG REST API: pathway ID → pathway name
kegg_pathway_names <- data.table::fread(
  "https://rest.kegg.jp/list/pathway/hsa",
  header = FALSE, col.names = c("pathway_id", "pathway_name")
) |>
  dplyr::mutate(pathway_name = sub(" - Homo sapiens \\(human\\)", "", pathway_name))

# Join: ENSG → Entrez → KEGG pathway → pathway name
kegg_pathways <- entrez_mapping |>
  dplyr::inner_join(kegg_gene_to_pathway, by = "entrezgene_id") |>
  dplyr::inner_join(kegg_pathway_names, by = "pathway_id") |>
  dplyr::left_join(ensp_map, by = "ensembl_gene_id") |>
  dplyr::left_join(
    gene_info |> dplyr::select(ensembl_id, gene),
    by = c("ensembl_gene_id" = "ensembl_id")
  ) |>
  dplyr::select(
    ensembl_id = ensembl_gene_id,
    gene_name = gene,
    ensembl_peptide_id,
    entrezgene_id,
    pathway_id,
    pathway_name
  ) |>
  dplyr::distinct() |>
  dplyr::arrange(pathway_id, ensembl_id)

message(glue::glue(
  "KEGG: {length(unique(kegg_pathways$pathway_id))} pathways, ",
  "{length(unique(kegg_pathways$ensembl_id))} genes"
))

vroom::vroom_write(
  kegg_pathways,
  glue::glue("{variant_annotation_dir}/kegg_pathways.tsv")
)
