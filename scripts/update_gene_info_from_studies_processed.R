source("../pipeline_steps/constants.R")
suppressMessages(require(biomaRt))

# select which mart to use, in this case ensembl
mart <- useMart("ensembl")
mart <- useDataset("hsapiens_gene_ensembl", mart = mart)

gene_info <- vroom::vroom(glue::glue("{variant_annotation_dir}/gene_info.tsv"), show_col_types = F)
studies_processed <- vroom::vroom(glue::glue("{latest_results_dir}/studies_processed.tsv.gz"), show_col_types = F)

new_ensg_ids <- unique(studies_processed$ensg[!studies_processed$ensg %in% gene_info$ensembl_id]) |> na.omit()

# new_genes_by_name <- getBM(filters="external_gene_name",
#     attributes= c("ensembl_gene_id", "external_gene_name", "description","gene_biotype", "chromosome_name", "start_position", "end_position", "strand"), # nolint: line_length_linter.
#     values= new_ensg_ids, mart= mart
# ) |>
#     dplyr::mutate(chromosome_name = as.numeric(chromosome_name))

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
ensp_mapping <- getBM(
  filters = c("ensembl_gene_id", "transcript_is_canonical"),
  attributes = c("ensembl_gene_id", "ensembl_peptide_id"),
  values = list(new_ensg_ids, TRUE),
  mart = mart
) |>
  dplyr::filter(ensembl_peptide_id != "") |>
  dplyr::distinct(ensembl_gene_id, .keep_all = TRUE)

new_genes <- new_genes |>
  dplyr::left_join(ensp_mapping, by = "ensembl_gene_id")

new_genes <- new_genes |>
  dplyr::rename(
    chr = "chromosome_name",
    start = "start_position",
    stop = "end_position",
    ensembl_id = "ensembl_gene_id",
    gene_name = "external_gene_name"
  ) |>
  dplyr::mutate(source = sub(".*\\[(.*)\\]", "\\1", description)) |>
  dplyr::mutate(description = sub(" \\[.*\\]", "", description)) |>
  dplyr::distinct(ensembl_id, .keep_all = TRUE)

gene_info <- dplyr::bind_rows(gene_info, new_genes)

# Backfill ENSP for existing rows that are missing it
missing_ensp <- gene_info$ensembl_id[is.na(gene_info$ensembl_peptide_id)]
if (length(missing_ensp) > 0) {
  backfill_ensp <- getBM(
    filters = c("ensembl_gene_id", "transcript_is_canonical"),
    attributes = c("ensembl_gene_id", "ensembl_peptide_id"),
    values = list(missing_ensp, TRUE),
    mart = mart
  ) |>
    dplyr::filter(ensembl_peptide_id != "") |>
    dplyr::distinct(ensembl_gene_id, .keep_all = TRUE)

  idx <- match(backfill_ensp$ensembl_gene_id, gene_info$ensembl_id)
  gene_info$ensembl_peptide_id[idx[!is.na(idx)]] <- backfill_ensp$ensembl_peptide_id[!is.na(idx)]
}

vroom::vroom_write(gene_info, glue::glue("{variant_annotation_dir}/gene_info.tsv"))

# --- KEGG pathway mapping ---
# Get Entrez IDs for all genes in gene_info via biomaRt
all_ensg <- unique(gene_info$ensembl_id[!is.na(gene_info$ensembl_id)])
entrez_mapping <- getBM(
  filters = "ensembl_gene_id",
  attributes = c("ensembl_gene_id", "entrezgene_id"),
  values = all_ensg,
  mart = mart
) |>
  dplyr::filter(!is.na(entrezgene_id)) |>
  dplyr::mutate(entrezgene_id = as.character(entrezgene_id)) |>
  dplyr::distinct(ensembl_gene_id, entrezgene_id)

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
  dplyr::left_join(
    gene_info |> dplyr::select(ensembl_id, gene_name, ensembl_peptide_id),
    by = c("ensembl_gene_id" = "ensembl_id")
  ) |>
  dplyr::select(
    ensembl_id = ensembl_gene_id,
    gene_name,
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
