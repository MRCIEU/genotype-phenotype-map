suppressMessages(require(biomaRt))

# select which mart to use, in this case ensembl
mart <- useMart("ensembl")
mart <- useDataset("hsapiens_gene_ensembl", mart = mart)

gene_info <- vroom::vroom('{variant_annotation_dir}/gene_info.tsv', show_col_types = F)
studies_processed <- vroom::vroom('{results_dir}/latest/studies_processed.tsv', show_col_types = F)

new_ensg_ids <- unique(studies_processed$ensg[!studies_processed$ensg %in% gene_info$ensembl_gene_id])

new_genes <- getBM(filters="ensembl_gene_id",
    attributes= c("ensembl_gene_id", "external_gene_name", "description","gene_biotype", "chromosome_name", "start_position", "end_position", "strand"), 
    values= new_ensg_ids, mart= mart
)

new_genes <- new_genes |>
    dplyr::rename(chr='chromosome_name',
        start='start_position',
        stop='end_position',
        ensembl_id='ensembl_gene_id',
        gene_name='external_gene_name'
    ) |>
    dplyr::mutate(source = sub('.*\\[(.*)\\]', '\\1', description)) |>
    dplyr::mutate(description = sub(' \\[.*\\]', '', description))

gene_info <- dplyr::bind_rows(gene_info, new_genes)
vroom::vroom_write(gene_info, '{variant_annotation_dir}/gene_info.tsv', show_col_types = F)

