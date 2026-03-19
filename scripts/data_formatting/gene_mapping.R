library(biomaRt)
library(STRINGdb)
library(dplyr)

# Since the QTLs already provide the direct mapping from SNP to Gene, your task is much simpler: you need a Functional Gene Map.
# This is essentially a "lookup table" where each Gene ID is a key that reveals its biological "neighborhood" (PPI) and its "job description" (GO).

# To achieve this without bloating your R package, you can pre-compute a Gene-to-Function Master Table and a Gene-to-Gene Interaction Table.

# 1. The Data Structure You Need

# You should pre-compute two specific tables on your OCI server:

# Table A: The Functional Glossary (gene_metadata.tsv)
# | gene_id | symbol | go_terms (Processes) | kegg_pathways |
# | :--- | :--- | :--- | :--- |
# | ENSG00000169174 | PCSK9 | Cholesterol metabolism; LDL clearance | Metabolic pathways |

# Table B: The Interaction Map (gene_interactions.tsv)
# | gene_a | gene_b | confidence_score | source |
# | :--- | :--- | :--- | :--- |
# | Gene_1 | Gene_2 | 900 | STRING (Physical) |
#
# 2. The Pre-computation Script (One-Time Run)

# This script generates those two tables. It uses biomaRt to get the "Functional Terms" and STRINGdb to get the "Implicated Genes" (Interactions).

main <- function() {
  # 1. Get your unique genes from GPMap
  all_genes <- unique(results$coloc_groups$gene_id)

  # 2. Build the Glossary (Gene -> GO Terms)
  mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  glossary <- getBM(
    attributes = c("ensembl_gene_id", "external_gene_name", "name_1006"),
    filters = "ensembl_gene_id",
    values = all_genes,
    mart = mart
  )

  gene_glossary <- glossary %>%
    group_by(ensembl_gene_id) %>%
    summarise(
      symbol = first(external_gene_name),
      functional_terms = paste(unique(head(name_1006, 5)), collapse = "|"),
      .groups = 'drop'
    )

  # 3. Build the Interaction Map (Gene -> Implicated Genes)
  string_db <- STRINGdb$new(version="12.0", species=9606, score_threshold=700) # High confidence
  mapped <- string_db$map(data.frame(gene=all_genes), "gene", removeUnmappedRows = TRUE)
  interactions <- string_db$get_interactions(mapped$STRING_id)

  # Convert STRING IDs back to Ensembl for your users
  id_map <- mapped %>% select(STRING_id, gene)
  gene_interactions <- interactions %>%
    inner_join(id_map, by = c("from" = "STRING_id")) %>%
    rename(gene_a = gene) %>%
    inner_join(id_map, by = c("to" = "STRING_id")) %>%
    rename(gene_b = gene) %>%
    select(gene_a, gene_b, combined_score)

  # 4. Save these as lean TSVs for your R package
  write.table(gene_glossary, "gene_glossary.tsv", sep="\t", row.names=F, quote=F)
  write.table(gene_interactions, "gene_interactions.tsv", sep="\t", row.names=F, quote=F)
}

main()
