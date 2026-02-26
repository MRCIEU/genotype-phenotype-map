# Pull exon boundaries for MANE transcript of each gene in build hg38

library(biomaRt)
library(dplyr)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Retrieve canonical transcripts
mane_transcripts <- getBM(
  attributes = c(
    "ensembl_gene_id", # Ensembl transcript ID
    "ensembl_transcript_id", # Ensembl transcript ID
    "chromosome_name", # Chromosome number
    "strand", # Strand information
    "transcript_mane_select", # MANE Select flag
    "transcript_is_canonical",
    "hgnc_symbol" # Gene symbol
  ),
  filters = c("transcript_is_canonical", "chromosome_name"),
  values = list(TRUE, 1:22), # Filter for canonical transcripts
  mart = ensembl
)

tx_ids <- unique(mane_transcripts$ensembl_transcript_id)

# Retrieve exon positons
mane_exons <- getBM(
  attributes = c(
    "ensembl_gene_id", # Ensembl transcript ID
    "ensembl_transcript_id", # Ensembl transcript ID
    "chromosome_name", # Chromosome number
    "strand", # Strand information
    "ensembl_exon_id", # Ensembl exon ID
    "rank",
    "exon_chrom_start", # MANE Select flag
    "exon_chrom_end" # Gene symbol
  ),
  filters = c("ensembl_transcript_id"),
  values = list(tx_ids), # Filter on canoncical tanscript ID
  mart = ensembl
)

exons_out <- mane_exons |>
  left_join(
    mane_transcripts |> filter(duplicated(ensembl_transcript_id) == FALSE),
    by = c("ensembl_gene_id", "ensembl_transcript_id", "chromosome_name", "strand")
  )

data.table::fwrite(exons_out, "/local-scratch/projects/genotype-phenotype-map/data/exon_ranges_hg38.tsv", sep = "\t")
