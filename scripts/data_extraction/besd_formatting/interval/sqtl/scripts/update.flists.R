library(data.table)
library(dplyr)

in_dir <- "/local-scratch/data/tmp_processes/interval/sqtl/flist_drafts/"

out_dir <- "/local-scratch/data/tmp_processes/interval/sqtl/updated_flists/"

ensembl_coords <- fread(
  "/home/rj18633/scratch/gp.map/data/parquet.to.besd/coordinates/ensembl_coords.txt",
  select = c("ensembl_gene_id", "external_gene_name", "start_position", "end_position", "strand")
)

ref <- fread("/local-scratch/data/interval/INTERVAL_sQTL_summary_statistics/INTERVAL_sQTL_phenotype_summary.tsv")

ref <- ref[, c(1, 2, 3)]
colnames(ref)[1] <- "ProbeID"

flists <- list.files(in_dir)

## loop over flist data

for (i in flists) {
  print(paste0("Updating ", i))
  flist <- fread(paste0(in_dir, i))

  ## emrge with ref to get gene name
  m <- merge(flist, ref)
  NROW(flist) == NROW(m)

  ##  keep first gene ref
  ## get strand orientation from probeID
  m$gene_id <- sapply(strsplit(as.character(m$gene_id), ","), `[`, 1)
  m$gene_name <- sapply(strsplit(as.character(m$gene_name), ","), `[`, 1)
  m$Orientation <- substr(sub(".*_", "", m$ProbeID), nchar(sub(".*_", "", m$ProbeID)), nchar(sub(".*_", "", m$ProbeID)))


  m1 <- m[, c("Chr", "ProbeID", "GeneticDistance", "ProbeBp", "gene_name", "Orientation", "PathOfEsd")]
  colnames(m1)[5] <- "Gene"
  write.table(m1, paste0(out_dir, i), col.names = T, row.names = F, quote = F, sep = "\t")
}

print(paste0("flist data update completed."))
