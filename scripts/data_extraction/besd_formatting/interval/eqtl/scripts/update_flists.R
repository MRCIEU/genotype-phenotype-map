library(data.table)
library(dplyr)

in_dir <- "flist_drafts/"
out_dir <- "updated_flists/"
ensembl_coords <- fread(
  "ensembl_coords.txt",
  select = c("ensembl_gene_id", "external_gene_name", "start_position", "end_position", "strand")
)

flists <- list.files(in_dir)

## loop over flist data
for (i in flists) {
  print(paste0("Updating", i))
  flist <- fread(paste0(in_dir, i))

  ## format gene names
  flist$ensembl_gene_id <- flist$ProbeID

  ## merge
  m <- merge(flist, ensembl_coords, by = "ensembl_gene_id")
  NROW(flist)
  NROW(m)

  ### double check probe position aligns to either start or end position
  NROW(m[which(m$ProbeBp == m$start_position), ])
  NROW(m[which(m$ProbeBp == m$end_position), ])

  NROW(m[which(m$ProbeBp == m$start_position), ]) + NROW(m[which(m$ProbeBp == m$end_position), ]) == NROW(m)
  ### 441 are different - why  - not an issue as probe position not used later

  ## not all transcripts map - keep those which are missing and put strand as "N", unknown
  missing <- flist[!which(flist$ensembl_gene_id %in% m$ensembl_gene_id), ] # 380


  m <- m[, c("Chr", "ProbeID", "ProbeBp", "external_gene_name", "strand", "PathOfEsd")]

  ## org
  missing_gene_symbol <- m[m$external_gene_name == "", ]
  ## missing_gene_symbol$external_gene_name <- missing_gene_symbol$ProbeID

  m1 <- m[!which(m$ProbeID %in% missing_gene_symbol$ProbeID), ]
  m2 <- m[which(m$ProbeID %in% missing_gene_symbol$ProbeID), ]
  m2$external_gene_name <- m2$ProbeID
  m3 <- missing

  m3 <- m3[, c(1:5)]
  m3$external_gene_name <- m3$ProbeID

  m3$strand <- "N"
  m3 <- m3[, c("Chr", "ProbeID", "ProbeBp", "external_gene_name", "strand", "PathOfEsd")]

  if ((NROW(m1) + NROW(m2) + NROW(m3)) == NROW(flist)) {
    print(paste0("All probes accounted for in ", i))
  } else {
    print(paste0("Some probes missing in ", i))
  }

  df <- rbind(m1, m2, m3)
  df <- df[order(df$Chr), ]


  df <- df %>% mutate(GeneticDistance = 0)
  df <- df %>% relocate(GeneticDistance, .after = 2)

  df$strand <- ifelse(df$strand == -1, "-", ifelse(df$strand == 1, "+", df$strand))
  df

  write.table(df, paste0(out_dir, i), col.names = T, row.names = F, quote = F, sep = "\t")
}

print(paste0("flist data update completed."))
