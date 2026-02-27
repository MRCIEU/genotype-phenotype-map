library(data.table)
library(dplyr)
library(parallel)

in_dir <- "INTERVAL_sQTL_summary_statistics/"
flist_dir <- "sqtl/flist_drafts/"
out_dir <- "esds/"

chrs <- list.files(in_dir)
chrs <- chrs[1:22]

# Define the per-file processing function
process_file <- function(data) {
  print(paste("Processing:", data))

  chromosome <- sub("\\.tsv$", "", data)
  file_path <- file.path(in_dir, data)

  file <- fread(file_path, select = c(1:4, 7:13))
  file$SNP <- paste(file$chr, file$pos_b38, sep = ":")

  colnames(file) <- c("phenotype_id", "variant_id", "splicemid_distance", "Freq", "pval_nominal",
                      "b", "se", "chr", "pos_b38", "A1", "A2", "SNP")

  # Standardise alleles
  standardise_alleles <- function(qtl) {
    qtl$Original <- TRUE
    to_flip <- qtl$A1 > qtl$A2
    if (any(to_flip)) {
      qtl$Freq[to_flip] <- 1 - qtl$Freq[to_flip]
      qtl$b[to_flip] <- -1 * qtl$b[to_flip]
      temp <- qtl$A2[to_flip]
      qtl$A2[to_flip] <- qtl$A1[to_flip]
      qtl$A1[to_flip] <- temp
      qtl$Original[to_flip] <- FALSE
    }
    qtl$SNP <- paste(qtl$chr, qtl$pos_b38, qtl$A1, qtl$A2, sep = "_")
    return(qtl)
  }

  dat_flipped <- standardise_alleles(file)

  dat_clean <- dat_flipped[,
    c("chr", "SNP", "pos_b38", "A1", "A2", "Freq", "b", "se", "pval_nominal", "splicemid_distance", "phenotype_id")
  ]
  colnames(dat_clean) <- c("Chr", "SNP", "Bp", "A1", "A2", "Freq",
                           "Beta", "se", "p", "splicemid_distance", "probe")
  dat_clean$Bp <- as.numeric(dat_clean$Bp)
  setDT(dat_clean)

  # Create esd output directory
  esd_dir <- file.path(out_dir, chromosome)
  if (!dir.exists(esd_dir)) dir.create(esd_dir, recursive = TRUE)

  # Create flist metadata
  flist_dat <- dat_clean[, .(Chr, probe, Bp, splicemid_distance)]
  flist_dat[, ProbeBp := Bp - splicemid_distance]
  flist_dat <- unique(flist_dat[, .(Chr, probe, ProbeBp)])
  flist_dat[, PathOfEsd := file.path(esd_dir, paste0(probe, ".esd"))]
  colnames(flist_dat) <- c("Chr", "ProbeID", "ProbeBp", "PathOfEsd")
  flist_dat$GeneticDistance <- 0

  write.table(
    unique(flist_dat),
    file.path(flist_dir, paste0("flist_", chromosome, ".draft.txt")),
    col.names = TRUE,
    row.names = FALSE,
    sep = "\t",
    quote = FALSE
  )

  # Split by gene and write ESDs
  genes <- split(dat_clean, by = "probe", keep.by = TRUE)
  genes <- lapply(genes, function(df) df[, 1:9])
  for (name in names(genes)) {
    out_file <- file.path(esd_dir, paste0(name, ".esd"))
    write.table(genes[[name]], out_file, col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
  }

  return(paste("Completed:", data))
}

# run in parallel
results <- mclapply(chrs, process_file, mc.cores = 12)
print("completed formatting esds! Next update flist drafts.")
