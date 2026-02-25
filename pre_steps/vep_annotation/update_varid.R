infile <- commandArgs(T)[1]

# Order alleles alphabetically
vars <- read.table(infile, col.names = c("chr", "pos1", "pos2", "alleles", "strand"))

vars$REF <- sub("\\/.*", "", vars$alleles)
vars$ALT <- sub(".*\\/", "", vars$alleles)

vars$name <- ifelse(vars$REF < vars$ALT,
  paste0(vars$chr, ":", vars$pos1, "_", vars$REF, "_", vars$ALT),
  paste0(vars$chr, ":", vars$pos1, "_", vars$ALT, "_", vars$REF)
)

vars <- vars[order(as.numeric(vars$chr), as.numeric(vars$pos1)), ]

write.table(vars[, c("chr", "pos1", "pos2", "alleles", "strand", "name")],
  file = infile,
  row.names = F, col.names = F, quote = F, sep = " "
)
