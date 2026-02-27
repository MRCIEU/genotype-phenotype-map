infile <- commandArgs(T)[1]

forvep <- data.table::fread(infile, col.names = c("chr", "pos1", "pos2", "alleles", "strand", "name"))
filename <- sub(".txt", "", infile)

forvep$REF <- sub("\\/.*", "", forvep$alleles)
forvep$ALT <- sub(".*\\/", "", forvep$alleles)

# Correct position for insertions
# An insertion (of any size) is indicated by start coordinate = end coordinate + 1

forvep$pos1 <- ifelse(nchar(forvep$ALT) > 1, forvep$pos1 + 1, forvep$pos1)
forvep$alleles <- ifelse(nchar(forvep$ALT) > 1, paste0("-/", substr(forvep$ALT, 2, nchar(forvep$ALT))), forvep$alleles)
forvep$REF <- ifelse(nchar(forvep$ALT) > 1, "-", forvep$REF)
forvep$ALT <- ifelse(nchar(forvep$ALT) > 1, substr(forvep$ALT, 2, nchar(forvep$ALT)), forvep$ALT)

# Correct position for deletions
forvep$pos1 <- ifelse(nchar(forvep$REF) > 1, forvep$pos1 + 1, forvep$pos1)
forvep$pos2 <- ifelse(nchar(forvep$REF) > 1, forvep$pos1 + nchar(forvep$REF) - 2, forvep$pos2)
forvep$alleles <- ifelse(nchar(forvep$REF) > 1, paste0(substr(forvep$REF, 2, nchar(forvep$REF)), "/-"), forvep$alleles)
forvep$ALT <- ifelse(nchar(forvep$REF) > 1, "-", forvep$ALT)
forvep$REF <- ifelse(nchar(forvep$REF) > 1, substr(forvep$REF, 2, nchar(forvep$REF)), forvep$REF)

write.table(forvep[, c("chr", "pos1", "pos2", "alleles", "strand", "name")],
  file = paste0(filename, "_posupdate.txt"),
  row.names = F, col.names = F, quote = F, sep = " "
)
