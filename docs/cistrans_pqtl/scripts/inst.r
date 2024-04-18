# ACOX1 ieu-b-38 SBP
# ADAMTS8 ieu-b-38 SBP
# BMP6 ieu-b-38 SBP
# CCN3 ieu-b-38 SBP
# CD5 ieu-b-18 MS
# CD58 ieu-b-18 MS
# CD40 ieu-b-18 MS
# PCSK9 ieu-a-300 LDL

library(dplyr)
library(here)
library(parallel)
library(data.table)

lookup_txt <- function(fn, pos) {
    tf <- tempfile()
    tf2 <- tempfile()
    write.table(unique(pos), file=tf, row=F, col=F, qu=F)
    cmd <- glue("zgrep -wf {tf} {fn} > {tf2}")
    system(cmd)
    fread(tf2)
}

lookup_txt2 <- function(fn, pos) {
    a <- fread(fn)
    subset(a, ID %in% pos)
}

extract_tar <- function(variants, tarfile) {
    cmd <- paste0("tar xf '", tarfile, "' -C .")
    system(cmd)
    dirname <- gsub(".tar", "", basename(tarfile))
    filelist <- list.files(dirname, full.names=TRUE)
    lapply(filelist, \(x) {
        message(x)
        lookup_txt2(x, variants)
    }) %>% bind_rows()
}

tarfiles <- list.files("/projects/MRC-IEU/research/projects/icep2/wp1/028/working/data/ukbb_pqtl", full.names=TRUE) %>% grep(".tar", .,  value=TRUE)
p <- c("ACOX1", "ADAMTS8", "BMP6", "CCN3", "CD5", "CD58", "CD40", "PCSK9")
tarfiles <- sapply(p, \(x) { grep(paste0(x, "_"), tarfiles, value=TRUE) })

inst <- readRDS(here("data", "ukb_pqtls.rds"))
p <- c("ACOX1", "ADAMTS8", "BMP6", "CCN3", "CD5", "CD58", "CD40", "PCSK9")
subset(inst, `Assay Target` %in% p) %>% 
    group_by(`Assay Target`, `cis/trans`) %>%
    summarise(n())
variants <- unique(inst$`Variant ID (CHROM:GENPOS (hg37):A0:A1:imp:v1)`)

o <- mclapply(tarfiles, \(x) extract_tar(variants, x))
