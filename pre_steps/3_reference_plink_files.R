library(glue)
library(dplyr)
library(data.table)
library(here)

bfiles <- "/local-scratch/projects/Lifecourse-GWAS/ukb/geno_input/"
eur <- fread("eur_ldregions.txt") %>% as_tibble()
outdir <- "/local-scratch/projects/genotype-phenotype-map/data/ldmat_gib/eur"
file.exists(outdir)
tfile <- tempfile()
keepfile <- here("ldmat", "ukb50k.ids")
plink <- "/local-scratch/projects/genotype-phenotype-map/bin/plink2"

for(i in 1:nrow(eur))
{
    message(i)
    chr <- eur$V1[i]
    pos1 <- eur$V2[i]
    pos2 <- eur$V3[i]
    out <- glue("{outdir}/{chr}/{pos1}-{pos2}")

    s <- glue("{chr} {pos1} {pos2} a")
    write.table(s, file=tfile, row=F, col=F, qu=F)

    bfile <- file.path(bfiles, paste0("data.chr", sprintf("%02d", chr)))

    # create subset
    glue("{plink} --bfile {bfile} --chr {chr} --extract range {tfile} --keep {keepfile} --make-bed --out {out} --keep-allele-order") %>% system()

    # update alleles

    bim <- data.table::fread(paste0(out, ".bim"))
    bim$switch <- bim$V5 > bim$V6
    temp <- bim$V5[bim$switch]
    bim$V5[bim$switch] <- bim$V6[bim$switch]
    bim$V6[bim$switch] <- temp
    # table(bim$switch)

    switchfile <- paste0(tfile, ".switch")
    data.table::fwrite(subset(bim, select=c(V2, V6)), file=switchfile, quote=FALSE, col.names=FALSE, sep=" ")

    glue(
        "{plink} --bfile {out} --rm-dup force-first --ref-allele {switchfile} 2 1 --make-bed --out {out} --keep-allele-order"
    ) %>% system()

    if(nchar(out) > 0) unlink(glue("{out}*~"))
}






