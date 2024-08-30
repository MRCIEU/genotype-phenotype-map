library(glue)
library(dplyr)
library(data.table)

bfiles <- "/local-scratch/projects/genotype-phenotype-map/data/ld_reference_panel/EUR/full_maf"
eur <- fread("pipeline_steps/data/ld_regions.tsv") %>% as_tibble() %>% filter(ancestry == "EUR")
outdir <- "/local-scratch/projects/genotype-phenotype-map/data/ldmat_gib/EUR"
file.exists(outdir)
tfile <- tempfile()
plink <- "/local-scratch/projects/genotype-phenotype-map/bin/plink2"

# Generate plink subsets

for(i in 1:nrow(eur))
{
    message(i)
    chr <- eur$chr[i]
    pos1 <- eur$start[i]
    pos2 <- eur$stop[i]
    out <- glue("{outdir}/{chr}/{pos1}_{pos2}")

    s <- glue("{chr} {pos1} {pos2} a")
    write.table(s, file=tfile, row=F, col=F, qu=F)

    #bfile <- file.path(bfiles, paste0("data.chr", sprintf("%02d", chr)))

    # create subset
    glue("{plink} --bfile {bfiles} --chr {chr} --extract range {tfile} --make-bed --out {out} --keep-allele-order --freq") %>% system()

    # update alleles

    #bim <- data.table::fread(paste0(out, ".bim"))
    #bim$switch <- bim$V5 > bim$V6
    #temp <- bim$V5[bim$switch]
    #bim$V5[bim$switch] <- bim$V6[bim$switch]
    #bim$V6[bim$switch] <- temp
    # table(bim$switch)

    #switchfile <- paste0(tfile, ".switch")
    #data.table::fwrite(subset(bim, select=c(V2, V6)), file=switchfile, quote=FALSE, col.names=FALSE, sep=" ")

    #glue(
    #    "{plink} --bfile {out} --rm-dup force-first --ref-allele {switchfile} 2 1 --make-bed --out {out} --keep-allele-order"
    #) %>% system()

    if(nchar(out) > 0) unlink(glue("{out}*~"))
}

# Generate LD mats and frequencies

for(i in 1:nrow(eur))
{
    message(i)
    chr <- eur$chr[i]
    pos1 <- eur$start[i]
    pos2 <- eur$stop[i]
    out <- glue("{outdir}/{chr}/{pos1}_{pos2}")

    glue("{plink} --bfile {out} --r-unphased square --out {out} --keep-allele-order") %>% system()
    # glue("gzip {out}.unphased.vcor1") %>% system()
    glue("{plink} --bfile {out} --out {out} --keep-allele-order --freq") %>% system()
}


# Create a single individual level dataset of all the regions
#mergelist <- paste0(outdir, "/", eur$chr, "/", eur$start, "_", eur$stop)
#mergefile <- tempfile()
#write.table(mergelist, file=mergefile, row=F, col=F, qu=F)
#out <- "/local-scratch/projects/genotype-phenotype-map/data/ldmat_gib/EUR"
#glue("plink1.9 --merge-list {mergefile} --make-bed --out {out} --keep-allele-order") %>% system()
