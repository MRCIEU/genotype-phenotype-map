library(glue)
library(dplyr)
library(data.table)
library(here)

bfiles <- "/local-scratch/projects/Lifecourse-GWAS/ukb/geno_input/"
eur <- fread(here("pre_steps", "eur_ldregions.txt")) %>% as_tibble()
outdir <- "/local-scratch/projects/genotype-phenotype-map/data/ldmat_gib/EUR"
file.exists(outdir)
tfile <- withr::local_tempfile()
plink <- "plink1.9"

# Create a single individual level dataset of all the regions
mergelist <- paste0(outdir, "/", eur$V1, "/", eur$V2, "-", eur$V3)
mergefile <- withr::local_tempfile()
write.table(mergelist, file = mergefile, row = F, col = F, qu = F)
out <- "/local-scratch/projects/genotype-phenotype-map/data/ldmat_gib/EUR/ukb50k"
glue("{plink} --merge-list {mergefile} --make-bed --out {out} --keep-allele-order") %>% system()
