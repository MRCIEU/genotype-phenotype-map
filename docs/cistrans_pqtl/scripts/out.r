library(tidyr)
library(ieugwasr)
library(dplyr)
library(here)

select_api("private")

ap <- tribble(
~protein, ~id, ~trait,
"ACOX1", "ieu-b-38", "SBP",
"ADAMTS8", "ieu-b-38", "SBP",
"BMP6", "ieu-b-38", "SBP",
"CCN3", "ieu-b-38", "SBP",
"CD5", "ieu-b-18", "MS",
"CD58", "ieu-b-18", "MS",
"CD40", "ieu-b-18", "MS",
"PCSK9", "ieu-a-300", "LDL",
)


pqtl <- readRDS(here("data", "ukb_pqtls.rds"))

subset(pqtl, `Assay Target` %in% p) %>% 
    group_by(`Assay Target`, `cis/trans`) %>%
    summarise(n())


expdat <- readRDS(here("data", "pqtl_protein_extract.rds"))

expdat


# cistrans



out_cistrans <- lapply(1:nrow(ap), \(i) {
    message(i)
    associations(subset(pqtl, `Assay Target` == ap$protein[i])$rsID, ap$id[i])
})

out_cistrans

mapping <- subset(pqtl, select=c(`Variant ID (CHROM:GENPOS (hg37):A0:A1:imp:v1)`, `rsID`)) %>%
filter(!duplicated(rsID))


exp_trans <- lapply(expdat, \(x) {
    x %>% 
        filter(!duplicated(ID)) %>%
        mutate(p=10^-LOG10P, fdr=p.adjust(p, "fdr")) %>% 
        filter(fdr < 0.05) %>%
        tidyr::separate(ID, into=c("chr", "pos", "a0", "a1", "imp", "v1"), sep=":", remove=FALSE) %>%
        left_join(., mapping, c("ID"="Variant ID (CHROM:GENPOS (hg37):A0:A1:imp:v1)"))
})

exp_trans %>% lapply(dim)
out_trans <- lapply(1:nrow(ap), \(i) {
    message(i)
    a <- exp_trans[[i]]$rsID
    a <- a[!is.na(a)]
    associations(a, ap$id[i])
})

temp <- expdat[[1]] %>% 
    filter(!duplicated(ID)) %>%
    mutate(p=10^-LOG10P, fdr=p.adjust(p, "fdr")) %>% 
    filter(fdr < 0.05) %>%
    tidyr::separate(ID, into=c("chr", "pos", "a0", "a1", "imp", "v1"), sep=":", remove=FALSE) %>%
    left_join(., mapping, c("ID"="Variant ID (CHROM:GENPOS (hg37):A0:A1:imp:v1)"))
temp
temp
mapping

save(out_cistrans, exp_trans, out_trans, pqtl, file=here("data", "extracted.rdata"))

