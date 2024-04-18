library(TwoSampleMR)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(here)

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

load(here("data", "extracted.rdata"))

pqtlf <- pqtl %>%
tidyr::separate(`Variant ID (CHROM:GENPOS (hg37):A0:A1:imp:v1)`, into=c("chr", "pos", "a0", "a1", "imp", "v1"), sep=":", remove=FALSE) %>%
format_data(.,
    snp_col="rsID",
    beta_col="BETA (discovery, wrt. A1)",
    se_col="SE (discovery)",
    eaf_col="A1FREQ (discovery)",
    effect_allele_col="a0",
    other_allele_col="a1",
    phenotype_col="Assay Target"
) %>% left_join(
    ., subset(pqtl, select=c(rsID, `cis/trans`, `Assay Target`)),
    by=c("exposure"="Assay Target", "SNP"="rsID")
)
dim(pqtlf)


ieugwasr_to_2smr <- function(d) {
    format_data(d, 
        type="outcome",
        snp_col="rsid",
        chr_col="chr", pos_col="position", beta_col="beta", se_col="se", pval_col="p", phenotype_col="id", effect_allele_col="ea", other_allele_col="nea", eaf_col="eaf"
    )
}

dat_dcis <- lapply(1:nrow(ap), \(i) {
    harmonise_data(
        subset(pqtlf, exposure == ap$protein[i] & `cis/trans` == "cis"),
        out_cistrans[[i]] %>% ieugwasr_to_2smr,
        action=1
    )
}) %>% bind_rows()

dat_dtrans <- lapply(1:nrow(ap), \(i) {
    harmonise_data(
        subset(pqtlf, exposure == ap$protein[i] & `cis/trans` == "trans"),
        out_cistrans[[i]] %>% ieugwasr_to_2smr,
        action=1
    )
}) %>% bind_rows()

dat_itrans <- lapply(1:nrow(ap), \(i) {
    a <- format_data(exp_trans[[i]] %>% mutate(exp=names(exp_trans)[i]), type="exposure", snp_col="rsID", beta_col="BETA", se_col="SE", effect_allele_col="a1", other_allele_col="a0", eaf_col="A1FREQ", pval_col="p", phenotype_col="exp")
    b <- out_trans[[i]] %>% ieugwasr_to_2smr
    harmonise_data(a, b, action=1)
}) %>% bind_rows()


res <- bind_rows(
    mr(dat_dcis) %>% mutate(instrument="cis"), 
    mr(dat_dtrans, method_list="mr_ivw") %>% mutate(instrument="discovery trans"), 
    mr(dat_itrans, method_list="mr_ivw") %>% mutate(instrument = "indirect cis")
)

ap$code <- paste(ap$protein, "->", ap$trait)

inner_join(res, ap, by=c("exposure"="protein")) %>%
ggplot(., aes(x=b, y=instrument)) +
    geom_vline(xintercept=0, linetype="dotted") +
    geom_point(size=3, aes(colour=instrument)) +
    geom_label_repel(size=2, aes(label=nsnp)) +
    geom_errorbarh(aes(colour=instrument, xmin=b-se*1.96, xmax=b+se*1.96), height=0) +
    facet_wrap(~ code, scale="free_x") +
    scale_colour_brewer(type="qual") +
    theme(legend.position="none")
ggsave(here("results", "initial_cistrans.pdf"))








