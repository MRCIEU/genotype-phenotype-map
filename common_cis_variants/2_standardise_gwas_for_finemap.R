library(argparser, quietly = TRUE)

parser <- arg_parser("Standardising GWAS for finemap")
parser <- add_argument(parser,
                       "--gwas_filename",
                       help = "GWAS filename",
                       type = "character"
)
parser <- add_argument(parser,
                       "--snp_list",
                       help = "List of RSIDs to filter on",
                       type = "character"
)

args <- parse_args(parser)
snps <- vroom::vroom(args$snp_list, col_names = F, delim=",")$X1

gwas <- vroom::vroom(args$gwas_filename, delim = " ") |>
  dplyr::filter(rsid %in% snps)
gwas <- gwas[!duplicated(gwas$rsid), ]

#to_flip <- (gwas$maf > 0.5)
#if (any(to_flip)) {
#  gwas$maf[to_flip] <- 1 - gwas$maf[to_flip]
#  gwas$beta[to_flip] <- -1 * gwas$beta[to_flip]
#
#  temp <- gwas$allele2[to_flip]
#  gwas$allele2[to_flip] <- gwas$allele1[to_flip]
#  gwas$allele1[to_flip] <- temp
#}

vroom::vroom_write(gwas, args$gwas_filename, delim=" ")

