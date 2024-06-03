load("R/gwas_formatting.R")
library(argparser, quietly = TRUE)

parser <- arg_parser("Standardising GWAS for finemap")
parser <- add_argument(parser,
                       "--gwas_filename",
                       help = "GWAS filename",
                       type = "character"
)
parser <- add_argument(parser,
                       "--ld_matrix",
                       help = "LD matrix (created by plink)",
                       type = "character"
)
parser <- add_argument(parser,
                       "--snp_list",
                       help = "List of RSIDs to filter on",
                       type = "character"
)

parser <- add_argument(parser,
                       "--gwas_n",
                       help = "Sample size of GWAS",
                       type = "numeric"
)


args <- parse_args(parser)

#args <- list(gwas_filename="/Users/wt23152/Documents/Projects/scratch/011/data/study/ukb-b-10003/EUR_6_32530029.z",
#             snp_list="/Users/wt23152/Documents/Projects/scratch/011/data/study/ukb-b-10003/plink/EUR_6_32530029.snplist",
#             ld_matrix="/Users/wt23152/Documents/Projects/scratch/011/data/study/ukb-b-10003/plink/EUR_6_32530029.ld",
#             gwas_n=400000
#)

snps <- vroom::vroom(args$snp_list, col_names = F, delim=",")$X1

gwas <- vroom::vroom(args$gwas_filename, delim = " ") |>
  dplyr::filter(rsid %in% snps)
gwas <- gwas[!duplicated(gwas$rsid), ]

ld <- vroom::vroom(args$ld_matrix, col_names=F)
ld_matrix <- matrix(as.vector(data.matrix(ld)), nrow=nrow(gwas), ncol=nrow(gwas))

# to_flip <- (gwas$maf > 0.5)
# if (any(to_flip)) {
#   gwas$maf[to_flip] <- 1 - gwas$maf[to_flip]
#   gwas$beta[to_flip] <- -1 * gwas$beta[to_flip]
#
#   temp <- gwas$allele2[to_flip]
#   gwas$allele2[to_flip] <- gwas$allele1[to_flip]
#   gwas$allele1[to_flip] <- temp
# }
# vroom::vroom_write(gwas, args$gwas_filename, delim=" ")

convert_beta_and_se_to_z_score <- function(gwas) {
  gwas$Z <- gwas$BETA / gwas$SE
  return(gwas)
}

susie <- function() {
  #L = 3 takes 1.4 minutes, whereas the default L = 10 takes 4.2 minutes.
  #presumably, L = 3 is sufficient for the small regions were taking
  start_time <- Sys.time()
  susie_results <- susieR::susie_rss(bhat=gwas$beta, shat=gwas$se, R=ld_matrix, n=args$gwas_n, L=3)
  print(Sys.time() - start_time)

  start_time <- Sys.time()
  susie_results_10 <- susieR::susie_rss(bhat=gwas$beta, shat=gwas$se, R=ld_matrix, n=args$gwas_n)
  print(Sys.time() - start_time)
}

carma <- function() {
  input_columns <- list(BETA="beta", SE="se", RSID="rsid", CHR="chromosome", BP="position", EA="allele1", OA="allele2", EAF="maf", LOG_P="lp")
  gwas <- standarise_gwas(gwas, input_columns=input_columns, N=args$gwas_n)
  gwas <- convert_beta_and_se_to_z_score(gwas)
  lambda_list<-list()
  lambda_list[[1]] <- 1
  #add annotation list data?
  CARMA::CARMA(gwas$Z, ld_matrix, lamda_list, outlier_switch=TRUE) #outlier_switch=T because of out sample LD matrix
}
