load("R/gwas_formatting.R")
library(argparser, quietly = TRUE)

parser <- arg_parser("Finemapping stuff")
parser <- add_argument(parser,
                       "--gwas_filename",
                       help = "GWAS filename",
                       type = "character"
)
parser <- add_argument(parser,
                       "--ld_block",
                       help = "LD block that the ",
                       type = "character"
)

parser <- add_argument(parser,
                       "--gwas_n",
                       help = "Sample size of GWAS",
                       type = "numeric"
)


args <- parse_args(parser)

args <- list(gwas_filename="/Users/wt23152/Documents/Projects/scratch/011/data/study/ukb-b-10003/EUR_6_32530029.z",
            snp_list="/Users/wt23152/Documents/Projects/scratch/011/data/study/ukb-b-10003/plink/EUR_6_32530029.snplist",
            ld_matrix="/Users/wt23152/Documents/Projects/scratch/011/data/study/ukb-b-10003/plink/EUR_6_32530029.ld",
            gwas_n=400000
)

gwas <- vroom::vroom(args$gwas_filename, delim = " ")
ld <- vroom::vroom(args$ld_matrix, col_names=F)
gwas <- gwas[match(gwas$rsid , ld$X0), ]
ld_for_gwas <- ld[match(ld$X0, gwas$rsid), 2:ncol(ld), ]


ld_matrix <- matrix(as.vector(data.matrix(ld), nrow=nrow(gwas), ncol=nrow(gwas)))

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

convert_negative_log_p_to_p <- function(gwas) {
  gwas$P <- 10^-gwas$LOG_P
  return(gwas)
}


#' Convert log Bayes Factor to summary stats
#'
#' @param lbf p-vector of log Bayes Factors for each SNP
#' @param n Overall sample size
#' @param af p-vector of allele frequencies for each SNP
#' @param prior_v Variance of prior distribution. SuSiE uses 50
#'
#' @return tibble with lbf, af, beta, se, z
#' @export
lbf_to_z_cont <- function(lbf, n, af, prior_v = 50) {
  se <- sqrt(1 / (2 * n * af * (1-af)))
  r <- prior_v / (prior_v + se^2)
  z <- sqrt((2 * lbf - log(sqrt(1-r)))/r)
  beta <- z * se
  return(data.frame(beta, se))
}

susie <- function() {
  #L = 3 takes 1.4 minutes, whereas the default L = 10 takes 4.2 minutes.
  #presumably, L = 3 is sufficient for the small regions were taking
  start_time <- Sys.time()
  susie_results <- susieR::susie_rss(bhat=gwas$beta, shat=gwas$se, R=ld_matrix, n=args$gwas_n, L=3)
  print(Sys.time() - start_time)
  #to act on a warning, can do a tryCatch(expr={}, warning=function(w){})

  conditioned_gwases <- apply(susie_results$lbf_variable, 1, function(lbf) {
    conditioned_gwas <- lbf_to_z_cont(lbf, args$gwas_n, gwas$maf)
    return(conditioned_gwas)
  })
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

conditioned_gwases <- apply(susie$l3$lbf_variable, 1, function(lbf) {
  conditioned_gwas <- lbf_to_z_cont(lbf, args$gwas_n, gwas$maf)
  gwas$beta <- conditioned_gwas$beta
  gwas$se <- conditioned_gwas$se
  gwas$p <- 10^-gwas$lp
  return(gwas)
})

lapply(conditioned_gwases, function(gwas) {qqman::manhattan(gwas, p='p', chr = 'chromosome', bp = 'position', snp = 'rsid')})