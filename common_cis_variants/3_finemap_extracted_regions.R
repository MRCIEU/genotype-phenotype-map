load("R/gwas_formatting.R")

parser <- argparser::arg_parser("Finemapping stuff")
parser <- argparser::add_argument(parser, "--gwas_filename", help = "GWAS filename", type = "character")
parser <- argparser::add_argument(parser, "--ld_block", help = "LD block that the ", type = "character")
parser <- argparser::add_argument(parser, "--gwas_n", help = "Sample size of GWAS", type = "numeric")

args <- parse_args(parser)

args <- list(gwas_filename="/Users/wt23152/Documents/Projects/scratch/011/data/study/ukb-b-10003/EUR_6_32530029.z",
            snp_list="/Users/wt23152/Documents/Projects/scratch/011/data/study/ukb-b-10003/plink/EUR_6_32530029.snplist",
            ld_region="/Users/wt23152/Documents/Projects/scratch/011/data/ld_block_matrices/EUR/10_93335047_95396367.tsv.gz",
            gwas_n=400000
)

args <- list(ld_block_dir="/Users/wt23152/Documents/Projects/scratch/011/data/ld_blocks/EUR/10/93335047_95396367/",
             ld_region="/Users/wt23152/Documents/Projects/scratch/011/data/ld_block_matrices/EUR/10_93335047_95396367.tsv.gz",
             gwas_n=400000
)


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
  p <- 2 * pnorm(z, lower.tail = F)
  return(data.frame(beta, se, p))
}

susie <- function(study, gwas, ld_matrix) {
  tryCatch(expr = {
    start_time <- Sys.time()
    susie_results <- susieR::susie_rss(bhat=gwas$beta, shat=gwas$se, R=ld_matrix, n=args$gwas_n, N=5)
    print(Sys.time() - start_time)

    saveRDS(susie_results, paste0(data_dir, 'finemap_tests/susie_', basename(study), '.rds'))

    conditioned_gwases <- apply(susie_results$lbf_variable, 1, function(lbf) {
      conditioned_gwas <- lbf_to_z_cont(lbf, args$gwas_n, gwas$maf)
      gwas$beta <- conditioned_gwas$beta
      gwas$se <- conditioned_gwas$se
      gwas$p <- conditioned_gwas$p
      return(gwas)
    })
    return(conditioned_gwases)

  }, warning = function(w) {
    warning(w)
    return(gwas)
  })
}

carma <- function(gwas, ld_matrix) {
  input_columns <- list(BETA="beta", SE="se", RSID="rsid", CHR="chromosome", BP="position", EA="allele1", OA="allele2", EAF="maf", LOG_P="lp")
  gwas <- standarise_gwas(gwas, input_columns=input_columns, N=args$gwas_n)
  gwas <- convert_beta_and_se_to_z_score(gwas)
  lambda_list<-list()
  lambda_list[[1]] <- 1
  #outlier_switch=T because of out sample LD matrix
  carma_result <- CARMA::CARMA(gwas$Z, ld_matrix, lamda_list, outlier_switch=TRUE)
  saveRDS(carma_result, paste0(data_dir, '/finemap_tests/carma_', basename(study), '.rds'))
}


ld_region <- vroom::vroom(args$ld_region, col_names=F)
colnames(ld_region) <- ld_region$X1
ld_region <- ld_region[, 2:(ncol(ld_region)-1)]

all_studies <- list.files(args$ld_block_dir, pattern = "*.z", full.names = T)
all_studies <- head(all_studies)

all_conditioned_gwases <- lapply(all_studies, function(study) {
  print(study)
  gwas <- vroom::vroom(study, delim = " ")
  if (typeof(gwas$maf) == 'character') {
    warning('maf is poo, skipping')
    return(gwas)
  }

  gwas <- dplyr::filter(gwas, rsid %in% colnames(ld_region))
  keep <- colnames(ld_region) %in% gwas$rsid
  ld_for_gwas <- ld_region[keep, keep]
  ld_matrix <- matrix(as.vector(data.matrix(ld_for_gwas)), nrow=nrow(gwas), ncol=nrow(gwas))
  conditioned_gwases <- susie(study, gwas, ld_matrix)
})

flatten_list_of_lists <- function(x) {
  if (is.data.frame(x)) return(list(x))
  if (!is.list(x)) return(x)
  unlist(lapply(x, flattenMixed), FALSE)
}

aall_conditioned_gwases <- flatten_list_of_lists(all_conditioned_gwases)

lapply(aall_conditioned_gwases, function(gwas) {
  range <- c(min(gwas$position), max(gwas$position))
  #yrange <- c(min(gwas$p), max(gwas$p))
  if(!"p" %in% colnames(gwas)) return()
  gwas$p[gwas$p == 0] <- .Machine$double.xmin
  gwas <- gwas[!is.na(gwas$p), ]
  qqman::manhattan(gwas, p='p', chr = 'chromosome', bp = 'position', snp = 'rsid', xlim=range)
})
