weak_signal_example <- list(
    ld_matrix_file = '/local-scratch/projects/genotype-phenotype-map/data/ld_reference_panel_hg38/EUR/12/103736758-105970320.unphased.vcor1',
    ld_matrix_info_file = '/local-scratch/projects/genotype-phenotype-map/data/ld_reference_panel_hg38/EUR/12/103736758-105970320.tsv',
    gwas_file = '/local-scratch/projects/genotype-phenotype-map/test/data/study/ebi-a-GCST90002304/standardised/EUR_12_104476236.tsv.gz'
)

strong_signal_example <- list(
    ld_matrix_file = '/local-scratch/projects/genotype-phenotype-map/data/ld_reference_panel_hg38/EUR/2/45817500-46829710.unphased.vcor1',
    ld_matrix_info_file = '/local-scratch/projects/genotype-phenotype-map/data/ld_reference_panel_hg38/EUR/2/45817500-46829710.tsv',
    gwas_file = '/local-scratch/projects/genotype-phenotype-map/test/data/study/ebi-a-GCST90002304/standardised/EUR_2_46126027.tsv.gz'
)

complex_ld_region_example <- list(
    ld_matrix_file = '/local-scratch/projects/genotype-phenotype-map/data/ld_reference_panel_hg38/EUR/6/31282730-32518958.unphased.vcor1',
    ld_matrix_info_file = '/local-scratch/projects/genotype-phenotype-map/data/ld_reference_panel_hg38/EUR/6/31282730-32518958.tsv',
    gwas_file = '/local-scratch/projects/genotype-phenotype-map/test/data/study/ebi-a-GCST90002304/standardised/EUR_6_31344761.tsv.gz'
)

example <- complex_ld_region_example 

ld_matrix <- vroom::vroom(example$ld_matrix_file, col_names = F, show_col_types = F) |> data.matrix()
ld_matrix_info <- vroom::vroom(example$ld_matrix_info_file, show_col_types = F)
gwas <- vroom::vroom(example$gwas_file, show_col_types = F)


#' Basic imputation function
#' 
#' @param R The correlation matrix - must be complete for the set of SNPs that need to be imputed
#' @param ss A data frame with columns betahat2 = vector of effect estimates in the same order as R and with NAs for variants that need to be imputed; se = as with betahat2 but for available standard errors, af = allele frequencies (no missing values allowed, so use reference panel if there are missing values)
#' @param index The positions of the SNPs that are causal and will be used to generate the simulated summary statistics. This can just be the top hit.
#' 
#' @return A list with the following elements:
#' - ss: The input data frame with the imputed values added
#' - b_adj: The adjustment factor for the effect sizes
#' - se_adj: The adjustment factor for the standard errors
#' - b_cor: The correlation between the true and imputed effect sizes - this is critical for evaluation of the performance of the imputation, it should be close to 1 e.g > 0.7 would be a reasonable threshold
#' - se_cor: The correlation between the true and imputed standard errors
outlier_detection <- function(r) {
    sd1 <- sd(r, na.rm=T)
    m <- median(r, na.rm=T)
    r[r < (m - 3 * sd1)] <- NA
    r[r > (m + 3 * sd1)] <- NA
    sd2 <- sd(r, na.rm=T)
    r[r < (m - 3 * sd2)] <- NA
    r[r > (m + 3 * sd2)] <- NA
    outliers <- is.na(r)
    return(outliers)
}


adjust <- function(truth, predicted) {
    outs <- outlier_detection(truth / predicted)
    reg <- lm(truth[!outs] ~ predicted[!outs])
    adj <- predicted * reg$coef[2] - reg$coef[1]
    corr <- cor(adj[!outs], truth[!outs], use="pair")
    return(list(adj=adj, outliers = outs, adj_slope=reg$coef[2], adj_intercept=reg$coef[1], corr=corr))
}


#' Basic imputation function
#' 
#' @param gwas A data frame with columns betahat2 = vector of effect estimates in the same order as ld_matrix and with NAs for variants that need to be imputed;
#'  se = as with betahat2 but for available standard errors, af = allele frequencies (no missing values allowed, so use reference panel if there are missing values)
#' @param ld_matrix The correlation matrix - must be complete for the set of SNPs that need to be imputed
#' @param index The positions of the SNPs that are causal and will be used to generate the simulated summary statistics. This can just be the top hit.
#' 
#' @return A list with the following elements:
#' - gwas: The input data frame with the imputed values added
#' - b_adj: The adjustment factor for the effect sizes
#' - se_adj: The adjustment factor for the standard errors
#' - b_cor: The correlation between the true and imputed effect sizes - this is critical for evaluation of the performance of the imputation,
#'      it should be close to 1 e.g > 0.7 would be a reasonable threshold
#' - se_cor: The correlation between the true and imputed standard errors
imp <- function(gwas, ld_matrix, index) {
    b <- gwas$BETA
    se <- gwas$SE
    af <- gwas$EAF
    to_impute <- is.na(b)
    num_to_impute <- sum(to_impute)

    if (num_to_impute == 0) {
      return(list(
          gwas = gwas, b_cor = NA, se_cor = NA, b_adj = NA, se_adj = NA, indices = NA, rows_imputed = 0
      ))
    }

    message("Checking data")
    nsnp <- length(b)
    stopifnot(ncol(ld_matrix) == nsnp)
    stopifnot(nrow(ld_matrix) == nsnp)
    stopifnot(length(af) == nsnp)
    stopifnot(length(se) == nsnp)
    stopifnot(all(index) %in% 1:nsnp)
    stopifnot(length(index) < nsnp)
    stopifnot(all(af > 0 & af < 1))
    stopifnot(all(!is.na(af)))
    stopifnot(all(se > 0, na.rm=TRUE))

    # Calculate the diagonal matrix of variances and the inverse
    message("Generating D")
    D <- diag(sqrt(2 * af * (1 - af)))
    Di <- diag(1 / diag(D))

    # Get the conditional estimates of the index SNP effects
    message("Conditional index variant estimates")
    if(length(index) == 1) {
        bhat2 <- b[index]
    } else {
        bhat2 <- D[index,index] %*% MASS::ginv(ld_matrix[index,index]) %*% Di[index,index] %*% b[index]
    }
    b2 <- rep(0, nsnp)
    b2[index] <- bhat2

    # Get the simulated effect sizes
    message("Simulated effects")
    betahat_sim <- as.numeric(Di %*% ld_matrix %*% D %*% b2)

    # Initialise the SE - this doesn't account for var(y) or sample size, but those are constants that can be obtained from regression re-scaling
    sehat <- (diag(Di))

    # Re-scale effect sizes and standard errors
    message("Rescaling")
    b_adj <- adjust(b, betahat_sim)
    se_adj <- adjust(se, sehat)

    # str(b_adj)
    # str(se_adj)

    # png("temp4.png")
    # plot(b_adj$adj[!b_adj$outliers], b[!b_adj$outliers])
    # abline(0, 1)
    # dev.off()

    # png("temp5.png")
    # plot(se_adj$adj[!se_adj$outliers], se[!se_adj$outliers])
    # abline(0, 1)
    # dev.off()

    gwas$BETA_IMPUTED <- b_adj$adj
    gwas$SE_IMPUTED <- se_adj$adj
    stopifnot(all(!is.na(gwas$BETA_IMPUTED)))
    stopifnot(all(!is.na(gwas$SE_IMPUTED)))

    message("Finalising")
    gwas$BETA[to_impute] <- gwas$BETA_IMPUTED[to_impute]
    gwas$SE[to_impute] <- gwas$SE_IMPUTED[to_impute]
    gwas$Z[to_impute] <- gwas$BETA_IMPUTED[to_impute] / gwas$SE_IMPUTED[to_impute]
    gwas$P[to_impute] <- 2 * pnorm(-abs(gwas$Z[to_impute]))
    gwas$imp <- to_impute


    # png("temp4.png")
    # plot(BETA ~ BETA_IMPUTED, gwas)
    # abline(0,1)
    # dev.off()

    return(
      list(
        gwas = gwas,
        b_adj = b_adj$adj_slope,
        b_adj_intercept = b_adj$adj_intercept,
        se_adj = se_adj$adj_slope,
        se_adj_intercept = se_adj$adj_intercept,
        b_cor = b_adj$corr,
        se_cor = se_adj$corr,
        indices = length(index),
        rows_imputed = num_to_impute,
        rows_region = nrow(gwas)
      )
    )
}

clump_ld_region <- function(z, R, zthresh = qnorm(1.5e-4, low=F), rthresh = 0.01) {
  z <- abs(z)
  z[z < zthresh] <- NA
  k <- c()

  while(!all(is.na(z))) {
    i <- which.max(z)
    k <- c(k, i)
    z[i] <- NA
    z[which(R[i,]^2 > rthresh)] <- NA
  }
  return(k)
}


gwas_to_impute <- dplyr::left_join(
    dplyr::select(ld_matrix_info, -EAF),
    dplyr::select(gwas, -CHR, -BP, -EA, -OA),
    by=dplyr::join_by(SNP)
)

rows_to_impute <- !ld_matrix_info$SNP %in% gwas$SNP
print(glue::glue('rows to impute: {sum(rows_to_impute)}'))
print(glue::glue('rows in gwas: {nrow(gwas)}'))
print(glue::glue('rows in ld matrix: {nrow(ld_matrix_info)}'))
print(glue::glue('gwas order matches ld matrix: {all(ld_matrix_info$SNP == gwas_to_impute$SNP)}'))

gwas_to_impute$EAF[rows_to_impute] <- ld_matrix_info$EAF[rows_to_impute]

clumped_snp_index <- clump_ld_region(gwas_to_impute$Z, ld_matrix)
print('clumped snps:')
print(gwas_to_impute[clumped_snp_index, ])
result <- imp(gwas_to_impute, ld_matrix, clumped_snp_index)

print(result$ss[, c(2,6,7,8,11,12)])
print(result$ss[is.na(gwas_to_impute$BETA), c(2,6,7,8,11,12)])
print(result)



#Gib's example:
# load(url("https://github.com/explodecomputer/lab-book/raw/refs/heads/main/posts/2024-09-18-conditional-summary-stats/1kg_region.rdata"))
