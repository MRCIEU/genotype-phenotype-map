ld_matrix_file <- '/local-scratch/projects/genotype-phenotype-map/data/ld_reference_panel_hg38/EUR/12/103736758-105970320.unphased.vcor1'
ld_matrix_info_file <- '/local-scratch/projects/genotype-phenotype-map/data/ld_reference_panel_hg38/EUR/12/103736758-105970320.tsv'
gwas_file <- '/local-scratch/projects/genotype-phenotype-map/test/data/study/ebi-a-GCST90002304/standardised/EUR_12_104476236.tsv.gz'

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
imp <- function(ss, R, index) {
    b <- ss$BETA
    se <- ss$SE
    af <- ss$EAF
    nsnp <- length(b)
    stopifnot(ncol(R) == nsnp)
    stopifnot(nrow(R) == nsnp)
    stopifnot(length(af) == nsnp)
    stopifnot(length(se) == nsnp)
    stopifnot(all(index) %in% 1:nsnp)
    stopifnot(length(index) < nsnp)
    stopifnot(all(af > 0 & af < 1))
    stopifnot(all(!is.na(af)))
    stopifnot(all(se > 0, na.rm=TRUE))
    if(all(!is.na(b))) {
        message("No missing values in b, imputation not required")
        b_cor=1
        se_cor=1
        mod1=1
        mod2=1
    } else {
        # Calculate the diagonal matrix of variances and the inverse
        D <- diag(sqrt(2 * af * (1 - af)))
        Di <- diag(1 / diag(D))

        # Get the conditional estimates of the index SNP effects
        if(length(index) == 1) {
            bhat2 <- b[index]
        } else {
            bhat2 <- D[index,index] %*% MASS::ginv(R[index,index]) %*% Di[index,index] %*% b[index]
        }
        b2 <- rep(0, nsnp)
        b2[index] <- bhat2

        # Get the simulated effect sizes
        betahat_sim <- as.numeric(Di %*% R %*% D %*% b2)

        # Initialise the SE - this doesn't account for var(y) or sample size, but those are constants that can be obtained from regression re-scaling
        sehat <- sqrt(diag(Di))

        # Re-scale effect sizes and standard errors
        # vb <- var(b, na.rm=TRUE)
        # vse <- var(se, na.rm=TRUE)
        # mod1 <- cov(b, betahat_sim, use="pair") / vb
        mod1 <- lm(betahat_sim ~ b)$coef[2]
        # mod2 <- cov(se, sehat, use="pair") / vse
        mod2 <- lm(sehat ~ se)$coef[2]

        # Performance metrics
        # b_cor = mod1 * sqrt(vb) / sd(betahat_sim, na.rm=TRUE)
        b_cor <- cor(b, betahat_sim, use="pair")
        # se_cor = mod2 * sqrt(vse) / sd(sehat, na.rm=TRUE)
        se_cor <- cor(se, sehat, use="pair")

        # Re-scale
        betahat_sim <- betahat_sim / mod1
        sehat <- sehat / mod2

        # Fill in missing values
        b[is.na(b)] <- betahat_sim[is.na(b)]
        se[is.na(se)] <- sehat[is.na(se)]

        stopifnot(all(!is.na(b)))
        stopifnot(all(!is.na(se)))
    }

    ss$betahatimp <- b
    ss$seimp <- se
    ss$zimp <- b / se
    ss$pimp <- 2 * pnorm(-abs(ss$zimp))

    # Output
    out <- list(
        ss = ss,
        b_adj = mod1,
        se_adj = mod2,
        b_cor = b_cor,
        se_cor = se_cor,
        n_ind = length(index)
    )
    return(out)
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

ld_matrix <- vroom::vroom(ld_matrix_file, col_names = F, show_col_types = F) |> data.matrix()
ld_region_from_reference_panel <- vroom::vroom(ld_matrix_info_file, show_col_types = F)
gwas <- vroom::vroom(gwas_file, show_col_types = F)

gwas_to_impute <- dplyr::left_join(
    dplyr::select(ld_region_from_reference_panel, -EAF),
    dplyr::select(gwas, -CHR, -BP, -EA, -OA),
    by=dplyr::join_by(SNP)
)

rows_to_impute <- !ld_region_from_reference_panel$SNP %in% gwas$SNP
print(glue::glue('rows to impute: {sum(rows_to_impute)}'))
print(glue::glue('rows in gwas: {nrow(gwas)}'))
print(glue::glue('rows in ld matrix: {nrow(ld_region_from_reference_panel)}'))
print(glue::glue('gwas order matches ld matrix: {all(ld_region_from_reference_panel$SNP == gwas_to_impute$SNP)}'))

gwas_to_impute$EAF[rows_to_impute] <- ld_region_from_reference_panel$EAF[rows_to_impute]


clumped_snp_index <- clump_ld_region(gwas_to_impute$Z, ld_matrix)
print('clumped snps:')
print(gwas_to_impute[clumped_snp_index, ])
result <- imp(gwas_to_impute, ld_matrix, clumped_snp_index)

print(result$ss[, c(2,6,7,8,11,12)])
print(result$ss[is.na(gwas_to_impute$BETA), c(2,6,7,8,11,12)])
print(result)

