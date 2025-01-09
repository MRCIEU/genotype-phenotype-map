weak_signal_example <- list(
    ld_matrix_file = '/local-scratch/projects/genotype-phenotype-map/data/ld_reference_panel_hg38/EUR/12/103736758-105970320.unphased.vcor1',
    ld_matrix_eig = '/local-scratch/projects/genotype-phenotype-map/data/ld_reference_panel_hg38/EUR/12/103736758-105970320.ldeig.rds',
    ld_matrix_info_file = '/local-scratch/projects/genotype-phenotype-map/data/ld_reference_panel_hg38/EUR/12/103736758-105970320.tsv',
    gwas_file = '/local-scratch/projects/genotype-phenotype-map/test/data/study/ebi-a-GCST90002304/standardised/EUR_12_104476236.tsv.gz'
)

strong_signal_example <- list(
    ld_matrix_file = '/local-scratch/projects/genotype-phenotype-map/data/ld_reference_panel_hg38/EUR/2/45817500-46829710.unphased.vcor1',
    ld_matrix_eig = '/local-scratch/projects/genotype-phenotype-map/data/ld_reference_panel_hg38/EUR/2/45817500-46829710.ldeig.rds',
    ld_matrix_info_file = '/local-scratch/projects/genotype-phenotype-map/data/ld_reference_panel_hg38/EUR/2/45817500-46829710.tsv',
    gwas_file = '/local-scratch/projects/genotype-phenotype-map/test/data/study/ebi-a-GCST90002304/standardised/EUR_2_46126027.tsv.gz'
)

complex_ld_region_example <- list(
    ld_matrix_file = '/local-scratch/projects/genotype-phenotype-map/data/ld_reference_panel_hg38/EUR/6/31282730-32518958.unphased.vcor1',
    ld_matrix_eig = '/local-scratch/projects/genotype-phenotype-map/data/ld_reference_panel_hg38/EUR/6/31282730-32518958.ldeig.rds',
    ld_matrix_info_file = '/local-scratch/projects/genotype-phenotype-map/data/ld_reference_panel_hg38/EUR/6/31282730-32518958.tsv',
    gwas_file = '/local-scratch/projects/genotype-phenotype-map/test/data/study/ebi-a-GCST90002304/standardised/EUR_6_31344761.tsv.gz'
)


example <- complex_ld_region_example 

ld_matrix <- vroom::vroom(example$ld_matrix_file, col_names = F, show_col_types = F) |> data.matrix()
ld_matrix_eig <- readRDS(example$ld_matrix_eig)
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
outlier_detection <- function(r, thresh=3) {
    sd1 <- sd(r, na.rm=T)
    m <- median(r, na.rm=T)
    r[r < (m - thresh * sd1)] <- NA
    r[r > (m + thresh * sd1)] <- NA
    sd2 <- sd(r, na.rm=T)
    r[r < (m - thresh * sd2)] <- NA
    r[r > (m + thresh * sd2)] <- NA
    outliers <- is.na(r)
    return(outliers)
}


adjust <- function(truth, predicted, eval_frac = 0.5) {
    outs <- outlier_detection(truth / predicted)
    reg <- lm(truth[!outs] ~ predicted[!outs])
    adj <- predicted * reg$coef[2] + reg$coef[1]
    corr <- cor(adj[!outs], truth[!outs], use="pair")
    iqr <- truth > quantile(truth, 1-(eval_frac/2), na.rm=T) | truth < quantile(truth, eval_frac/2, na.rm=T)
    corrw <- cor(adj[!outs & iqr], truth[!outs & iqr], use="pair")

    return(list(adj=adj, outliers = outs, adj_slope=reg$coef[2], adj_intercept=reg$coef[1], corr=corr, corrw=corrw))
}

#' Basic imputation function
#' 
#' @param gwas A data frame with columns BETA = vector of effect estimates in the same order as ld_matrix and with NAs for variants that need to be imputed;
#'  SE = as with BETA but for available standard errors, EAF = allele frequencies (no missing values allowed, so use reference panel if there are missing values)
#' @param pc list of values and vectors from PCA of LD matrix
#' @param thresh Fraction of variance of LD matrix to use for projection.
#' @param eval_frac Fraction of largest betas to use to test for agreement of imputation and observed betas
#' 
#' @return A list with the following elements:
#' - gwas: The input data frame with the imputed values added
#' - b_adj: The adjustment factor for the effect sizes
#' - se_adj: The adjustment factor for the standard errors
#' - b_cor: The correlation between the true and imputed effect sizes - this is critical for evaluation of the performance of the imputation,
#'      it should be close to 1 e.g > 0.7 would be a reasonable threshold
#' - se_cor: The correlation between the true and imputed standard errors
perform_imputation_eig <- function(gwas, pc, thresh=0.9, eval_frac=0.25) {
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
    stopifnot(nrow(pc$vectors) == nsnp)
    stopifnot(length(af) == nsnp)
    stopifnot(length(se) == nsnp)
    stopifnot(all(af > 0 & af < 1))
    stopifnot(all(!is.na(af)))
    stopifnot(all(se > 0, na.rm=TRUE))

    # Initialise the SE - this doesn't account for var(y) or sample size, but those are constants that can be obtained from regression re-scaling

    D <- diag(sqrt(2 * af * (1 - af)))
    Di <- diag(1 / diag(D))
    sehat <- (diag(Di))
    se_adj <- adjust(se, sehat)
    gwas$SE_IMPUTED <- se_adj$adj
    stopifnot(all(!is.na(gwas$SE_IMPUTED)))

    # Sometimes SE is very far away from SE_IMPUTED.
    # This could cause problems if the beta is instable but still used for imputation
    # Set those betas to NA to be imputed
    # This is a form of smoothing of the data
    # Those excluded betas should now be a function of neighbouring LD
    gwas$SE[to_impute] <- gwas$SE_IMPUTED[to_impute]
    se_outliers <- set_se_outliers_missing(gwas$SE, gwas$SE_IMPUTED)
    b[se_outliers] <- NA
    to_impute <- is.na(b) | se_outliers
    num_se_outliers <- sum(se_outliers)

    # Readjust SE
    gwas$SE[se_outliers] <- gwas$SE_IMPUTED[se_outliers]
    se_adj2 <- adjust(gwas$SE, gwas$SE_IMPUTED)
    gwas$SE_IMPUTED <- se_adj2$adj
    gwas$SE[to_impute] <- gwas$SE_IMPUTED[to_impute]

    # Perform beta imputation
    imp <- eig_imp(pc, thresh, b)
    betahat_sim <- imp$dat$X


    # Re-scale effect sizes and standard errors
    message("Rescaling")
    b_adj <- adjust(b, betahat_sim)

    gwas$BETA_IMPUTED <- b_adj$adj
    stopifnot(all(!is.na(gwas$BETA_IMPUTED)))
    
    message("Finalising")
    gwas$BETA[to_impute] <- gwas$BETA_IMPUTED[to_impute]
    gwas$SE[to_impute] <- gwas$SE_IMPUTED[to_impute]
    gwas$Z[to_impute] <- gwas$BETA_IMPUTED[to_impute] / gwas$SE_IMPUTED[to_impute]
    gwas$P[to_impute] <- 2 * pnorm(-abs(gwas$Z[to_impute]))
    gwas$imp <- to_impute

    return(
      list(
        gwas = gwas,
        b_adj = b_adj$adj_slope,
        b_adj_intercept = b_adj$adj_intercept,
        se_adj = se_adj$adj_slope,
        se_adj_intercept = se_adj$adj_intercept,
        b_cor = b_adj$corr,
        b_corr_top = b_adj$corrw,
        se_cor = se_adj$corr,
        se_corr_top = se_adj$corrw,
        rows_imputed = num_to_impute,
        num_se_outliers = num_se_outliers,
        rows_region = nrow(gwas),
        n_components = imp$ncomp
      )
    )
}

eig_imp <- function(pc, thresh, X) {
    message("Getting threshold")
    i <- which(cumsum(pc$values) / sum(pc$values) >= thresh)[1]
    E <- pc$vectors[,1:i]
    mask <- is.na(X)
    dat <- dplyr::bind_cols(tibble::tibble(X), tibble::as_tibble(E))
    message("Fitting model")
    mod <- lm(X ~ ., data=dat)
    rsq <- summary(mod)$adj.r.squared
    message("Imputing")
    p <- predict(mod, dat)
    co <- cor(p, X, use="pair")

    # png("temp2.png")
    # plot(p, X)
    # abline(0,1)
    # dev.off()

    temp <- tibble::tibble(X=p, mask, pos=1:length(X))
    return(list(dat=temp, rsq=rsq, cor=co, ncomp=i))
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

set_se_outliers_missing <- function(se, se_imputed, outthresh = 3) {
    mod <- lm(se ~ se_imputed)
    cd <- cooks.distance(mod)
    i <- outlier_detection(cd, outthresh)
    return(i)
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
result <- perform_imputation_eig(gwas_to_impute, ld_matrix_eig)

print(result$ss[, c(2,6,7,8,11,12)])
print(result$ss[is.na(gwas_to_impute$BETA), c(2,6,7,8,11,12)])
print(result)



#Gib's example:
# load(url("https://github.com/explodecomputer/lab-book/raw/refs/heads/main/posts/2024-09-18-conditional-summary-stats/1kg_region.rdata"))
