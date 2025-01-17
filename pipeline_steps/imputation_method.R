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
#' - z_adj: The adjustment factor for the effect sizes
#' - se_adj: The adjustment factor for the standard errors
#' - b_cor: The correlation between the true and imputed effect sizes - this is critical for evaluation of the performance of the imputation,
#'      it should be close to 1 e.g > 0.7 would be a reasonable threshold
#' - se_cor: The correlation between the true and imputed standard errors
perform_imputation <- function(file, gwas, pc, thresh=0.9, eval_frac=0.5) {
  unaltered_gwas <- gwas
  b <- gwas$BETA
  se <- gwas$SE
  z <- b/se
  af <- gwas$EAF
  to_impute <- is.na(b)
  num_to_impute <- sum(to_impute)

  if (num_to_impute == 0) {
    return(list(
      gwas = unaltered_gwas, b_cor = NA, se_cor = NA, z_adj = NA, se_adj = NA, indices = NA, rows_imputed = 0
    ))
  }

  nsnp <- length(b)
  stopifnot(nrow(pc$vectors) == nsnp)

  if (any(af <= 0 | af >= 1)) {
    stop(glue::glue('{file} has funky EAF values :()'))
  }
  stopifnot(all(!is.na(af)))
  if (any(se <= 0, na.rm=TRUE)) {
    stop(glue::glue('{file} has funky SE values :()'))
  }
  stopifnot(all(se > 0, na.rm=TRUE))

  # Initialise the SE - this doesn't account for var(y) or sample size, but those are constants that can be obtained from regression re-scaling
  D <- diag(sqrt(2 * af * (1 - af)))
  Di <- diag(1 / diag(D))
  sehat <- (diag(Di))
  se_adj <- adjust(se, sehat)
  if (any(is.na(se_adj$adj))) {
    warning(glue::glue('{file} could not be adjusted, skipping imputation'))
    return(list(
      gwas = unaltered_gwas, b_cor = NA, se_cor = NA, z_adj = NA, se_adj = NA, indices = NA, rows_imputed = 0
    ))
  }

  gwas$SE_IMPUTED <- se_adj$adj


  # Sometimes SE is very far away from SE_IMPUTED.
  # This could cause problems if the beta is instable but still used for imputation
  # Set those betas to NA to be imputed
  # This is a form of smoothing of the data
  # Those excluded betas should now be a function of neighbouring LD
  gwas$SE[to_impute] <- gwas$SE_IMPUTED[to_impute]
  se_outliers <- set_se_outliers_missing(gwas$SE, gwas$SE_IMPUTED)
  b[se_outliers] <- NA
  z[se_outliers] <- NA
  to_impute <- is.na(z) | se_outliers
  num_se_outliers <- sum(se_outliers)

  # Readjust SE
  gwas$SE[se_outliers] <- gwas$SE_IMPUTED[se_outliers]
  se_adj2 <- adjust(gwas$SE, gwas$SE_IMPUTED)
  gwas$SE_IMPUTED <- se_adj2$adj
  gwas$SE[to_impute] <- gwas$SE_IMPUTED[to_impute]

  # Perform beta imputation
  imp <- eig_imp(pc, thresh, z)
  z_sim <- imp$dat$X

  # Re-scale effect sizes and standard errors
  z_adj <- adjust(z, z_sim)

  gwas$Z_IMPUTED <- z_adj$adj
  stopifnot(all(!is.na(gwas$Z_IMPUTED)))

  gwas$BETA_IMPUTED <- gwas$Z_IMPUTED * gwas$SE_IMPUTED

  gwas$BETA[to_impute] <- gwas$BETA_IMPUTED[to_impute]
  gwas$SE[to_impute] <- gwas$SE_IMPUTED[to_impute]
  gwas$Z[to_impute] <- gwas$BETA_IMPUTED[to_impute] / gwas$SE_IMPUTED[to_impute]
  gwas$P[to_impute] <- 2 * pnorm(-abs(gwas$Z[to_impute]))
  gwas$IMPUTED <- to_impute

  return(
    list(
      gwas = gwas,
      z_adj = z_adj$adj_slope,
      z_adj_intercept = z_adj$adj_intercept,
      se_adj = se_adj$adj_slope,
      se_adj_intercept = se_adj$adj_intercept,
      b_cor = z_adj$corr,
      b_corr_top = z_adj$corrw,
      se_cor = se_adj$corr,
      se_corr_top = se_adj$corrw,
      rows_imputed = num_to_impute,
      num_se_outliers = num_se_outliers,
      rows_region = nrow(gwas),
      n_components = imp$ncomp
    )
  )
}

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

  reg <- lm(truth[!outs] ~ poly(predicted[!outs], 3, raw=T))
  adj <- predicted * reg$coef[2] + predicted^2 * reg$coef[3] + predicted^3 * reg$coef[4] + reg$coef[1]
  if (is.na(reg$coef[4])) {
    reg <- lm(truth[!outs] ~ poly(predicted[!outs], 2, raw=T))
    adj <- predicted * reg$coef[2] + predicted^2 * reg$coef[3] + reg$coef[1]
  }

  corr <- cor(adj[!outs], truth[!outs], use="pair")
  iqr <- truth > quantile(truth, 1-(eval_frac/2), na.rm=T) | truth < quantile(truth, eval_frac/2, na.rm=T)

  if (length(adj[!outs & iqr]) == 0) {
    corrw <- NA
  } else {
    corrw <- cor(adj[!outs & iqr], truth[!outs & iqr], use="pair")
  }

  return(list(adj=adj, outliers = outs, adj_slope=reg$coef[2], adj_intercept=reg$coef[1], corr=corr, corrw=corrw))
}

eig_imp <- function(pc, thresh, X) {
  i <- which(cumsum(pc$values) / sum(pc$values) >= thresh)[1]
  E <- pc$vectors[,1:i]
  mask <- is.na(X)
  dat <- dplyr::bind_cols(tibble::tibble(X), tibble::as_tibble(E))
  mod <- lm(X ~ ., data=dat)
  rsq <- summary(mod)$adj.r.squared
  p <- predict(mod, dat)
  co <- cor(p, X, use="pair")

  temp <- tibble::tibble(X=p, mask, pos=1:length(X))
  return(list(dat=temp, rsq=rsq, cor=co, ncomp=i))
}

set_se_outliers_missing <- function(se, se_imputed, outthresh = 3) {
  mod <- lm(se ~ se_imputed)
  cd <- cooks.distance(mod)
  i <- outlier_detection(cd, outthresh)
  return(i)
}

#' filter_imputation_results: remove results with over-inflated p-vlaues
#' 
#' @param imputed gwas
#' @param ld correlation matrix
#' 
#' First, filters the imputation results to inside the range of the original gwas
#' Second, finds imputed variants with low pvalues, non-imputed rows that are correlated
#' If imputed variant is more significant the non-imputed rows, drop that variant
#' @return filtered gwas
filter_imputation_results <- function(gwas, ld_matrix, min_bp, max_bp) {
  only_keep_inside_gwas_range <- gwas$BP > min_bp & gwas$BP < max_bp 
  gwas <- gwas[only_keep_inside_gwas_range, ]

  snps_to_remove <- c()
  snps_to_investigate <- which(gwas$P < lowest_p_value_threshold & gwas$IMPUTED == T)
  for (snp_location in snps_to_investigate) {
    ld_correlations <- which(c(ld_matrix[snp_location, ]) > p_value_filter_correlation_threshold)
    ld_correlations <- ld_correlations[ld_correlations != snp_location]
    gwas_correlations <- gwas[(1:nrow(gwas) %in% ld_correlations) & gwas$IMPUTED == F, ]

    snp <- gwas[snp_location, ]

    if (length(gwas_correlations$P) > 0 && min(gwas_correlations$P * 0.1) > snp$P) {
      snps_to_remove <- c(snps_to_remove, snp_location)
    }
  }

  if (length(snps_to_remove) > 0) gwas <- gwas[-snps_to_remove, ]
  return(list(gwas = gwas,
    significant_rows_imputed = length(snps_to_investigate),
    significant_rows_filtered = length(snps_to_remove))
  )
}