p_value_filter_correlation_threshold <- 0.6

#' Basic imputation function
#'
#' @param gwas A data frame with columns BETA = vector of effect estimates in the same order as ld_matrix and with NAs for variants that need to be imputed; # nolint: line_length_linter.
#'  SE = as with BETA but for available standard errors, EAF = allele frequencies (no missing values allowed, so use reference panel if there are missing values) # nolint: line_length_linter.
#' @param pc list of values and vectors from PCA of LD matrix
#' @param thresh Fraction of variance of LD matrix to use for projection.
#' @param eval_frac Fraction of largest betas to use to test for agreement of imputation and observed betas
#'
#' @return A list with the following elements:
#' - gwas: The input data frame with the imputed values added
#' - z_adj_coef1: The adjustment factor for the effect sizes
#' - se_adj_coef1: The adjustment factor for the standard errors
#' - b_cor: The correlation between the true and imputed effect sizes - this is critical for evaluation of the performance of the imputation, # nolint: line_length_linter.
#'      it should be close to 1 e.g > 0.7 would be a reasonable threshold
#' - se_cor: The correlation between the true and imputed standard errors
perform_imputation <- function(file, gwas, pc, thresh = 0.9, eval_frac = 0.5) {
  unaltered_gwas <- gwas
  b <- gwas$BETA
  se <- gwas$SE
  z <- b / se
  af <- gwas$EAF
  to_impute <- is.na(b)
  num_to_impute <- sum(to_impute)

  if (num_to_impute == 0) {
    return(list(
      gwas = unaltered_gwas, b_cor = NA, se_cor = NA,
      z_adj_coef1 = NA, se_adj_coef1 = NA,
      indices = NA, rows_imputed = 0
    ))
  }

  nsnp <- length(b)
  stopifnot(nrow(pc$vectors) == nsnp)

  if (any(af <= 0 | af >= 1)) {
    stop(glue::glue("{file} has funky EAF values :()"))
  }
  stopifnot(all(!is.na(af)))
  if (any(se <= 0, na.rm = TRUE)) {
    stop(glue::glue("{file} has funky SE values :()"))
  }
  stopifnot(all(se > 0, na.rm = TRUE))

  # Initialise the SE - this doesn't account for var(y) or sample size, but those are constants that can be obtained from regression re-scaling # nolint: line_length_linter.
  D <- diag(sqrt(2 * af * (1 - af)))
  Di <- diag(1 / diag(D))
  sehat <- (diag(Di))
  se_adj <- adjust(se, sehat, npoly = 1)
  if (any(is.na(se_adj$adj))) {
    warning(glue::glue("{file} could not be adjusted, skipping imputation"))
    return(list(
      gwas = unaltered_gwas, b_cor = NA, se_cor = NA,
      z_adj_coef1 = NA, se_adj_coef1 = NA,
      indices = NA, rows_imputed = 0
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
  se_adj2 <- adjust(gwas$SE, gwas$SE_IMPUTED, npoly = 1)
  gwas$SE_IMPUTED <- se_adj2$adj
  gwas$SE[to_impute] <- gwas$SE_IMPUTED[to_impute]

  # Perform z imputation
  imp <- eig_imp_lasso(pc, thresh, z)
  if (is.null(imp)) {
    return(list(
      gwas = unaltered_gwas, b_cor = NA, se_cor = NA,
      z_adj_coef1 = NA, se_adj_coef1 = NA,
      indices = NA, rows_imputed = 0
    ))
  }
  z_sim <- imp$dat$X

  # Re-scale effect sizes and standard errors
  z_adj <- adjust(z, z_sim, npoly = 3)

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
      z_adj_coef1 = z_adj$adj_coef1,
      z_adj_coef2 = z_adj$adj_coef2,
      z_adj_coef3 = z_adj$adj_coef3,
      se_adj_coef1 = se_adj$adj_coef1,
      se_adj_coef2 = se_adj$adj_coef2,
      se_adj_coef3 = se_adj$adj_coef3,
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


outlier_detection <- function(r, thresh = 3) {
  sd1 <- sd(r, na.rm = T)
  m <- median(r, na.rm = T)
  r[r < (m - thresh * sd1)] <- NA
  r[r > (m + thresh * sd1)] <- NA

  sd2 <- sd(r, na.rm = T)
  r[r < (m - thresh * sd2)] <- NA
  r[r > (m + thresh * sd2)] <- NA
  outliers <- is.na(r)
  return(outliers)
}

rescale_se <- function(truth_se, predicted_se) {
  if (min(predicted_se, na.rm = T) > 0) {
    return(predicted_se)
  }

  truth_min <- min(truth_se, na.rm = T)
  truth_max <- max(truth_se, na.rm = T)
  truth_range <- truth_max - truth_min
  predicted_min <- min(predicted_se, na.rm = T)
  predicted_max <- max(predicted_se, na.rm = T)
  predicted_range <- predicted_max - predicted_min
  se_rescaled <- predicted_se / predicted_range * truth_range + truth_min
  return(se_rescaled)
}

adjust <- function(truth, predicted, eval_frac = 0.5, debug = NULL, npoly = 3) {
  outs <- outlier_detection(truth / predicted)

  truth_clean <- truth[!outs]
  predicted_clean <- predicted[!outs]

  reg <- lm(truth_clean ~ 0 + poly(predicted_clean, npoly, raw = T))
  adj <- rep(0, length(predicted))
  for (i in 1:npoly) {
    adj <- adj + predicted^i * reg$coef[i]
    if (is.na(reg$coef[i])) {
      warning(glue::glue("Coefficient {i} is NA, stopping adjustment at {i-1}"))
      break
    }
  }

  if (!is.null(debug)) {
    # fit <- lm(mpg ~ hp + I(hp^2), data = mtcars)
    prd <- data.frame(
      predicted_clean = seq(from = range(predicted_clean)[1], to = range(predicted_clean)[2], length.out = 100)
    )
    err <- predict(reg, newdata = prd, se.fit = TRUE)
    prd$lci <- err$fit - 1.96 * err$se.fit
    prd$fit <- err$fit
    prd$uci <- err$fit + 1.96 * err$se.fit

    p1 <- ggplot2::ggplot(prd, aes(x = predicted_clean, y = fit)) +
      theme_bw() +
      geom_line() +
      geom_point(data = data.frame(predicted_clean, truth_clean), aes(x = predicted_clean, y = truth_clean)) +
      geom_smooth(aes(ymin = lci, ymax = uci), stat = "identity")
    ggsave(p1, filename = debug, width = 8, height = 6)
  }

  corr <- cor(adj[!outs], truth_clean, use = "pair")
  iqr <- truth > quantile(truth, 1 - (eval_frac / 2), na.rm = T) | truth < quantile(truth, eval_frac / 2, na.rm = T)

  if (length(adj[!outs & iqr]) == 0) {
    corrw <- NA
  } else {
    corrw <- cor(adj[!outs & iqr], truth[!outs & iqr], use = "pair")
  }

  return(list(
    adj = adj, outliers = outs,
    adj_coef1 = reg$coef[1], adj_coef2 = reg$coef[2], adj_coef3 = reg$coef[3],
    corr = corr, corrw = corrw
  ))
}

eig_imp <- function(pc, thresh, X) {
  i <- which(cumsum(pc$values) / sum(pc$values) >= thresh)[1]
  E <- pc$vectors[, 1:i]
  mask <- is.na(X)
  dat <- dplyr::bind_cols(tibble::tibble(X), tibble::as_tibble(E))
  mod <- lm(X ~ ., data = dat)
  rsq <- summary(mod)$adj.r.squared
  p <- predict(mod, dat)
  co <- cor(p, X, use = "pair")

  temp <- tibble::tibble(X = p, mask, pos = seq_along(X))
  return(list(dat = temp, rsq = rsq, cor = co, ncomp = i))
}

eig_imp_lasso <- function(pc, thresh, X) {
  i <- which(cumsum(pc$values) / sum(pc$values) >= thresh)[1]
  E <- pc$vectors[, 1:i]
  mask <- is.na(X)
  E1 <- E[!mask, ]
  X1 <- X[!mask]
  # Check for constant or too few values
  if (length(unique(X1)) <= 1) {
    warning("eig_imp_lasso: y is constant or has only one value; skipping glmnet")
    return(NULL)
  }
  mod <- glmnet::glmnet(x = E1, y = X1, alpha = 0.5, family = "gaussian")
  cv_fit <- glmnet::cv.glmnet(E1, X1, alpha = 0.5, family = "gaussian")
  best_lambda <- cv_fit$lambda.min
  p <- predict(mod, E, s = best_lambda)
  co <- cor(p, X, use = "pair")

  temp <- tibble::tibble(X = p, mask, pos = seq_along(X))
  return(list(dat = temp, cor = co, ncomp = i))
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
#' If imputed variant is more significant the non-imputed rows, drop that variant.
#' Also, if the imputed variant has an illegal value, drop that variant.
#' @return filtered gwas
filter_imputation_results <- function(gwas, ld_matrix, min_bp, max_bp) {
  only_keep_inside_gwas_range <- gwas$BP > min_bp & gwas$BP < max_bp
  gwas <- gwas[only_keep_inside_gwas_range, ]

  snps_to_remove <- which(gwas$SE <= 0)
  if (length(snps_to_remove) > 0) {
    message(glue::glue("Removing {length(snps_to_remove)} SNPs with SE <= 0"))
  }

  min_original_p_value <- min(gwas$P[gwas$IMPUTED == F], na.rm = T)
  snps_to_investigate <- which(gwas$IMPUTED == T & gwas$P < min(min_original_p_value, lowest_p_value_threshold))
  for (snp_location in snps_to_investigate) {
    ld_correlations <- which(c(ld_matrix[snp_location, ]) > p_value_filter_correlation_threshold)
    ld_correlations <- ld_correlations[ld_correlations != snp_location]
    gwas_correlations <- gwas[(seq_len(nrow(gwas)) %in% ld_correlations) & gwas$IMPUTED == F, ]

    snp <- gwas[snp_location, ]

    if (length(gwas_correlations$P) > 0 && min(gwas_correlations$P, na.rm = T) * 0.1 > snp$P) {
      snps_to_remove <- c(snps_to_remove, snp_location)
    } else if (length(gwas_correlations$P) == 0 && min_original_p_value > snp$P) {
      snps_to_remove <- c(snps_to_remove, snp_location)
    }
  }

  removed_gwas <- NA
  updated_gwas <- gwas
  if (length(snps_to_remove) > 0) {
    updated_gwas <- gwas[-snps_to_remove, ]
    removed_gwas <- gwas[snps_to_remove, , drop = F]
  }

  return(list(
    gwas = updated_gwas,
    removed_gwas = removed_gwas,
    significant_rows_imputed = length(snps_to_investigate),
    significant_rows_filtered = length(snps_to_remove),
    bps_to_remove = gwas$BP[snps_to_remove]
  ))
}
