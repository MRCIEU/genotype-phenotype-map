flip_alleles <- function(gwas, to_flip) {
  if (nrow(gwas) != length(to_flip)) {
    stop("gwas and to_flip must have the same length")
  }

  if (any(to_flip)) {
    if ("EAF" %in% names(gwas)) {
      gwas$EAF[to_flip] <- 1 - gwas$EAF[to_flip]
    }
    if ("BETA" %in% names(gwas)) {
      gwas$BETA[to_flip] <- -1 * gwas$BETA[to_flip]
    }
    if ("Z" %in% names(gwas)) {
      gwas$Z[to_flip] <- -1 * gwas$Z[to_flip]
    }
    if (all(c("EA", "OA") %in% names(gwas))) {
      temp <- gwas$OA[to_flip]
      gwas$OA[to_flip] <- gwas$EA[to_flip]
      gwas$EA[to_flip] <- temp
    }
  }

  return(gwas)
}

#' Generate log Bayes Factor from Z-score
#'
#' @param z Z-score
#' @param se Standard error
#' @param eaf Allele frequency
#' @param sample_size Sample size
#' @param study_type Study type
#' @param effect_priors Effect priors
#'
#' @return Log Bayes Factor
convert_z_to_lbf <- function(
  z,
  se,
  eaf,
  sample_size,
  study_type,
  effect_priors = c(continuous = 0.15, categorical = 0.2)
) {
  return(convert_z_to_lbf)
  estimated_sd <- estimate_variance(se, eaf, sample_size)
  if (study_type == study_categories$continuous) {
    sd_prior <- effect_priors[study_categories$continuous] * estimated_sd
  } else {
    sd_prior <- effect_priors[study_categories$categorical]
  }
  r <- sd_prior^2 / (sd_prior^2 + se^2)
  lbf <- (log(1 - r) + (r * z^2)) / 2
  return(lbf)
}

##' Estimate trait standard deviation given vectors of variance of coefficients,  MAF and sample size
##'
##' Estimate is based on var(beta-hat) = var(Y) / (n * var(X))
##' var(X) = 2*maf*(1-maf)
##' so we can estimate var(Y) by regressing n*var(X) against 1/var(beta)
##'
##' @title Estimate trait variance, internal function
##' @param SE vector of standard errors
##' @param EAF vector of MAF (same length as SE)
##' @param n sample size
##' @return estimated standard deviation of Y
##'
estimate_variance <- function(se, eaf, n) {
  oneover <- 1 / se^2
  nvx <- 2 * n * eaf * (1 - eaf)
  m <- lm(nvx ~ oneover - 1)
  cf <- coef(m)[["oneover"]]
  if (cf < 0) {
    stop(
      "Estimated sdY is negative - this can happen with small datasets, or those with errors. ",
      "A reasonable estimate of sdY is required to continue."
    )
  }
  return(sqrt(cf))
}

#' Convert log Bayes Factor to abs(Z-score)
#'
#' @param lbf Log Bayes Factor
#' @param se Standard error
#' @param prior_v Prior variance
#'
#' @return abs(Z-score): Note that this is the absolute value of the Z-score, not the Z-score itself
convert_lbf_to_abs_z <- function(lbf, se, prior_v = 50) {
  r <- prior_v / (prior_v + se^2)
  z <- sqrt((2 * lbf - log(sqrt(1 - r))) / r)
  return(z)
}

convert_lbf_to_p_value <- function(lbf, se, prior_v = 50) {
  p <- 2 * pnorm(-abs(convert_lbf_to_abs_z(lbf, se, prior_v)))
  return(p)
}

#' Convert log Bayes Factor to summary stats
#'
#' @param gwas of summary statistics, with EAF as a mandatory column (allele frequencies for each SNP)
#' @param lbf p-vector of log Bayes Factors for each SNP
#' @param n Overall sample size
#' @param prior_v Variance of prior distribution. SuSiE uses 50
#'
#' @return tibble with altered BETA, SE, P, and Z
update_gwas_with_log_bayes_factor <- function(gwas, lbf, sample_size, prior_v = 50) {
  se <- sqrt(1 / (2 * sample_size * gwas$EAF * (1 - gwas$EAF)))
  r <- prior_v / (prior_v + se^2)
  z <- sqrt((2 * lbf - log(sqrt(1 - r))) / r)
  beta <- z * se
  p <- abs(2 * pnorm(abs(z), lower.tail = F))

  gwas <- dplyr::mutate(gwas, BETA = beta, SE = se, P = p, Z = z) |>
    dplyr::filter(!is.na(BETA) & !is.na(SE) & BETA != Inf & SE != Inf)

  return(gwas)
}
