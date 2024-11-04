source('constants.R')

imputation_correlation_threshold <- 0.7
p_value_filter_correlation_threshold <- 0.6

parser <- argparser::arg_parser('Impute GWASes for pipeline')
parser <- argparser::add_argument(parser, '--ld_block', help = 'LD block that the ', type = 'character')
parser <- argparser::add_argument(parser, '--completed_output_file', help = 'Completed output file', type = 'character')
args <- argparser::parse_args(parser)

main <- function() {
  ld_info <- ld_block_dirs(args$ld_block)
  ld_matrix_info <- vroom::vroom(glue::glue('{ld_info$ld_reference_panel_prefix}.tsv'), show_col_types = F)
  ld_matrix <- vroom::vroom(glue::glue('{ld_info$ld_reference_panel_prefix}.unphased.vcor1'), col_names=F, show_col_types = F)
  ld_matrix_eig <- readRDS(glue::glue('{ld_info$ld_reference_panel_prefix}.ldeig.rds'))

  standardised_studies_file <- glue::glue('{ld_info$ld_block_data}/standardised_studies.tsv')
  standardised_studies  <- vroom::vroom(standardised_studies_file , show_col_types = F)

  imputed_studies_file <- glue::glue('{ld_info$ld_block_data}/imputed_studies.tsv')
  if (file.exists(imputed_studies_file)) {
    existing_imputed_studies <- vroom::vroom(imputed_studies_file, show_col_types = F, col_types = imputed_column_types)
  } else {
    existing_imputed_studies <- empty_imputed_studies()
  }

  if (nrow(standardised_studies) > 0) {
    imputed_studies <- apply(standardised_studies, 1, function (study) {
      start_time <- Sys.time()
      imputed_file <- sub('standardised', 'imputed', study[['file']])

      if (imputed_file %in% existing_imputed_studies$file) {
        return()
      }

      gwas <- vroom::vroom(study['file'], show_col_types = F)

      gwas_to_impute <- dplyr::left_join(
        dplyr::select(ld_matrix_info, -EAF),
        dplyr::select(gwas, -CHR, -BP, -EA, -OA),
        by=dplyr::join_by(SNP)
      ) 

      rows_to_impute <- !ld_matrix_info$SNP %in% gwas$SNP
      gwas_to_impute$EAF[rows_to_impute] <- ld_matrix_info$EAF[rows_to_impute]

      result <- perform_imputation(gwas_to_impute, ld_matrix_eig)

      pre_filter_file <- sub('.tsv.gz', '_pre_filter.tsv.gz', imputed_file)
      vroom::vroom_write(result$gwas, pre_filter_file)

      filtered_results <- filter_imputation_results(result$gwas, ld_matrix, min(gwas$BP), max(gwas$BP))

      if(result$b_cor >= imputation_correlation_threshold) {
        vroom::vroom_write(filtered_results$gwas, imputed_file)
      } else {
        vroom::vroom_write(gwas, imputed_file)
      }

      time_taken <- hms::as_hms(difftime(Sys.time(), start_time)) 

      imputation_info <- data.frame(
        study=study[['study']],
        file=imputed_file,
        ancestry=study['ancestry'],
        chr = as.character(study[['chr']]),
        bp = as.numeric(study[['bp']]),
        p_value_threshold = as.numeric(study[['p_value_threshold']]),
        category=study['category'],
        sample_size=as.numeric(study['sample_size']),
        cis_trans=study['cis_trans'],
        rows_imputed=result$rows_imputed,
        significant_rows_imputed=filtered_results$significant_rows_imputed,
        significant_rows_filtered=filtered_results$significant_rows_filtered,
        b_cor=result$b_cor,
        se_cor=result$se_cor,
        z_adj=result$z_adj,
        se_adj=result$se_adj,
        time_taken=as.character(time_taken)
      )

      return(imputation_info)
    }) |> dplyr::bind_rows()
  }

  if (nrow(imputed_studies) > 0) {
    imputed_studies <- dplyr::bind_rows(existing_imputed_studies, imputed_studies) |>
      dplyr::distinct()

    vroom::vroom_write(imputed_studies, imputed_studies_file)
  }

  vroom::vroom_write(data.frame(), args$completed_output_file)
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
#' - z_adj: The adjustment factor for the effect sizes
#' - se_adj: The adjustment factor for the standard errors
#' - b_cor: The correlation between the true and imputed effect sizes - this is critical for evaluation of the performance of the imputation,
#'      it should be close to 1 e.g > 0.7 would be a reasonable threshold
#' - se_cor: The correlation between the true and imputed standard errors
perform_imputation <- function(gwas, pc, thresh=0.9, eval_frac=0.25) {
    b <- gwas$BETA
    se <- gwas$SE
    z <- b/se
    af <- gwas$EAF
    to_impute <- is.na(b)
    num_to_impute <- sum(to_impute)

    if (num_to_impute == 0) {
      return(list(
          gwas = gwas, b_cor = NA, se_cor = NA, z_adj = NA, se_adj = NA, indices = NA, rows_imputed = 0
      ))
    }

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

adjust <- function(truth, predicted, eval_frac = 0.5, npoly=3) {
    outs <- outlier_detection(truth / predicted)
    reg <- lm(truth[!outs] ~ poly(predicted[!outs], npoly, raw=T))
    # adj <- predict(reg, newdata=data.frame(predicted=predicted))
    adj <- predicted * reg$coef[2] + predicted^2 * reg$coef[3] + predicted^3 * reg$coef[4] + reg$coef[1]
    corr <- cor(adj[!outs], truth[!outs], use="pair")
    iqr <- truth > quantile(truth, 1-(eval_frac/2), na.rm=T) | truth < quantile(truth, eval_frac/2, na.rm=T)
    corrw <- cor(adj[!outs & iqr], truth[!outs & iqr], use="pair")

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

#' Filter imputed results, and remove results with over-inflated p-vlaues
#' @param imputed gwas
#' @param ld correlation matrix
#' 
#' First, filters the imputation results to inside the range of the original gwas
#' Second, finds imputed variants with low pvalues, non-imputed rows that are correlated
#' If imputed variant is more significant the non-imputed rows, drop that variant
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
    if (min(gwas_correlations$P * 0.1) > snp$P) {
      snps_to_remove <- c(snps_to_remove, snp_location)
    }
  }

  if (length(snps_to_remove) > 0) gwas <- gwas[-snps_to_remove, ]
  return(list(gwas=gwas, significant_rows_imputed=length(snps_to_investigate), significant_rows_filtered=length(snps_to_remove)))
}

empty_imputed_studies <- function() {
  return(
    data.frame(
      study=character(),
      file=character(),
      ancestry=character(),
      chr=character(),
      bp=numeric(),
      p_value_threshold=numeric(),
      category=character(),
      sample_size=numeric(),
      cis_trans=character(),
      rows_imputed=numeric(),
      b_cor=numeric(),
      se_cor=numeric(),
      z_adj=numeric(),
      se_adj=numeric(),
      time_taken=character()
    )
  )
}

main()
