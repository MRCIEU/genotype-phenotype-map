source('constants.R')

#https://explodecomputer.github.io/lab-book/posts/2024-09-21-summary-imputation/
imputation_correlation_threshold <- 0

parser <- argparser::arg_parser('Impute GWASes for pipeline')
parser <- argparser::add_argument(parser, '--ld_block', help = 'LD block that the ', type = 'character')
parser <- argparser::add_argument(parser, '--completed_output_file', help = 'Completed output file', type = 'character')
args <- argparser::parse_args(parser)

main <- function() {
  ld_info <- ld_block_dirs(args$ld_block)
  ld_matrix <- vroom::vroom(paste0(ld_info$ld_reference_panel_prefix, '.unphased.vcor1'), col_names = F, show_col_types = F) |>
    data.matrix()
  ld_region_from_reference_panel <- vroom::vroom(paste0(ld_info$ld_reference_panel_prefix, '.tsv'), show_col_types = F)

  standardised_studies_file <- paste0(ld_info$ld_block_data, '/standardised_studies.tsv')
  standardised_studies  <- vroom::vroom(standardised_studies_file , show_col_types = F)

  imputed_studies_file <- paste0(ld_info$ld_block_data, '/imputed_studies.tsv')
  if (file.exists(imputed_studies_file)) {
    existing_imputed_studies <- vroom::vroom(imputed_studies_file,
                                                  show_col_types = F,
                                                  col_types = vroom::cols(
                                                    chr = vroom::col_character(),
                                                    bp = vroom::col_number(),
                                                    sample_size = vroom::col_number(),
                                                    p_value_threshold = vroom::col_number(),
                                                    time_taken=vroom::col_character()
                                                  )
    )
  } else {
    existing_imputed_studies <- empty_imputed_studies()
  }

  if (nrow(standardised_studies) > 0) {
    imputed_studies <- apply(standardised_studies, 1, function (study) {
      start_time <- Sys.time()
      imputed_file <- sub('standardised', 'imputed', study[['file']])

      # if (imputed_file %in% existing_imputed_studies$file) {
      #   return()
      # }

      gwas <- vroom::vroom(study['file'], show_col_types = F)
      clumped_snps <- vroom::vroom(glue::glue('{extracted_study_dir}{study["study"]}/clumped_snps.tsv'))

      gwas_to_impute <- dplyr::left_join(
        dplyr::select(ld_region_from_reference_panel, -EAF),
        dplyr::select(gwas, -CHR, -BP, -EA, -OA),
        by=dplyr::join_by(SNP)
      ) 

      rows_to_impute <- !ld_region_from_reference_panel$SNP %in% gwas$SNP
      print(glue::glue('rows to impute: {sum(rows_to_impute)}'))
      print(glue::glue('rows in gwas: {nrow(gwas)}'))
      print(glue::glue('rows in ld matrix: {nrow(ld_region_from_reference_panel)}'))
      gwas_to_impute$EAF[rows_to_impute] <- ld_region_from_reference_panel$EAF[rows_to_impute]
      print(gwas[, c(9, 6, 7)])
      print(ld_region_from_reference_panel[, c(1,2,5)])
      print(gwas_to_impute[, c(2, 5, 6, 7)])

      print(glue::glue('gwas order matches ld matrix: {all(ld_region_from_reference_panel$SNP == gwas_to_impute$SNP)}'))

      # if we want to use all the clumped hits from the extraction phase...
      # given clumped snps, find the corresponding row numbers in the gwas
      # clumped_snp_index <- dplyr::mutate(gwas_to_impute, rn = dplyr::row_number()) |>
      #   dplyr::filter(BP %in% clumped_snps$BP) |>
      #   dplyr::group_by(BP) |>
      #   dplyr::arrange(P) |>
      #   dplyr::filter(dplyr::row_number() == 1)

      # if we want to use just the top hit in the region...
      # print(min(gwas_to_impute$P, na.rm=T))
      # clumped_snp_index <- which(gwas_to_impute$P == min(gwas_to_impute$P, na.rm = T))

      # if we want to re-clump from the ld ref panel
      clumped_snp_index <- clump_ld_region(gwas_to_impute$Z, ld_matrix)
      print('clumped snps:')
      print(gwas_to_impute[clumped_snp_index, ])

      result <- perform_imputation(gwas_to_impute, ld_matrix, clumped_snp_index)
      print(result$gwas[, c(2,6,7,8,11,12)])
      print(result$gwas[is.na(gwas_to_impute$BETA), c(2,6,7,8,11,12)])
      print(result)

      if(result$b_cor >= imputation_correlation_threshold) {
        vroom::vroom_write(result$gwas, imputed_file)
      } else {
        vroom::vroom_write(gwas, imputed_file)
      }

      time_taken <- hms::as_hms(difftime(Sys.time(), start_time)) 
      print(time_taken)

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
        b_cor=result$b_cor,
        se_cor=result$se_cor,
        b_adj=result$b_adj,
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
perform_imputation <- function(gwas, ld_matrix, index) {
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
    D <- diag(sqrt(2 * af * (1 - af)))
    Di <- diag(1 / diag(D))

    # Get the conditional estimates of the index SNP effects
    if(length(index) == 1) {
        bhat2 <- b[index]
    } else {
        bhat2 <- D[index,index] %*% MASS::ginv(ld_matrix[index,index]) %*% Di[index,index] %*% b[index]
    }
    b2 <- rep(0, nsnp)
    b2[index] <- bhat2

    # Get the simulated effect sizes
    betahat_sim <- as.numeric(Di %*% ld_matrix %*% D %*% b2)

    # Initialise the SE - this doesn't account for var(y) or sample size, but those are constants that can be obtained from regression re-scaling
    sehat <- sqrt(diag(Di))

    # Re-scale effect sizes and standard errors
    beta_scale <- lm(betahat_sim ~ b)$coef[2]
    se_scale <- lm(sehat ~ se)$coef[2]

    # Performance metrics
    b_cor <- cor(b, betahat_sim, use="pair")
    se_cor <- cor(se, sehat, use="pair")

    # Re-scale
    betahat_sim <- betahat_sim / beta_scale
    sehat <- sehat / se_scale

    gwas$BETA_IMPUTED <- betahat_sim
    gwas$SE_IMPUTED <- sehat
    stopifnot(all(!is.na(gwas$BETA_IMPUTED)))
    stopifnot(all(!is.na(gwas$SE_IMPUTED)))

    gwas$BETA[to_impute] <- gwas$BETA_IMPUTED[to_impute]
    gwas$SE[to_impute] <- gwas$SE_IMPUTED[to_impute]
    gwas$Z[to_impute] <- gwas$BETA_IMPUTED[to_impute] / gwas$SE_IMPUTED[to_impute]
    gwas$P[to_impute] <- 2 * pnorm(-abs(gwas$Z[to_impute]))

    return(
      list(
        gwas = gwas,
        b_adj = beta_scale,
        se_adj = se_scale,
        b_cor = b_cor,
        se_cor = se_cor,
        indices = length(index),
        rows_imputed = num_to_impute
      )
    )
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
      b_adj=numeric(),
      se_adj=numeric(),
      time_taken=character()
    )
  )
}

main()