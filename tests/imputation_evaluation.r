library(dplyr)
library(MASS)
library(data.table)
library(glue)
library(tidyr)
library(ggplot2)
library(ggpointdensity)
library(parallel)
library(exvatools)

clump_gwas <- function(z, R, zthresh = qnorm(1.5e-4, low=F), rthresh = 0.01) {
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

make_drd <- function(ld_matrix, af) {
    D <- diag(sqrt(2 * af * (1 - af)))
    Di <- diag(1 / diag(D))
    DRD <- Di %*% ld_matrix %*% D
    DRD
}

make_drd2 <- function(ld_matrix, af) {
    D <- diag(sqrt(2 * af * (1 - af)))
    Di <- diag(1 / diag(D))
    # DRD <- Di %*% ld_matrix %*% D
    DRD <- exvatools::multd(exvatools::dmult(Di, ld_matrix), D)
    DRD
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
perform_imputation <- function(gwas, ld_matrix, index, eval_frac=0.25) {
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
    # betahat_sim <- as.numeric(Di %*% ld_matrix %*% D %*% b2)
    # betahat_sim <- as.numeric(DRD %*% b2)
    drd <- exvatools::multd(exvatools::dmult(Di, ld_matrix), D)
    betahat_sim <- as.numeric(drd %*% b2)

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
        b_corr_top = b_adj$corrw,
        se_cor = se_adj$corr,
        se_corr_top = se_adj$corrw,
        indices = length(index),
        rows_imputed = num_to_impute,
        rows_region = nrow(gwas)
      )
    )
}

organise_trait <- function(ssfile, ldvars, af) {
    sss <- fread(ssfile)
    print(table(sss$SNP %in% ldvars))
    sss <- subset(sss, SNP %in% ldvars)

    smiss <- tibble(SNP = ldvars[!ldvars %in% sss$SNP], BETA=NA, SE=NA, Z=NA, P=NA, RSID=SNP) %>% 
        tidyr::separate(SNP, into=c("CHR", "BP", "EA", "OA"), remove=FALSE) %>%
        mutate(CHR=as.numeric(CHR), BP=as.numeric(BP))

    s <- bind_rows(sss, smiss)
    m <- match(ldvars, s$SNP)
    stopifnot(all(s$SNP[m] == ldvars))

    s <- s[m,]
    stopifnot(all(s$SNP == ldvars))

    i <- which(is.na(s$EAF))
    if(length(i) > 0) {
        stopifnot(all(s$SNP[i] %in% af$ID))
        m <- match(s$SNP[i], af$ID)
        stopifnot(all(s$SNP[i] == af$ID[m]))
        s$EAF[i] <- af$ALT_FREQS[m]
    }
    s
}

generate_ldmatrix <- function(bfile, out="test") {
    glue("plink2 --bfile {bfile} --r-unphased square ref-based --maf 0.01 --out {out}") %>% system()
    ld <- fread(paste0(out, ".unphased.vcor1")) %>% as.matrix()
    af <- fread(glue("{bfile}.afreq"))
    ldvars <- scan(paste0(out, ".unphased.vcor1.vars"), what="character")
    return(list(ld=ld, af=af, ldvars=ldvars))
}

get_region <- function(sfile) {
    x <- sfile %>% gsub(".tsv.gz", "", .) %>% basename %>% {strsplit(., "_")[[1]]}
    d <- file.path("/local-scratch/projects/genotype-phenotype-map/data/ld_reference_panel_hg38", x[1], x[2])
    pos <- as.numeric(x[3])
    fn <- list.files(d, full.names=TRUE) %>% grep(".bed$", ., value=TRUE) %>% gsub(".bed", "", .)
    r <- tibble(fn, r=basename(fn)) %>% tidyr::separate(r, into=c("start", "end"), remove=FALSE) %>% mutate(start=as.numeric(start), end=as.numeric(end))
    f <- subset(r, pos >= start & pos <= end)$fn
    ld <- fread(paste0(f, ".unphased.vcor1")) %>% as.matrix()
    ldvars <- scan(paste0(f, ".unphased.vcor1.vars"), what="character")
    af <- fread(paste0(f, ".afreq"))
    pc <- readRDS(paste0(f, ".ldeig.rds"))
    return(list(ld=ld, ldvars=ldvars, af=af, pc=pc))
}

eig_imp_iterative <- function(pc, thresh, niter, X) {
    i <- which(cumsum(pc$sdev) / sum(pc$sdev) >= thresh)[1]
    E <- pc$loadings[,1:i]
    mask <- is.na(X)
    X[is.na(X)] <- mean(X, na.rm=T)
    for(j in 1:niter) {
        X_old <- X
        Y <- t(E) %*% X
        X2 <- E %*% Y
        X[mask] <- X2[mask]
        maxdif <- max(abs(X_old - X))
        print(maxdif)
    }

    temp <- tibble(X, mask, pos=1:length(X))
    temp
}

eig_imp <- function(pc, thresh, X) {
    message("Getting threshold")
    i <- which(cumsum(pc$values) / sum(pc$values) >= thresh)[1]
    E <- pc$vectors[,1:i]
    mask <- is.na(X)
    dat <- bind_cols(tibble(X), as_tibble(E))
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

    temp <- tibble(X=p, mask, pos=1:length(X))
    return(list(dat=temp, rsq=rsq, cor=co, ncomp=i))
}

eig_imp_ve <- function(pc, thresh, X) {
    message("Getting threshold")
    ve <- sapply(1:nrow(E), \(i) cumsum(E[i,]^2))
    i <- which(apply(ve, 1, median) > thresh)[1]
    # i <- which(cumsum(pc$values) / sum(pc$values) >= thresh)[1]
    E <- pc$vectors[,1:i]
    mask <- is.na(X)
    dat <- bind_cols(tibble(X), as_tibble(E))
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

    temp <- tibble(X=p, mask, pos=1:length(X))
    return(list(dat=temp, rsq=rsq, cor=co, ncomp=i))
}

set_se_outliers_missing <- function(se, se_imputed, outthresh = 3) {
    mod <- lm(se ~ se_imputed)
    cd <- cooks.distance(mod)
    i <- outlier_detection(cd, outthresh)
    return(i)
}

#' Basic imputation function
#' 
#' @param gwas A data frame with columns betahat2 = vector of effect estimates in the same order as ld_matrix and with NAs for variants that need to be imputed;
#'  se = as with betahat2 but for available standard errors, af = allele frequencies (no missing values allowed, so use reference panel if there are missing values)
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

setup_data_for_tests <- function(sfile, maxsnps=5000, remove_missing=TRUE) {
    ld <- get_region(sfile)
    str(ld)
    pc <- ld$pc
    af <- ld$af
    ldvars <- ld$ldvars
    ld <- ld$ld
    s <- organise_trait(sfile, ldvars, af)
    if(nrow(s) > maxsnps) {
        s$BETA[maxsnps:nrow(s)] <- NA
    }
    if(remove_missing) {
        i <- is.na(s$BETA)
        print(table(i))
        ld <- ld[!i, !i]
        ldvars <- ldvars[!i]
        af <- af[af$ID %in% ldvars,]
        s <- s[!i,]
        stopifnot(all(af$ID == s$SNP))
        dim(pc$vectors)
        dim(s)
        head(colnames(pc$vectors))
        pc$vectors <- pc$vectors[!i,]
    }
    return(list(s=s, af=af, ld=ld, ldvars=ldvars, pc=pc))
}

add_missing_data <- function(testdat, frac, seed=1234) {
    set.seed(seed)
    m <- rbinom(nrow(testdat$s), 1, frac) %>% as.logical()
    testdat$s$BETA_TRUE <- testdat$s$BETA
    testdat$s$SE_TRUE <- testdat$s$SE
    testdat$s$Z_TRUE <- testdat$s$Z
    testdat$s$P_TRUE <- testdat$s$P
    testdat$s$BETA[m] <- NA
    testdat$s$SE[m] <- NA
    testdat$s$Z[m] <- NA
    testdat$s$P[m] <- NA
    return(testdat)
}

run_test_th <- function(testdat, param) {
    param$bcor <- NA
    i <- 1
    for(i in 1:nrow(param)) {
        td <- add_missing_data(testdat, param$frac[i])
        ind <- clump_gwas(td$s$Z, td$ld, param$zthresh[i], param$rthresh[i])
        message("Number of index variants: ", length(ind))
        o1 <- perform_imputation(td$s, td$ld, ind)
        m <- is.na(td$s$BETA)
        ev <- cor_out(o1$gwas$BETA[m], o1$gwas$BETA_TRUE[m])
        param$nout[i] <- ev$nout
        param$bcor[i] <- ev$c
        param$nmiss[i] <- sum(m)
    }
    return(param)
}


cor_out <- function(a, b, pl=NULL) {
    cd <- cooks.distance(lm(a ~ b))
    i <- outlier_detection(cd)
    a[i] <- NA
    b[i] <- NA
    if(!is.null(pl)) {
        png(pl)
        plot(a, b)
        abline(0,1)
        dev.off()
    }
    return(list(c=cor(a, b, use="pair"), nout=sum(i)))
}


run_test_eig <- function(testdat, param) {
    param$bcor <- NA
    for(i in 1:nrow(param)) {
        td <- add_missing_data(testdat, param$frac[i])
        # o1 <- perform_imputation(td$s, td$ld, clump_gwas(td$s$Z, td$ld, param$zthresh[i], param$rthresh[i]))
        o1 <- perform_imputation_eig(td$s, td$pc, param$eigthresh[i])
        m <- is.na(td$s$BETA)
        ev <- cor_out(o1$gwas$BETA[m], o1$gwas$BETA_TRUE[m])
        param$nout[i] <- ev$nout
        param$bcor[i] <- ev$c
        param$nmiss[i] <- sum(m)
    }
    return(param)
}


###########################################

# Get list of files to test

a <- list.files("/local-scratch/projects/genotype-phenotype-map/test/data/study/ebi-a-GCST90002304/standardised/", full.names=TRUE) %>% sample(100, replace=F)



###########################################

# Perform analysis on one region

sfile <- a[1]
testdat <- setup_data_for_tests(A, maxsnps = 100000, remove_missing=FALSE)
str(testdat)
o <- perform_imputation_eig(testdat$s, testdat$pc, 0.9)
o

ggsave(
    ggplot(o$gwas %>% filter(!imp), aes(BETA, BETA_IMPUTED)) +
    geom_point(alpha=0.5) +
    geom_abline(),
    file="temp.png"
)

ggsave(
    ggplot(o$gwas %>% filter(!imp), aes(SE, SE_IMPUTED)) +
    geom_point(alpha=0.5) +
    geom_abline(),
    file="temp.png"
)


###########################################

# Enumerate over all regions, methods and paramters

param <- expand.grid(
    sfile = a,
    method = "th",
    rthresh = c(0.01, 0.1),
    zthresh = qnorm(c(1e-4), low=F),
    frac = c(0.1, 0.5)
) %>% mutate(sfile = as.character(sfile))

res_th <- mclapply(a, \(A) {
    testdat <- try(setup_data_for_tests(A))
    if(inherits(testdat, "try-error")) {return(NULL)}
    pt <- run_test_th(testdat, subset(param, sfile==A))    
}, mc.cores=50)

res_th <- lapply(res_th, \(x) {
    if(inherits(x, "data.frame")) {
         return (x)
    } else { 
        return(NULL)
    }
}) %>% bind_rows() %>% as_tibble


param <- expand.grid(
    sfile = a,
    method = "eig",
    eigthresh = c(0.5, 0.6, 0.7, 0.8, 0.9, 0.99),
    frac = c(0.1, 0.5)
) %>% mutate(sfile = as.character(sfile))

res_eig <- mclapply(a, \(A) {
    testdat <- try(setup_data_for_tests(A))
    if(class(testdat) == "try-error") { return(NULL) }
    pt <- run_test_eig(testdat, subset(param, sfile==A))    
}, mc.cores=50) %>% bind_rows()

res <- bind_rows(res_eig, res_th) %>% as_tibble()


#########################################



res_best <- res %>% arrange(desc(bcor)) %>% group_by(method, sfile, frac) %>% slice_head(n=1)
res_best

table(res_best$eigthresh)
table(res_best$rthresh)
table(res_best$zthresh)

# Which is the best eig threshold

res_best2 <- res %>% arrange(desc(bcor)) %>% group_by(sfile, frac) %>% slice_head(n=1)
res_best2
table(res_best2$method)

summary(res_best2$bcor)


#########################################

ggsave(
    ggplot(res, aes(x=method, y=bcor)) +
    geom_boxplot(aes(fill=as.factor(eigthresh))) +
    facet_grid(. ~ frac, label=label_both),
    file="method_eval.png"
)

