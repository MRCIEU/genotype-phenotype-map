library(dplyr)
library(MASS)
library(data.table)
library(glue)
library(tidyr)
library(ggplot2)
library(ggpointdensity)
library(parallel)
library(exvatools)

clump_gwas <- function(z, R, zthresh = qnorm(1.5e-4, low = F), rthresh = 0.01) {
  z <- abs(z)
  z[z < zthresh] <- NA
  k <- c()

  while (!all(is.na(z))) {
    i <- which.max(z)
    k <- c(k, i)
    z[i] <- NA
    z[which(R[i, ]^2 > rthresh)] <- NA
  }
  return(k)
}

outlier_detection <- function(r) {
  sd1 <- sd(r, na.rm = T)
  m <- median(r, na.rm = T)
  r[r < (m - 3 * sd1)] <- NA
  r[r > (m + 3 * sd1)] <- NA
  sd2 <- sd(r, na.rm = T)
  r[r < (m - 3 * sd2)] <- NA
  r[r > (m + 3 * sd2)] <- NA
  outliers <- is.na(r)
  return(outliers)
}

adjust <- function(truth, predicted, eval_frac = 0.5) {
  outs <- outlier_detection(truth / predicted)
  reg <- lm(truth[!outs] ~ predicted[!outs])
  adj <- predicted * reg$coef[2] - reg$coef[1]
  corr <- cor(adj[!outs], truth[!outs], use = "pair")
  iqr <- truth > quantile(truth, 1 - (eval_frac / 2), na.rm = T) | truth < quantile(truth, eval_frac / 2, na.rm = T)
  corrw <- cor(adj[!outs & iqr], truth[!outs & iqr], use = "pair")

  return(list(
    adj = adj,
    outliers = outs,
    adj_slope = reg$coef[2],
    adj_intercept = reg$coef[1],
    corr = corr,
    corrw = corrw
  ))
}

make_drd <- function(ld_matrix, af) {
  D <- diag(sqrt(2 * af * (1 - af)))
  Di <- diag(1 / diag(D))
  DRD <- Di %*% ld_matrix %*% D
  return(DRD)
}

make_drd2 <- function(ld_matrix, af) {
  D <- diag(sqrt(2 * af * (1 - af)))
  Di <- diag(1 / diag(D))
  # DRD <- Di %*% ld_matrix %*% D
  DRD <- exvatools::multd(exvatools::dmult(Di, ld_matrix), D)
  return(DRD)
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
perform_imputation <- function(gwas, ld_matrix, index, eval_frac = 0.25) {
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
  stopifnot(all(se > 0, na.rm = TRUE))

  # Calculate the diagonal matrix of variances and the inverse
  message("Generating D")
  D <- diag(sqrt(2 * af * (1 - af)))
  Di <- diag(1 / diag(D))

  # Get the conditional estimates of the index SNP effects
  message("Conditional index variant estimates")
  if (length(index) == 1) {
    bhat2 <- b[index]
  } else {
    bhat2 <- D[index, index] %*% MASS::ginv(ld_matrix[index, index]) %*% Di[index, index] %*% b[index]
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

  smiss <- tibble(SNP = ldvars[!ldvars %in% sss$SNP], BETA = NA, SE = NA, Z = NA, P = NA, RSID = SNP) %>%
    tidyr::separate(SNP, into = c("CHR", "BP", "EA", "OA"), remove = FALSE) %>%
    mutate(CHR = as.numeric(CHR), BP = as.numeric(BP))

  s <- bind_rows(sss, smiss)
  m <- match(ldvars, s$SNP)
  stopifnot(all(s$SNP[m] == ldvars))

  s <- s[m, ]
  stopifnot(all(s$SNP == ldvars))

  i <- which(is.na(s$EAF))
  if (length(i) > 0) {
    stopifnot(all(s$SNP[i] %in% af$ID))
    m <- match(s$SNP[i], af$ID)
    stopifnot(all(s$SNP[i] == af$ID[m]))
    s$EAF[i] <- af$ALT_FREQS[m]
  }
  return(s)
}

generate_ldmatrix <- function(bfile, out = "test") {
  glue("plink2 --bfile {bfile} --r-unphased square ref-based --maf 0.01 --out {out}") %>% system()
  ld <- fread(paste0(out, ".unphased.vcor1.gz")) %>% as.matrix()
  af <- fread(glue("{bfile}.afreq"))
  ldvars <- scan(paste0(out, ".unphased.vcor1.vars"), what = "character")
  return(list(ld = ld, af = af, ldvars = ldvars))
}

get_region <- function(sfile) {
  x <- sfile %>%
    gsub(".tsv.gz", "", .) %>%
    basename() %>%
    {
      strsplit(., "_")[[1]]
    }
  d <- file.path("/local-scratch/projects/genotype-phenotype-map/data/ld_reference_panel_hg38", x[1], x[2])
  pos <- as.numeric(x[3])
  fn <- list.files(d, full.names = TRUE) %>%
    grep(".bed$", ., value = TRUE) %>%
    gsub(".bed", "", .)
  r <- tibble(fn, r = basename(fn)) %>%
    tidyr::separate(r, into = c("start", "end"), remove = FALSE) %>%
    mutate(start = as.numeric(start), end = as.numeric(end))
  f <- subset(r, pos >= start & pos <= end)$fn
  ld <- fread(paste0(f, ".unphased.vcor1.gz")) %>% as.matrix()
  ldvars <- scan(paste0(f, ".unphased.vcor1.vars"), what = "character")
  af <- fread(paste0(f, ".afreq"))
  pc <- readRDS(paste0(f, ".ldeig.rds"))
  return(list(ld = ld, ldvars = ldvars, af = af, pc = pc))
}

# Example analysis

bfile <- "/local-scratch/projects/genotype-phenotype-map/data/ld_reference_panel_hg38/EUR/9/31590323-32738472"
sfile <- "/local-scratch/projects/genotype-phenotype-map/test/data/study/ebi-a-GCST90002304/standardised/EUR_9_32450189.tsv.gz" # nolint: line_length_linter.

l <- generate_ldmatrix(bfile)
ld <- l$ld
af <- l$af
ldvars <- l$ldvars

ld <- get_region(sfile)
af <- ld$af
ldvars <- ld$ldvars
ld <- ld$ld

s <- organise_trait(sfile, ldvars, af)
o <- perform_imputation(s, ld, clump_gwas(s$Z, ld))
o

# Plot performance
png("test.png")
plot(y = o$gwas$BETA[!o$gwas$imp], x = o$gwas$BETA_IMPUTED[!o$gwas$imp])
abline(a = 0, b = 1)
dev.off()


eig <- eigen(ld)
o <- eig_imp2(eig, 0.99, 30, s$BETA)


# Try different parameters
param <- expand.grid(
  zthresh = qnorm(c(1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8), low = F),
  rthresh = c(0.01, 0.05, 0.1, 0.2, 0.5),
  nind = NA,
  bcor = NA,
  bcortop = NA
)

or <- mclapply(seq_len(nrow(param)), \(i) {
  ind <- clump_gwas(s$Z, ld, zthresh = param$zthresh[i], rthresh = param$rthresh[i])
  o2 <- perform_imputation(s, ld, ind)
  return(param[i, ] %>% mutate(nind = length(ind), bcor = o2$b_cor, bcortop = o2$b_corr_top))
}, mc.cores = 30) %>% bind_rows()

p1 <- ggplot(or, aes(x = as.factor(zthresh), y = as.factor(rthresh))) +
  geom_tile(aes(fill = bcortop)) +
  theme(axis.text.x = element_text(angle = 90))
ggsave(p1, file = "temp2.png")

## Eigendecomposition method

pc <- princomp(ld)

eig <- eigen(ld)


pc <- readRDS("pc2.rds")
class(pc)
dim(pc[[2]])


i <- which(cumsum(pc$sdev) / sum(pc$sdev) >= 0.99)[1]
i

which(cumsum(eig$values) / sum(eig$values) >= 0.99)[1]


E <- pc$loadings[, 1:i]
dim(E)

X <- s$BETA
X[is.na(X)] <- mean(X, na.rm = T)

for (j in 1:30) {
  X_old <- X
  Y <- t(E) %*% X
  X2 <- E %*% Y
  mask <- is.na(s$BETA)
  X[mask] <- X2[mask]
  maxdif <- max(abs(X_old - X))
  print(maxdif)
}

temp <- tibble(X, mask, pos = s$BP)

p1 <- ggplot(temp, aes(pos, X)) +
  geom_point(aes(colour = mask))
ggsave(p1, file = "temp2.png")

# Let's generate known missing values

thresh <- 0.99
eig_imp <- function(pc, thresh, niter, X) {
  i <- which(cumsum(pc$sdev) / sum(pc$sdev) >= thresh)[1]
  E <- pc$loadings[, 1:i]
  mask <- is.na(X)
  X[is.na(X)] <- mean(X, na.rm = T)
  for (j in 1:niter) {
    X_old <- X
    Y <- t(E) %*% X
    X2 <- E %*% Y
    X[mask] <- X2[mask]
    maxdif <- max(abs(X_old - X))
    print(maxdif)
  }

  temp <- tibble(X, mask, pos = 1:length(X))
  return(temp)
}

eig_imp2 <- function(pc, thresh, niter, X) {
  i <- which(cumsum(pc$values) / sum(pc$values) >= thresh)[1]
  E <- pc$vectors[, 1:i]
  mask <- is.na(X)
  X[is.na(X)] <- mean(X, na.rm = T)
  for (j in 1:niter) {
    X_old <- X
    Y <- t(E) %*% X
    X2 <- E %*% Y
    X[mask] <- X2[mask]
    maxdif <- max(abs(X_old - X))
    print(maxdif)
  }

  temp <- tibble(X, mask, pos = 1:length(X))
  return(temp)
}


mask2 <- rbinom(length(X), 1, 0.1) %>% as.logical()
X2 <- X
X2[mask2] <- NA
table(is.na(X2))
temp <- eig_imp(pc, 0.99, 30, X2)
temp$Xorig <- X

cor(temp$X[temp$mask], temp$Xorig[temp$mask])

mask2 <- rbinom(length(X), 1, 0.3) %>% as.logical()
X2 <- X
X2[mask2] <- NA
table(is.na(X2))
temp <- eig_imp2(eig, 1, 30, X2)
temp$Xorig <- X

summary()

summary(rowSums((pc$loadings^2)[, 1:300]))
summary(rowSums((pc$loadings^2)))


cor(temp$X[temp$mask], temp$Xorig[temp$mask])

png("temp6.png")
plot(temp$X[temp$mask], temp$Xorig[temp$mask])
abline(0, 1)
dev.off()

saveRDS(pc, file = "pc.rds")

a <- list(sdev = pc$sdev, loadings = pc$loadings)
str(a)
saveRDS(a, file = "pc2.rds")
a2 <- readRDS("pc2.rds")
str(a2)


bfile <- "/local-scratch/projects/genotype-phenotype-map/data/ld_reference_panel_hg38/EUR/9/31590323-32738472"
sfile <- "/local-scratch/projects/genotype-phenotype-map/test/data/study/ebi-a-GCST90002304/standardised/EUR_9_32450189.tsv.gz" # nolint: line_length_linter.

l <- generate_ldmatrix(bfile)
ld <- l$ld
af <- l$af
ldvars <- l$ldvars

s <- organise_trait(sfile, ldvars, af)

pc <- readRDS("pc2.rds")
str(pc)
names(pc) <- c("values", "vectors")

thresh <- 0.99
X <- s$BETA
table(is.na(X))
eig_imp3 <- function(pc, thresh, X) {
  message("Getting threshold")
  i <- which(cumsum(pc$values) / sum(pc$values) >= thresh)[1]
  E <- pc$vectors[, 1:i]
  mask <- is.na(X)
  dat <- bind_cols(tibble(X), as_tibble(E))
  message("Fitting model")
  mod <- lm(X ~ ., data = dat)
  rsq <- summary(mod)$adj.r.squared
  message("Imputing")
  p <- predict(mod, dat)
  co <- cor(p, X, use = "pair")

  # png("temp2.png")
  # plot(p, X)
  # abline(0,1)
  # dev.off()

  temp <- tibble(X = p, mask, pos = 1:length(X))
  return(list(dat = temp, rsq = rsq, cor = co, ncomp = i))
}

eig_imp4 <- function(pc, thresh, X) {
  message("Getting threshold")
  ve <- sapply(seq_len(nrow(E)), \(i) return(cumsum(E[i, ]^2)))
  i <- which(apply(ve, 1, median) > thresh)[1]
  # i <- which(cumsum(pc$values) / sum(pc$values) >= thresh)[1]
  E <- pc$vectors[, 1:i]
  mask <- is.na(X)
  dat <- bind_cols(tibble(X), as_tibble(E))
  message("Fitting model")
  mod <- lm(X ~ ., data = dat)
  rsq <- summary(mod)$adj.r.squared
  message("Imputing")
  p <- predict(mod, dat)
  co <- cor(p, X, use = "pair")

  # png("temp2.png")
  # plot(p, X)
  # abline(0,1)
  # dev.off()

  temp <- tibble(X = p, mask, pos = 1:length(X))
  return(list(dat = temp, rsq = rsq, cor = co, ncomp = i))
}

temp <- eig_imp3(pc, 0.99, s$BETA, 0.8)


table(temp$dat$mask)
table(is.na(temp$dat$X))


Y <- s$BETA
Y[is.na(Y)] <- temp$dat$X[is.na(Y)]
X <- Y
m <- rbinom(nrow(s), 1, 0.15) %>% as.logical()
X[m] <- NA
temp2 <- eig_imp3(pc, 0.99, X, 0.8)

temp2
cor(temp2$dat$X[m], Y[m])


temp3 <- eig_imp2(pc, 0.99, 30, X)
cor(temp3$X, temp2$dat$X)
cor(temp3$X[m], Y[m])

s2 <- s
s2$BETA <- X

o <- perform_imputation(s2, ld, clump_gwas(s$Z, ld, zthresh = qnorm(1e-4, low = F)))
cor(o$gwas$BETA[m], Y[m])


bfile <- "/local-scratch/projects/genotype-phenotype-map/data/ld_reference_panel_hg38/EUR/9/31590323-32738472"
sfile <- "/local-scratch/projects/genotype-phenotype-map/test/data/study/ebi-a-GCST90002304/standardised/EUR_9_32450189.tsv.gz" # nolint: line_length_linter.

l <- generate_ldmatrix(bfile)
ld <- l$ld
af <- l$af
ldvars <- l$ldvars

s <- organise_trait(sfile, ldvars, af)

i <- is.na(s$BETA)
ld <- ld[!i, !i]
ldvars <- ldvars[!i]
af <- af[af$ID %in% ldvars, ]
s <- s[!i, ]
stopifnot(all(af$ID == s$SNP))


pc <- eigen(ld)

m <- rbinom(nrow(s), 1, 0.15) %>% as.logical()
s2 <- s
s2$BETA[m] <- NA
s2$SE[m] <- NA
s2$Z[m] <- NA
s2$P[m] <- NA

o1 <- perform_imputation(s2, ld, clump_gwas(s2$Z, ld))
o2 <- eig_imp3(pc, 0.99, s2$BETA)
# o3 <- eig_imp4(pc, 0.3, s2$BETA)

cor(o1$gwas$BETA[m], s$BETA[m])
cor(o2$dat$X[m], s$BETA[m])
cor(o3$dat$X[m], s$BETA[m])





a <- list.files(
  "/local-scratch/projects/genotype-phenotype-map/test/data/study/ebi-a-GCST90002304/standardised/",
  full.names = TRUE
) %>%
  sample(100, replace = F)

sfile <- a[1]
sfile <- "/local-scratch/projects/genotype-phenotype-map/test/data/study/ebi-a-GCST90002304/standardised//EUR_3_169379112.tsv.gz" # nolint: line_length_linter.
setup_data_for_tests <- function(sfile) {
  ld <- get_region(sfile)
  str(ld)
  pc <- ld$pc
  af <- ld$af
  ldvars <- ld$ldvars
  ld <- ld$ld
  s <- organise_trait(sfile, ldvars, af)
  if (nrow(s) > 5000) {
    s$BETA[5000:nrow(s)] <- NA
  }
  i <- is.na(s$BETA)
  print(table(i))
  ld <- ld[!i, !i]
  ldvars <- ldvars[!i]
  af <- af[af$ID %in% ldvars, ]
  s <- s[!i, ]
  stopifnot(all(af$ID == s$SNP))
  dim(pc$vectors)
  dim(s)
  head(colnames(pc$vectors))
  pc$vectors <- pc$vectors[!i, ]
  return(list(s = s, af = af, ld = ld, ldvars = ldvars, pc = pc))
}

add_missing_data <- function(testdat, frac, seed = 1234) {
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

testdat <- setup_data_for_tests(sfile)
str(testdat)

run_test_th <- function(testdat, param) {
  param$bcor <- NA
  i <- 1
  for (i in seq_len(nrow(param))) {
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


cor_out <- function(a, b, pl = NULL) {
  cd <- cooks.distance(lm(a ~ b))
  i <- outlier_detection(cd)
  a[i] <- NA
  b[i] <- NA
  if (!is.null(pl)) {
    png(pl)
    plot(a, b)
    abline(0, 1)
    dev.off()
  }
  return(list(c = cor(a, b, use = "pair"), nout = sum(i)))
}


run_test_eig <- function(testdat, param) {
  param$bcor <- NA
  for (i in seq_len(nrow(param))) {
    td <- add_missing_data(testdat, param$frac[i])
    # o1 <- perform_imputation(td$s, td$ld, clump_gwas(td$s$Z, td$ld, param$zthresh[i], param$rthresh[i]))
    o1 <- eig_imp3(td$pc, param$eigthresh[i], td$s$BETA)
    m <- is.na(td$s$BETA)
    ev <- cor_out(o1$dat$X[m], td$s$BETA_TRUE[m])
    param$nout[i] <- ev$nout
    param$bcor[i] <- ev$c
    param$nmiss[i] <- sum(m)
  }
  return(param)
}


param <- expand.grid(
  sfile = a,
  method = "th",
  rthresh = c(0.01, 0.1),
  zthresh = qnorm(c(1e-4), low = F),
  frac = c(0.1)
) %>% mutate(sfile = as.character(sfile))

param <- subset(param, sfile == a[1])
testdat <- setup_data_for_tests(param$sfile[1])


res_th <- mclapply(a, \(A) {
  testdat <- try(setup_data_for_tests(A))
  if (class(testdat) == "try-error") {
    return(NULL)
  }
  pt <- run_test_th(testdat, subset(param, sfile == A))
  return(pt)
}, mc.cores = 50)

res_th <- lapply(res_th, \(x) {
  if (!inherits(x, "try-error")) {
    return(x)
  } else {
    return(NULL)
  }
}) %>% bind_rows()


param <- expand.grid(
  sfile = a,
  method = "eig",
  eigthresh = c(0.5, 0.6, 0.7, 0.8, 0.9, 0.99),
  frac = c(0.1, 0.5)
) %>% mutate(sfile = as.character(sfile))

res_eig <- mclapply(a, \(A) {
  testdat <- try(setup_data_for_tests(A))
  if (class(testdat) == "try-error") {
    return(NULL)
  }
  pt <- run_test_eig(testdat, subset(param, sfile == A))
  return(pt)
}, mc.cores = 50)

res_eig <- bind_rows(res_eig)

res <- bind_rows(res_eig, res_th) %>% as_tibble()


res_best <- res %>%
  arrange(desc(bcor)) %>%
  group_by(method, sfile, frac) %>%
  slice_head(n = 1)
res_best

table(res_best$eigthresh)
table(res_best$rthresh)
table(res_best$zthresh)

# Which is the best eig threshold

res_best2 <- res %>%
  arrange(desc(bcor)) %>%
  group_by(sfile, frac) %>%
  slice_head(n = 1)
res_best2
table(res_best2$method)



head(res_eig)

pt <- run_test_th(testdat, subset(param, sfile == a[1]))
pt2 <- run_test_eig(testdat, subset(param, sfile == a[1]))


pt
pt2


p1 <- ggplot(z, aes(a, b, col = abs(log(r)))) +
  geom_point()
ggsave(p1, file = "temp.png")


getr <- function(a, b) {
  ind <- a < 0
  a <- abs(a)
  b[ind] <- b[ind] * -1
  return(a / b)
}

z <- tibble(a = o1$dat$X[m], b = td$s$BETA_TRUE[m], r = getr(a, b))

z$r <- z$a / z$b2

z$cd <- cooks.distance(lm(a ~ b, z))


p1 <- ggplot(z, aes(a, b, col = cd)) +
  geom_point()
ggsave(p1, file = "temp.png")

z[outlier_detection(z$cd), ]


z %>% subset(b > 0.2)

mean(z$r)
