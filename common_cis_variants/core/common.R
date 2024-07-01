#TODO: look through these, cause I think most are not being used...  Delete later


#From Gib's handy Z -> beta/se labbook: https://explodecomputer.github.io/lab-book/posts/2024-06-25-z-to-beta/
# Setup data to have Z scores and some known betas
setup_data <- function(a) {
  a$z <- a$beta / a$se
  a$beta_known <- NA
  a$se_known <- NA
  index <- sample(1:nrow(a), nrow(a)*0.05)
  a$beta_known[index] <- a$beta[index]
  a$se_known[index] <- a$se[index]
  return(a)
}

make_beta <- function(dat) {
  dat$senew <- 1 / sqrt(2 * dat$eaf * (1 - dat$eaf) * dat$n)
  dat$betanew <- dat$z * dat$senew
  correction <- lm(dat$betanew ~ dat$beta_known)$coef[2]
  dat$betanew <- dat$betanew / correction
  dat$senew <- dat$senew / correction
  return(dat)
}

standardise_gwas <- function(gwas, N=0, input_gwas_columns=list(), remove_extra_columns=F) {
  gwas <- gwas |>
    change_column_names(input_gwas_columns, remove_extra_columns) |>
    standardise_columns(N) |>
    filter_incomplete_rows()

  return(gwas)
}
change_column_names <- function(gwas, columns = list(), remove_extra_columns = F) {
  for (name in names(columns)) {
    names(gwas)[names(gwas) == columns[name]] <- name
  }

  if (remove_extra_columns) {
    columns_to_remove <- setdiff(colnames(gwas), names(columns))
    gwas <- gwas[,-which(colnames(gwas) %in% columns_to_remove)]
  }

  return(gwas)
}

convert_beta_and_se_to_z_score <- function(beta, se) {
  return(beta / se)
}

convert_negative_log_p_to_p <- function(lp) {
  return(10^-lp)
}

load_ld_region <- function(ld_region_file) {
  ld_region <- vroom::vroom(ld_region_file, col_names=F)
  colnames(ld_region) <- ld_region$X1
  ld_region <- ld_region[, 2:(ncol(ld_region)-1)]
  return(ld_region)
}

trim_ld_region_for_study <- function(ld_region, rsids) {
  keep <- colnames(ld_region) %in% gwas$rsid
  ld_for_gwas <- ld_region[keep, keep]
  ld_matrix <- matrix(as.vector(data.matrix(ld_for_gwas)), nrow=nrow(gwas), ncol=nrow(gwas))

  ld_region <- vroom::vroom(ld_region_file, col_names=F)
  colnames(ld_region) <- ld_region$X1
  ld_region <- ld_region[, 2:(ncol(ld_region)-1)]
  return(ld_region)
}
