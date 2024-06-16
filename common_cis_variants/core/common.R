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
