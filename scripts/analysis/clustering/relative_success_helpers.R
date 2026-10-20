# Shared statistics for the drug target-indication (T-I) pair analyses
# (02_gpmap_support_for_ti_pairs.Rmd, 03_ti_pair_validation_of_coloc.Rmd; paper Figure 6).
#
# Unit of analysis: one T-I pair. These are observational records of drug programmes,
# so there are no biological/technical replicates.

# Relative success (RS) from a 2x2 contingency table built as
#   table(<genetic support: FALSE/TRUE>, <success: FALSE/TRUE>)
# i.e. row 1 = comparison group (no genetic support), row 2 = supported; column 2 = success.
#
# RS = P(success | supported) / P(success | not supported): the risk ratio from
# epitools::riskratio(method = "wald"), with a Wald 95% CI on the log scale.
# Also returns the counts for both groups (for the figure legends) and two-sided exact
# P value for the 2x2 table (Fisher's exact test).
# No adjustment for multiple comparisons is applied.
# Returns a one-row data frame of NAs when the table is not 2x2.
rs_from_table <- function(tab) {
  empty <- data.frame(
    n_supp_launched = NA_integer_, n_supp = NA_integer_,
    n_unsupp_launched = NA_integer_, n_unsupp = NA_integer_,
    RS = NA_real_, CI_lower = NA_real_, CI_upper = NA_real_,
    p_fisher = NA_real_
  )
  if (any(dim(tab) != 2)) {
    return(empty)
  }
  rr <- epitools::riskratio(tab, method = "wald")
  est <- rr$measure[2, ]
  data.frame(
    n_supp_launched = as.integer(tab[2, 2]),
    n_supp = as.integer(sum(tab[2, ])),
    n_unsupp_launched = as.integer(tab[1, 2]),
    n_unsupp = as.integer(sum(tab[1, ])),
    RS = unname(est["estimate"]),
    CI_lower = unname(est["lower"]),
    CI_upper = unname(est["upper"]),
    p_fisher = unname(rr$p.value[2, "fisher.exact"])
  )
}

# Write the numbers behind a Figure 6 panel (legend n values / Source Data) to
# <analysis_dir>/<name><trans_suffix>.tsv
write_panel_table <- function(df, name, analysis_dir, trans_suffix = "") {
  dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)
  out_file <- file.path(analysis_dir, paste0(name, trans_suffix, ".tsv"))
  data.table::fwrite(as.data.frame(df), out_file, sep = "\t")
  message("Wrote ", out_file)
  invisible(out_file)
}
