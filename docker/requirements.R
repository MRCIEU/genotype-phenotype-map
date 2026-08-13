options(repos = c(CRAN = "https://cloud.r-project.org"))
install.packages("BiocManager")
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  stop("FATAL ERROR: BiocManager failed to install.")
}

install.packages(c("remotes", "devtools"))

cran_install <- c(
  "testthat", "susieR", "Rfast", "duckdb", "validate", "redux", "sendmailR", "igraph", "svglite", "readxl",
  "janitor", "styler", "binom", "epitools", "scales", "ggrepel", "pheatmap", "UpSetR", "ieugwasr"
)
install.packages(cran_install)

# Pin lintr to match CI (.github/workflows/main.yml) and local make lint behaviour.
remotes::install_version("lintr", version = "3.3.0.1")

# Install packages from GitHub and specific versions
remotes::install_version("RcppEigen", version = "0.3.3.9.3")
devtools::install_github("jrs95/hyprcoloc", upgrade = "never")
devtools::install_github("MRCIEU/gwasglue", upgrade = "never")
devtools::install_github("MRCIEU/gpmapr", upgrade = "never")
devtools::install_github("MRCIEU/TwoSampleMR", upgrade = "never")
devtools::install_github("WSpiller/MVMR", upgrade = "never")

biocmanager_install <- c("Homo.sapiens", "GenomicRanges", "biomaRt")
BiocManager::install(biocmanager_install)
