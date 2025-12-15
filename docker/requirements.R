options(repos = c(CRAN = "https://cloud.r-project.org"))
install.packages('BiocManager')
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  stop("FATAL ERROR: BiocManager failed to install.")
}

install.packages(c('remotes', 'devtools'))

cran_install <- c('testthat', 'susieR', 'Rfast', 'duckdb', 'validate', 'redux', 'sendmailR', 'igraph', 'svglite', 'readxl', 'janitor')
install.packages(cran_install)

# Install packages from GitHub and specific versions
remotes::install_version('RcppEigen', version = '0.3.3.9.3')
devtools::install_github('jrs95/hyprcoloc', upgrade = 'never')
devtools::install_github('MRCIEU/gwasglue', upgrade = 'never')

biocmanager_install <- c('Homo.sapiens', 'GenomicRanges', 'biomaRt')
BiocManager::install(biocmanager_install)
