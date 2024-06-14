cran_install <- c("susieR", "Rfast")
install.packages(cran_install, repos = "http://cran.us.r-project.org")

devtools::install_github("ZikunY/CARMA")
remotes::install_version("RcppEigen", version = "0.3.3.9.3")
devtools::install_github("jrs95/hyprcoloc", upgrade = "never")

#biocmanager_install <- c("Homo.sapiens", "GenomicRanges")
#BiocManager::install(biocmanager_install)
