#cran_install <- c("susieR", "Rfast")
# biocmanager_install <- c("Homo.sapiens", "GenomicRanges")
devtools::install_github("ZikunY/CARMA")
remotes::install_version("RcppEigen", version = "0.3.3.9.3")
devtools::install_github("jrs95/hyprcoloc", upgrade = "never")

#install.packages(cran_install, repos = "http://cran.us.r-project.org")
#BiocManager::install(biocmanager_install)
