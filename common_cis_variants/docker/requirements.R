cran_install <- c("BiocManager")
biocmanager_install <- c("GenomicRanges")
github_install <- c("jrs95/hyprcoloc")

installed_correctly <- function(packages) {
  for (package in packages) {
    if (!library(package, character.only=TRUE, logical.return=TRUE) ) {
      quit(status=1, save='no')
    }
  }
}

install.packages(cran_install, repos = "http://cran.us.r-project.org")
installed_correctly(cran_install)

BiocManager::install(biocmanager_install)
installed_correctly(biocmanager_install)

#remotes::install_version("RcppEigen", version = "0.3.3.9.3")
#remotes::install_github(github_install)
#installed_correctly(github_install)
