cran_install <- c('susieR', 'Rfast', 'BiocManager', 'duckdb', 'validate', 'redux', 'sendmailR', 'igraph', 'svglite', 'readxl', 'janitor')
install.packages(cran_install, repos = 'http://cran.us.r-project.org')

remotes::install_version('RcppEigen', version = '0.3.3.9.3')
devtools::install_github('jrs95/hyprcoloc', upgrade = 'never')
devtools::install_github('MRCIEU/gwasglue', upgrade = 'never')

biocmanager_install <- c('Homo.sapiens', 'GenomicRanges', 'biomaRt')
BiocManager::install(biocmanager_install)
