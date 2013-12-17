source("/shared/code/R/.Rprofile")

install.packages(c("GSA", "devtools", "data.table"))

source("http://bioconductor.org/biocLite.R")
biocLite(c("globaltest", "affy", "limma", "hgu133plus2.db", "hgu133a2.db", "org.Hs.eg.db", "impute"))

source("http://depot.sagebase.org/CRAN.R")
pkgInstall("synapseClient")

require(devtools)
install_github("rGithubClient", "brian-bot", ref="rGithubClient-0.8")
