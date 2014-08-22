source("/shared/code/R/.Rprofile")

require(synapseClient)

## GET CONSENSUS RESULTS
grpResId <- "syn2533623"
c <- synGet(grpResId)
cms <- read.csv(getFileLocation(c), as.is=T)
d <- sapply(strsplit(cms$NewCMS4_unclear, ".", fixed=T), "[", 1)
cms$dataset <- d
cms$dataset[ cms$dataset == "tcgacrc_merged" ] <- "tcga_rnaseqAll"
cms$dataset[ cms$dataset == "petacc3" ] <- "petacc"

dss <- names(table(cms$dataset))

sgeFile <- "/shared/code/R/runTheseConsensus.sh"
sgeText <- c("#!/bin/bash")
for(d in dss){
  sgeCommand <- paste("qsub -V -wd /shared/code/repos/crcsc/evals -N ", d, " -b y -o /shared/tmp/eoFiles -e /shared/tmp/eoFiles /usr/bin/Rscript evalGenesetsConsensus.R ", d, sep="")
  sgeText <- c(sgeText, sgeCommand)
}

writeLines(sgeText, sgeFile)
