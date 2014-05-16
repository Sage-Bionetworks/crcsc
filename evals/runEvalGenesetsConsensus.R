source("/shared/code/R/.Rprofile")

require(synapseClient)

## GET CONSENSUS RESULTS
grpResId <- "syn2469968"
c <- synGet(grpResId)
cms <- read.csv(getFileLocation(c), as.is=T)
d <- sapply(strsplit(cms$dataset.sample, ".", fixed=T), "[", 1)
cms$dataset <- d
cms <- cms[cms$dataset != "tcga_rnaseq", ]

dss <- names(table(cms$dataset))

sgeFile <- "/shared/code/R/runTheseConsensus.sh"
sgeText <- c("#!/bin/bash")
for(d in dss){
  sgeCommand <- paste("qsub -V -wd /shared/code/repos/crcsc/evals -N ", g, "-", d, " -b y -o /shared/tmp/eoFiles -e /shared/tmp/eoFiles /usr/bin/Rscript evalGenesetsConsensus.R ", d, sep="")
  sgeText <- c(sgeText, sgeCommand)
}

writeLines(sgeText, sgeFile)
