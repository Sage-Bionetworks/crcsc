source("/shared/code/R/.Rprofile")

require(rGithubClient)
crcRepo <- getRepo("Sage-Bionetworks/crcsc")
sourceRepoFile(crcRepo, "evals/getDataFuncs.R")

groups <- paste("Group", c("A", "B", "C", "D", "E", "F", "G"), sep="")
dss <- names(allDatasets)

sgeFile <- "/shared/code/R/runThese.sh"
sgeText <- c("#!/bin/bash")
for(g in groups){
  for(d in dss){
    sgeCommand <- paste("qsub -V -wd /shared/code/repos/crcsc/evals -N ", g, "-", d, " -b y -o /shared/tmp/eoFiles -e /shared/tmp/eoFiles /usr/bin/Rscript evalGenesets.R ", d, " ", g, sep="")
    sgeText <- c(sgeText, sgeCommand)
  }
}

writeLines(sgeText, sgeFile)
