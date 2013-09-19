require(rGithubClient)
require(ggplot2)
require(synapseClient)
synapseLogin()

options(stringsAsFactors=F)

## GET INFO ABOUT THE GITHUB REPOSITORY
crcscRepo <- getRepo("Sage-Bionetworks/crcsc")
thisScript <- getPermlink(crcscRepo, "bootstrap/exampleTcgaAnalysis.R")

## SPECIFY SYNAPSE IDS
folderId <- "syn2231831"
exprId <- "syn2161141"
clinId <- "syn2165691"

## FUNCTION TO READ IN FILES AND SET ROWNAMES AND COLUMN NAMES
readCrcscTables <- function(synId, num=F){
  tmp <- synGet(synId)
  t <- read.delim(getFileLocation(tmp), as.is=T, header=F)
  rns <- t[-1, 1]
  cns <- t[1, -1]
  t <- t[-1, -1]
  if(num){
    t <- apply(t, 2, as.numeric)
  }
  rownames(t) <- rns
  colnames(t) <- cns
  return(t)
}

expr <- readCrcscTables(exprId, num=T)
clin <- readCrcscTables(clinId)

genderPvals  <- apply(expr, 1, function(x){
  fit  <- summary(lm(x ~ clin$gender))
  return(fit$coefficients[2, 4])
})

plotPath  <- file.path(tempdir(), "pvalueHistogram.png")
png(plotPath)
hist(genderPvals, main="", xlab="TCGA P-value histogram for gender")
dev.off()

plotFile  <- File(path=plotPath, parentId=folderId)
plotFile  <- synStore(plotFile, used=list(list(entity=exprId, wasExecuted=F),
                                          list(entity=clinId, wasExecuted=F),
                                          list(url=thisScript, name=basename(thisScript), wasExecuted=T)))


