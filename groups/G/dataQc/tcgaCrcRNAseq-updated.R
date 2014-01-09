## FUNCTIONS TO EXTRACT DATA OBJECTS FROM SYNAPSE AND QC
#####
## ANALYST: BRIAN M. BOT
#####
require(synapseClient)
require(rGithubClient)
require(Biobase)
require(limma)
require(corpcor)
require(ggplot2)

## SOURCE IN HELPER FUNCTIONS FROM GITHUB
crcscRepo <- getRepo("/Sage-Bionetworks/crcsc")
code1 <- getPermlink(crcscRepo, "groups/G/dataQc/JGnorm.R")
sourceRepoFile(crcscRepo, "groups/G/dataQc/JGnorm.R")

## GET THE LOCATION OF THIS FILE ON GITHUB
rUrl <- getPermlink(crcscRepo, "groups/G/dataQc/tcgaCrcRNAseq-updated.R")

## READ IN SYNAPSE FILES AS FORMATED BY TCGA LIVE
loadTCGAExprFile <- function(f){
  as.matrix(read.table(getFileLocation(f), sep="\t", header=TRUE, as.is=TRUE, row.names=1, comment="", quote="", check.names=FALSE))
}

## TCGA PATIENT IDS PARSING
extractTcgaPatientIds <- function(tcgaIds){
  fixIds <- gsub("\\.","-", as.matrix(tcgaIds))
  patientIds <- sapply(strsplit(fixIds, "-", fixed=T), function(x){
    paste(x[1:3], collapse="-")
  })
  return(patientIds)
}
barcode2tumorType <- function(tcgaIds){
  fixIds <- gsub("\\.","-", as.matrix(tcgaIds))
  patientIds <- substr(sapply(strsplit(fixIds, "-", fixed=T), "[[", 4), 1, 2)
  return(patientIds)
}

## COMBINE PROBES TO GENES BY FIRST SV
combineProbesToGene <- function(expr, genes, method="svd"){
  
  if(is.list(genes)) genes <- unlist(genes)
  
  stopifnot(dim(expr)[1] ==  length(genes))
  ugenes <- unique(genes)
  ugenes <- sort(ugenes[!is.na(ugenes)])
  M <- matrix(NaN, ncol=dim(expr)[2], nrow=length(ugenes),
              dimnames=list(ugenes, colnames(expr)))
  
  for(gene in ugenes){
    subExpr <- as.matrix(expr[which(genes == gene),])
    if(dim(subExpr)[2] == 1){
      M[gene, ] <- subExpr
    }else{
      tmp <- svd(subExpr - rowMeans(subExpr))$v[,1]
      tmpC <- mean(cor(tmp, t(subExpr)))
      multiplier <- ifelse(tmpC < 0, -1, 1)
      M[gene,] <- tmp * multiplier
    }
  }
  return(M)
}

## CONVENIENCE FUNCTION FOR SVD EVALUATIONS
fs <- function(x){
  require(corpcor)
  u <- fast.svd(t(scale(t(x), scale = FALSE)), tol = 0)
  u$d <- u$d^2/sum(u$d^2)
  return(u)
}


## SYNAPSE FOLDER FOR THE TCGA DATA
synFolder <- "syn2023932"

synCoadhi <- synGet("syn2320092")
coadhi <- loadTCGAExprFile(synCoadhi)
coadhi <- coadhi[, barcode2tumorType(colnames(coadhi)) == "01"]
synCoadga <- synGet("syn2320079")
coadga <- loadTCGAExprFile(synCoadga)
coadga <- coadga[, barcode2tumorType(colnames(coadga)) == "01"]

synReadhi <- synGet("syn2320098")
readhi <- loadTCGAExprFile(synReadhi)
readhi <- readhi[, barcode2tumorType(colnames(readhi)) == "01"]
synReadga <- synGet("syn2320147")
readga <- loadTCGAExprFile(synReadga)
readga <- readga[, barcode2tumorType(colnames(readga)) == "01"]

hi <- log2(cbind(coadhi, readhi) + 1)
ga <- log2(cbind(coadga, readga) + 1)

if( all(rownames(hi) == rownames(ga)) ){
  theseFeatures <- rownames(hi)
} else{
  stop("rownames do not match")
}
## GET RID OF GENES WITH NO GENE SYMBOL
rns <- sapply(strsplit(rownames(hi), "|", fixed=T), "[[", 1)
idx <- rns != "?"
hi <- combineProbesToGene(hi[idx, ], rns[idx])
colnames(hi) <- extractTcgaPatientIds(colnames(hi))
ga <- combineProbesToGene(ga[idx, ], rns[idx])
colnames(ga) <- extractTcgaPatientIds(colnames(ga))


#####
## GA
#####
## WRITE OUT AN ACTIVITY THAT CAPTURES WHAT WAS USED IN OUR ANALYSIS
actGa <- Activity(name="RNA-seq QC", used=list(synCoadga, synReadga, list(url=code1, name=basename(code1), wasExecuted=FALSE), list(url=rUrl, name=basename(rUrl), wasExecuted=TRUE)))
actGa <- synStore(actGa)

## EXPRESSION FILE
ga <- as.data.frame(ga)
tmpNames <- colnames(ga)
ga$feature <- rownames(ga)
ga <- ga[, c("feature", tmpNames)]
tcgaCrcExprGaFile <- file.path(tempdir(), "TCGACRC_expression-ga.tsv")
write.table(ga, file=tcgaCrcExprGaFile, sep="\t", quote=FALSE, row.names=FALSE)

exprGaFile <- File(path=tcgaCrcExprGaFile, parentId=synFolder)
generatedBy(exprGaFile) <- actGa
exprGaFile <- synStore(exprGaFile)

#####
## HI
#####
## WRITE OUT AN ACTIVITY THAT CAPTURES WHAT WAS USED IN OUR ANALYSIS
actHi <- Activity(name="RNA-seq QC", used=list(synCoadhi, synReadhi, list(url=code1, name=basename(code1), wasExecuted=FALSE), list(url=rUrl, name=basename(rUrl), wasExecuted=TRUE)))
actHi <- synStore(actHi)

## EXPRESSION FILE
hi <- as.data.frame(hi)
tmpNames <- colnames(hi)
hi$feature <- rownames(hi)
hi <- hi[, c("feature", tmpNames)]
tcgaCrcExprHiFile <- file.path(tempdir(), "TCGACRC_expression-hi.tsv")
write.table(hi, file=tcgaCrcExprHiFile, sep="\t", quote=FALSE, row.names=FALSE)

exprHiFile <- File(path=tcgaCrcExprHiFile, parentId=synFolder)
generatedBy(exprHiFile) <- actHi
exprHiFile <- synStore(exprHiFile)

