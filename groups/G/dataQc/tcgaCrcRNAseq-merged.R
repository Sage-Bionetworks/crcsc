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
rUrl <- getPermlink(crcscRepo, "groups/G/dataQc/tcgaCrcRNAseq-merged.R")

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

## IF SAMPLE RUN ON BOTH, KEEP HI-SEQ SAMPLE
coadga <- coadga[, !(colnames(coadga) %in% colnames(coadhi))]
readga <- readga[, !(colnames(readga) %in% colnames(readhi))]

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
ga <- combineProbesToGene(ga[idx, ], rns[idx])

## MERGE TOGETHER
expr <- cbind(hi,ga)

thesePatients <- extractTcgaPatientIds(colnames(expr))

if( all(duplicated(thesePatients) == FALSE) ){
  colnames(expr) <- thesePatients
  colnames(hi) <- extractTcgaPatientIds(colnames(hi))
  colnames(ga) <- extractTcgaPatientIds(colnames(ga))
} else{
  stop("duplicated patients")
}


## SVD ON EXPRESSION MATRIX -- ASSESS OVERALL STRUCTURE AND POSSIBLE LATENT STRUCTURE
s <- fs(expr)
platform <- c(rep("hiseq", ncol(hi)), rep("ga", ncol(ga)))

qplot(1:length(s$d), s$d,
      xlab="eigen gene",
      ylab="% variance explained")
qplot(s$v[, 1], s$v[, 2], colour=platform,
      xlab="1st svd",
      ylab="2nd svd")


# adjust using Combat
v1 <- apply(hi, 1, var)
v2 <- apply(ga, 1, var)
mask <- !(v1 == 0 & v2 == 0)
fit <- eb(hi[mask,],ga[mask,])
hi.adj <- fit$x
ga.adj <- fit$y

# recheck PCs
exprAdj <- cbind(hi.adj,ga.adj)

## SVD ON EXPRESSION MATRIX -- ASSESS OVERALL STRUCTURE AND POSSIBLE LATENT STRUCTURE
sAdj <- fs(exprAdj)

qplot(1:length(sAdj$d), sAdj$d,
      xlab="eigen gene",
      ylab="% variance explained")
qplot(sAdj$v[, 1], sAdj$v[, 2], colour=platform,
      xlab="1st svd",
      ylab="2nd svd")

## WRITE OUT AN ACTIVITY THAT CAPTURES WHAT WAS USED IN OUR ANALYSIS
act <- Activity(name="RNA-seq QC", used=list(synCoadhi, synCoadga, synReadhi, synReadga, list(url=code1, name=basename(code1), wasExecuted=FALSE), list(url=rUrl, name=basename(rUrl), wasExecuted=TRUE)))
act <- synStore(act)

## EXPRESSION FILE
exprAdj <- as.data.frame(exprAdj)
tmpNames <- colnames(exprAdj)
exprAdj$feature <- rownames(exprAdj)
exprAdj <- expr[, c("feature", tmpNames)]
tcgaCrcExprFile <- file.path(tempdir(), "TCGACRC_expression-merged.tsv")
write.table(expr, file=tcgaCrcExprFile, sep="\t", quote=FALSE, row.names=FALSE)

exprFile <- File(path=tcgaCrcExprFile, parentId=synFolder)
generatedBy(exprFile) <- act
exprFile <- synStore(exprFile)

