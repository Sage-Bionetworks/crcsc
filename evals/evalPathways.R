#####
## FOR (EVENTUALLY) RUNNING IN BATCH MODE ON AWS WITH ARGUMENTS DESCRIBED BELOW
#####
## SOURCE IN SHARED .Rprofile WHICH CONTAINS SYNAPSE LOGIN HOOK,
## SETS COMMON SYNAPSE CACHE FOR ALL WORKERS, AND SETS COMMON LIBPATH
# source("/shared/code/R/.Rprofile")
#####
## TAKES FOR ARGUMENTS (PASSED FROM sgeKickoff.R)
#####
##   dataset: dataset to analyze
##   group:   group to pull results from
#####
# myArgs <- commandArgs(trailingOnly=T)
# ds <- myArgs[1]
# group <- myArgs[2]

ds <- "tcga_rnaseq"
group <- "GroupG"

options(stringsAsFactors=F)

require(synapseClient)
require(rGithubClient)
require(affy)
require(pamr)
require(hgu133plus2.db)
require(hgu133a2.db)
require(org.Hs.eg.db)

## GENE SET METHODS TO BE USED
require(globaltest)
require(GSA)

# password will be request after calling this
# synapseLogin()

## SOURCE IN BACKGROUND FUNCTIONS FROM JG
crcRepo <- getRepo("Sage-Bionetworks/crcsc")
sourceRepoFile(crcRepo, "groups/G/pipeline/JGLibrary.R")
code1 <- getPermlink(crcRepo, "groups/G/pipeline/JGLibrary.R")
sourceRepoFile(crcRepo, "groups/G/pipeline/subtypePipelineFuncs.R")
code2 <- getPermlink(crcRepo, "groups/G/pipeline/subtypePipelineFuncs.R")

## SOURCE CODE TO READ IN DATA
sourceRepoFile(crcRepo, "evals/getDataFuncs.R")
code3 <- getPermlink(crcRepo, "evals/getExprPhenoData.R")




#####
## GET ALL NECESSARY DATA TO RUN GENESET ANALYSIS FOR THIS GROUP AND DATASET
#####

## GET GROUP RESULTS FOR SPECIFIED GROUP
groupResults <- getGroupResults(group)

if( !any(names(groupResults) == ds) ){
  stop(paste("Group ", group, " did not provide results for ", dataset, sep=""))
}

pmat <- groupResults[[ds]]
nSubtypes <- ncol(pmat)
st <- apply(pmat, 1, function(x){
  ww <- x==max(x)
  if( sum(x==max(x)) == 1 ){
    ## IF THERE IS ONLY ONE MAX - RETURN THAT SUBTYPE
    as.numeric(ww)
  } else{
    ## IF THERE ARE MORE THAN ONE - RETURN RANDOMLY SELECTED SUBTYPE
    r <- rnorm(sum(ww))
    rr <- 1:length(ww) == (which(ww)[r==max(r)])
    ww*rr
  }
})
st <- t(st)
if( ds == "tcga_rnaseq" ){
  rownames(st) <- gsub(".", "-", rownames(st), fixed=T)
}
colnames(st) <- colnames(pmat)

## GET THE EXPRESSION DATA FOR THIS DATASET
## SUBSET TO AND ORDER LIKE THE SAMPLES IN THE SUBTYPE MATRIX
d <- getExprSet(ds)
d <- d[, rownames(st)]
d <- d[apply(exprs(d), 1, sd) != 0, ]

## GET THE GENESETS
genesets <- read.delim(getFileLocation(synGet("syn2319124")), as.is=T, header=F, row.names=1)
genesets <- apply(genesets, 1, function(x){
  x <- x[-1]
  x <- x[x != ""]
  names(x) <- NULL
  x <- unlist(symbolMap(x))
  x <- x[!is.na(x)]
  intersect(x, featureNames(d))
})





#####
## RUN GENESET EVALUATION
#####
# respList <- list()

## GLOBAL TEST
gtResults <- sapply(as.list(1:nSubtypes), function(i){
  resp <- st[, i]
  op <- sapply(genesets, function(gs){
    gtModel <- gt(response=factor(resp), alternative=d[ gs, ])
    gtModel@result[1, "p-value"]
  })
  names(op) <- names(genesets)
  op
})
colnames(gtResults) <- colnames(st)

## TUKEY NON COMPETITIVE TEST
set.seed(20140101)
tukResults <- sapply(as.list(1:nSubtypes), function(i){
  resp <- st[, i]
  
  op <- sapply(genesets, function(gs){
    theseGenes <- exprs(d)[gs, ]
    mgd <- apply(theseGenes, 1, function(y){
      summary(lm(y~resp))$coefficients[2, "Pr(>|t|)"]
    })
    mgd <- sum(mgd<0.05)
    
    perms <- numeric()
    for(j in 1:1000){
      ran <- rnorm(length(resp))
      a <- apply(theseGenes, 1, function(y){
        summary(lm(y~resp[order(ran)]))$coefficients[2, "Pr(>|t|)"]
      })
      perms <- c(perms, sum(a<0.05))
      if(j/25 == floor(j/25)){
        cat(gs, " - permutation ", j, "\n")
      }
    }
    pval <- sum(perms>mgd)/length(perms)
    pval
  })
  names(op) <- names(genesets)
  op
})
colnames(tukResults) <- colnames(st)



## GSA
## RESULTS AVAILABLE FOR BOTH HI AND LO
gsaResults <- lapply(as.list(1:nSubtypes), function(i){
  ## GSA REQUIRES 1 AND 2 INSTEAD OF 0 AND 1
  resp <- st[, i] + 1
  op <- GSA(x=exprs(d), y=resp, genesets=genesets, genenames=featureNames(d), resp.type="Two class unpaired", nperms=1000, minsize=3)
  op
})

gsaHiResults <- sapply(gsaResults, function(r){
  r$pvalues.hi
})
rownames(gsaHiResults) <- names(genesets)
colnames(gsaHiResults) <- colnames(st)

gsaLoResults <- sapply(gsaResults, function(r){
  r$pvalues.lo
})
rownames(gsaLoResults) <- names(genesets)
colnames(gsaLoResults) <- colnames(st)



## KS TEST
## DIRECTION DOES NOT MATTER - ONLY DIFFERING SIGNIFICANCE DISTRIBUTIONS
diffExpResults <- lapply(as.list(1:nSubtypes), function(i){
  resp <- st[, i]
  op <- apply(exprs(d), 1, function(y){
    summary(lm(y~resp))$coefficients[2, c("Estimate", "Pr(>|t|)")]
  })
  op
})
diffExpFCs <- sapply(diffExpResults, function(x){
  2^as.numeric(x[1, ])
})
rownames(diffExpFCs) <- featureNames(d)
colnames(diffExpFCs) <- colnames(st)

diffExpPvalues <- sapply(diffExpResults, function(x){
  as.numeric(x[2, ])
})
rownames(diffExpPvalues) <- featureNames(d)
colnames(diffExpPvalues) <- colnames(st)


ksResults <- sapply(as.list(1:nSubtypes), function(i){
  op <- sapply(genesets, function(gs){
    ks.test(x=diffExpPvalues[which(rownames(diffExpPvalues) %in% gs), i],
            y=diffExpPvalues[-which(rownames(diffExpPvalues) %in% gs), i],
            alternative="less")$p.value
  })
})
colnames(ksResults) <- colnames(st)


