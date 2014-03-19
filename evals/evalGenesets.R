#####
## FOR (EVENTUALLY) RUNNING IN BATCH MODE ON AWS WITH ARGUMENTS DESCRIBED BELOW
#####
## SOURCE IN SHARED .Rprofile WHICH CONTAINS SYNAPSE LOGIN HOOK,
## SETS COMMON SYNAPSE CACHE FOR ALL WORKERS, AND SETS COMMON LIBPATH
source("/shared/code/R/.Rprofile")
#####
## TAKES FOR ARGUMENTS (PASSED FROM sgeKickoff.R)
#####
##   dataset: dataset to analyze
##   group:   group to pull results from
#####
myArgs <- commandArgs(trailingOnly=T)
ds <- myArgs[1]
group <- myArgs[2]

# ds <- "tcga_rnaseq"
# group <- "GroupG"

options(stringsAsFactors=F)

require(synapseClient)
require(rGithubClient)
require(affy)
require(limma)
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
sourceRepoFile(crcRepo, "evals/evalFuncs.R")
code3 <- getPermlink(crcRepo, "evals/evalFuncs.R")


## SOURCE CODE TO READ IN DATA
sourceRepoFile(crcRepo, "evals/getDataFuncs.R")
code4 <- getPermlink(crcRepo, "evals/getDataFuncs.R")

## THIS SCRIPT
thisCode <- getPermlink(crcRepo, "evals/evalGenesets.R")


#####
## GET ALL NECESSARY DATA TO RUN GENESET ANALYSIS FOR THIS GROUP AND DATASET
#####

## GET GROUP RESULTS FOR SPECIFIED GROUP
grpResId <- getGroupResultId(group, ds)
pmat <- getGroupResult(grpResId, group)
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
rownames(st) <- clean.names(rownames(st))
colnames(st) <- colnames(pmat)

## GET THE EXPRESSION DATA FOR THIS DATASET
## SUBSET TO AND ORDER LIKE THE SAMPLES IN THE SUBTYPE MATRIX
d <- getExprSet(ds)
sampleNames(d) <- clean.names(sampleNames(d))
d <- d[, as.character(rownames(st)) ]
d <- d[apply(exprs(d), 1, sd) != 0, ]

## GET THE GENESETS
genesets <- load.gmt.data(getFileLocation(synGet("syn2321865")))
genesets <- lapply(genesets, function(x){
  x <- x[ x != "" ]
  x <- unlist(symbolMap(x))
  x <- x[ !is.na(x) ]
  intersect(x, featureNames(d))
})


#####
## FIRST JUST RUN LMFIT ON EXPRESSION DATA FOR EACH SUBTYPE
#####
diffExprResults <- sapply(as.list(1:nSubtypes), function(i){
  resp <- st[, i]
  fit <- lmFit(d, design=model.matrix(~ factor(resp)))
  fit <- eBayes(fit)
})

diffExprPvalues <- sapply(diffExprResults, function(x){
  x$p.value[, "factor(resp)1"]
})
rownames(diffExprPvalues) <- featureNames(d)
colnames(diffExprPvalues) <- colnames(st)

diffExprFCs <- sapply(diffExprResults, function(x){
  2^x$coefficients[, "factor(resp)1"]
})
rownames(diffExprFCs) <- featureNames(d)
colnames(diffExprFCs) <- colnames(st)

pvalFile <- file.path(tempdir(), paste("diffExprPvalues-", group, "-", ds, ".tsv", sep=""))
write.table(diffExprPvalues, file=pvalFile, quote=F, sep="\t", col.names=NA)
pvalSyn <- synStore(File(path=pvalFile, parentId="syn2322802", group=group, dataset=ds, method="eBayes", stat="pvalue", evalDate=as.character(Sys.Date())), 
                    activity=Activity(name="differential expression",
                                      used=list(
                                        list(name=basename(code1), url=code1, wasExecuted=F),
                                        list(name=basename(code2), url=code2, wasExecuted=F),
                                        list(name=basename(code3), url=code3, wasExecuted=F),
                                        list(name=basename(code4), url=code4, wasExecuted=F),
                                        list(entity=synGet(allDatasets[[ds]]$exprSynId, downloadFile=F), wasExecuted=F),
                                        list(entity=synGet(grpResId, downloadFile=F), wasExecuted=F),
                                        list(name=basename(thisCode), url=thisCode, wasExecuted=T)
                                      )))

fcFile <- file.path(tempdir(), paste("diffExprFCs-", group, "-", ds, ".tsv", sep=""))
write.table(diffExprFCs, file=fcFile, quote=F, sep="\t", col.names=NA)
fcSyn <- synStore(File(path=fcFile, parentId="syn2322802", group=group, dataset=ds, method="eBayes", stat="fc", evalDate=as.character(Sys.Date())), 
                  activity=Activity(name="differential expression",
                                    used=list(
                                      list(name=basename(code1), url=code1, wasExecuted=F),
                                      list(name=basename(code2), url=code2, wasExecuted=F),
                                      list(name=basename(code3), url=code3, wasExecuted=F),
                                      list(name=basename(code4), url=code4, wasExecuted=F),
                                      list(entity=synGet(allDatasets[[ds]]$exprSynId, downloadFile=F), wasExecuted=F),
                                      list(entity=synGet(grpResId, downloadFile=F), wasExecuted=F),
                                      list(name=basename(thisCode), url=thisCode, wasExecuted=T)
                                    )))



#####
## RUN GENESET EVALUATION
#####

## GLOBAL TEST (NON COMPETITIVE)
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

gtFile <- file.path(tempdir(), paste("gt-", group, "-", ds, ".tsv", sep=""))
write.table(gtResults, file=gtFile, quote=F, sep="\t", col.names=NA)
gtSyn <- synStore(File(path=gtFile, parentId="syn2322802", group=group, dataset=ds, method="globaltest", evalDate=as.character(Sys.Date())), 
                  activity=Activity(name="geneset evaluation",
                                    used=list(
                                      list(name=basename(code1), url=code1, wasExecuted=F),
                                      list(name=basename(code2), url=code2, wasExecuted=F),
                                      list(name=basename(code3), url=code3, wasExecuted=F),
                                      list(name=basename(code4), url=code4, wasExecuted=F),
                                      list(entity=synGet(allDatasets[[ds]]$exprSynId, downloadFile=F), wasExecuted=F),
                                      list(entity=synGet("syn2321865", downloadFile=F), wasExecuted=F),
                                      list(entity=synGet(grpResId, downloadFile=F), wasExecuted=F),
                                      list(name=basename(thisCode), url=thisCode, wasExecuted=T)
                                    )))


## GSA
## RESULTS AVAILABLE FOR BOTH HI AND LO
gsaResults <- lapply(as.list(1:nSubtypes), function(i){
  ## GSA REQUIRES 1 AND 2 INSTEAD OF 0 AND 1
  resp <- st[, i] + 1
  op <- GSA(x=exprs(d), y=resp, genesets=genesets, genenames=featureNames(d), resp.type="Two class unpaired", nperms=10000, minsize=3)
  op
})

gsaHiResults <- sapply(gsaResults, function(r){
  r$pvalues.hi
})
rownames(gsaHiResults) <- names(genesets)
colnames(gsaHiResults) <- colnames(st)

# gsaLoResults <- sapply(gsaResults, function(r){
#   r$pvalues.lo
# })
# rownames(gsaLoResults) <- names(genesets)
# colnames(gsaLoResults) <- colnames(st)

gsaFile <- file.path(tempdir(), paste("gsa-", group, "-", ds, ".tsv", sep=""))
write.table(gsaHiResults, file=gsaFile, quote=F, sep="\t", col.names=NA)
gsaSyn <- synStore(File(path=gsaFile, parentId="syn2322802", group=group, dataset=ds, method="gsa", evalDate=as.character(Sys.Date())), 
                   activity=Activity(name="geneset evaluation",
                                     used=list(
                                       list(name=basename(code1), url=code1, wasExecuted=F),
                                       list(name=basename(code2), url=code2, wasExecuted=F),
                                       list(name=basename(code3), url=code3, wasExecuted=F),
                                       list(name=basename(code4), url=code4, wasExecuted=F),
                                       list(entity=synGet(allDatasets[[ds]]$exprSynId, downloadFile=F), wasExecuted=F),
                                       list(entity=synGet("syn2321865", downloadFile=F), wasExecuted=F),
                                       list(entity=synGet(grpResId, downloadFile=F), wasExecuted=F),
                                       list(name=basename(thisCode), url=thisCode, wasExecuted=T)
                                     )))



## KS TEST
## DIRECTION DOES NOT MATTER - ONLY DIFFERING SIGNIFICANCE DISTRIBUTIONS
ksResults <- sapply(as.list(1:nSubtypes), function(i){
  op <- sapply(genesets, function(gs){
    ks.test(x=diffExprPvalues[which(rownames(diffExprPvalues) %in% gs), i],
            y=diffExprPvalues[-which(rownames(diffExprPvalues) %in% gs), i],
            alternative="less")$p.value
  })
})
colnames(ksResults) <- colnames(st)

ksFile <- file.path(tempdir(), paste("ks-", group, "-", ds, ".tsv", sep=""))
write.table(ksResults, file=ksFile, quote=F, sep="\t", col.names=NA)
ksSyn <- synStore(File(path=ksFile, parentId="syn2322802", group=group, dataset=ds, method="ks", evalDate=as.character(Sys.Date())), 
                  activity=Activity(name="geneset evaluation",
                                    used=list(
                                      list(name=basename(code1), url=code1, wasExecuted=F),
                                      list(name=basename(code2), url=code2, wasExecuted=F),
                                      list(name=basename(code3), url=code3, wasExecuted=F),
                                      list(name=basename(code4), url=code4, wasExecuted=F),
                                      list(entity=synGet(allDatasets[[ds]]$exprSynId, downloadFile=F), wasExecuted=F),
                                      list(entity=synGet("syn2321865", downloadFile=F), wasExecuted=F),
                                      list(entity=synGet(grpResId, downloadFile=F), wasExecuted=F),
                                      list(name=basename(thisCode), url=thisCode, wasExecuted=T)
                                    )))

diffExprResults <- sapply(as.list(1:nSubtypes), function(i){
  resp <- st[, i]
  fit <- lmFit(d, design=model.matrix(~ factor(resp)))
  fit <- eBayes(fit)
})

diffExprPvalues <- sapply(diffExprResults, function(x){
  x$p.value[, "factor(resp)1"]
})
rownames(diffExprPvalues) <- featureNames(d)
colnames(diffExprPvalues) <- colnames(st)

## TUKEY NON COMPETITIVE TEST
set.seed(20140101)
tukResults <- sapply(as.list(1:nSubtypes), function(i){
  resp <- st[, i]
  
  op <- sapply(genesets, function(gs){
    theseGenes <- d[gs, ]
    fit <- lmFit(theseGenes, design=model.matrix(~ factor(resp)))
    fit <- eBayes(fit)
    
    mgd <- fit$p.value[, "factor(resp)1"]
    mgd <- sum(mgd<0.05)
    
    perms <- numeric()
    for(j in 1:10000){
      ran <- rnorm(length(resp))
      resp2 <- resp[order(ran)]
      fit2 <- lmFit(theseGenes, design=model.matrix(~ factor(resp2)))
      fit2 <- eBayes(fit2)
      a <- fit2$p.value[, "factor(resp2)1"]
      
      perms <- c(perms, sum(a<0.05))
      if(j/25 == floor(j/25)){
        cat("permutation ", j, "\n")
      }
    }
    pval <- sum(perms>mgd)/length(perms)
    pval
  })
  names(op) <- names(genesets)
  op
})
colnames(tukResults) <- colnames(st)

tukFile <- file.path(tempdir(), paste("tuk-", group, "-", ds, ".tsv", sep=""))
write.table(tukResults, file=tukFile, quote=F, sep="\t", col.names=NA)
tukSyn <- synStore(File(path=tukFile, parentId="syn2322802", group=group, dataset=ds, method="tukey", evalDate=as.character(Sys.Date())), 
                   activity=Activity(name="geneset evaluation",
                                     used=list(
                                       list(name=basename(code1), url=code1, wasExecuted=F),
                                       list(name=basename(code2), url=code2, wasExecuted=F),
                                       list(name=basename(code3), url=code3, wasExecuted=F),
                                       list(name=basename(code4), url=code4, wasExecuted=F),
                                       list(entity=synGet(allDatasets[[ds]]$exprSynId, downloadFile=F), wasExecuted=F),
                                       list(entity=synGet("syn2321865", downloadFile=F), wasExecuted=F),
                                       list(entity=synGet(grpResId, downloadFile=F), wasExecuted=F),
                                       list(name=basename(thisCode), url=thisCode, wasExecuted=T)
                                     )))


