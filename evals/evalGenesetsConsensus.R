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
#####
myArgs <- commandArgs(trailingOnly=T)
ds <- myArgs[1]
# ds <- "tcga_rnaseqAll"
group <- "cms4"


options(stringsAsFactors=F)

require(synapseClient)
require(rGithubClient)
require(affy)
require(limma)
require(hgu133plus2.db)
require(hgu133a2.db)
require(org.Hs.eg.db)

## GENE SET METHODS TO BE USED
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
thisCode <- getPermlink(crcRepo, "evals/evalGenesetsConsensus.R")


#####
## GET ALL NECESSARY DATA TO RUN GENESET ANALYSIS FOR THIS GROUP AND DATASET
#####

## GET CONSENSUS RESULTS
grpResId <- "syn2469968"
c <- synGet("grpResId")
cms <- read.csv(getFileLocation(c), as.is=T)
d <- sapply(strsplit(cms$dataset.sample, ".", fixed=T), "[", 1)
cms$dataset <- d
cms <- cms[cms$dataset != "tcga_rnaseq", ]
samp <- sapply(strsplit(cms$dataset.sample, ".", fixed=T), "[", 2)
cms$sample <- samp
rownames(cms) <- samp

cms <- cms[ cms$dataset == ds, ]

theseCfs <- names(table(cms$cms4))
tmp <- lapply(as.list(theseCfs), function(x){
  as.numeric(cms$cms4 == x)
})
st <- do.call(cbind, tmp)
colnames(st) <- theseCfs
rownames(st) <- rownames(cms)
nSubtypes <- ncol(st)

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
pvalSyn <- synStore(File(path=pvalFile, parentId="syn2476109", group=group, dataset=ds, method="eBayes", stat="pvalue", evalDate=as.character(Sys.Date())), 
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
fcSyn <- synStore(File(path=fcFile, parentId="syn2476109", group=group, dataset=ds, method="eBayes", stat="fc", evalDate=as.character(Sys.Date())), 
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

gsaFile <- file.path(tempdir(), paste("gsa-", group, "-", ds, ".tsv", sep=""))
write.table(gsaHiResults, file=gsaFile, quote=F, sep="\t", col.names=NA)
gsaSyn <- synStore(File(path=gsaFile, parentId="syn2476109", group=group, dataset=ds, method="gsa", evalDate=as.character(Sys.Date())), 
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

