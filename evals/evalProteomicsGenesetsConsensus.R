options(stringsAsFactors=F)

require(synapseClient)
require(rGithubClient)
require(limma)
require(hgu133plus2.db)
require(hgu133a2.db)
require(org.Hs.eg.db)

## GENE SET METHODS TO BE USED
require(GSA)

ds <- "tcga_rnaseqAll"
group <- "cms4"


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
thisCode <- getPermlink(crcRepo, "evals/evalProteomicsGenesetsConsensus.R")


#####
## GET ALL NECESSARY DATA TO RUN GENESET ANALYSIS FOR THIS GROUP AND DATASET
#####

## GET CONSENSUS RESULTS
grpResId <- "syn2755232"
c <- synGet(grpResId)
cms <- read.csv(getFileLocation(c), as.is=T)
## REMOVE SAMPLES WITHOUT A CLEAR SUBTYPE
cms <- cms[ !(cms$CMS4network_plus_classifier_in_noncore_samples %in% c("UNK", "NOLBL")), ]
cms <- cms[ cms$dataset == ds, ]

theseCfs <- names(table(cms$CMS4network_plus_classifier_in_noncore_samples))
tmp <- lapply(as.list(theseCfs), function(x){
  as.numeric(cms$CMS4network_plus_classifier_in_noncore_samples == x)
})
st <- do.call(cbind, tmp)
colnames(st) <- theseCfs
rownames(st) <- cms$sample
nSubtypes <- ncol(st)

## GET THE PROTEOMICS DATA
## SUBSET TO AND ORDER LIKE THE SAMPLES IN THE SUBTYPE MATRIX
dd <- synGet("syn2767938")
d <- read.csv(getFileLocation(dd), as.is=T, row.names=1, check.names=F)
colnames(d) <- sapply(strsplit(colnames(d), "-", fixed=T), function(x){paste(x[1:3], collapse="-")})
## THERE ARE DUPLICATES FOR FIVE PATIENTS - JUST TAKE THE FIRST ONE
d <- d[, !duplicated(colnames(d))]
d <- d[, colnames(d) %in% rownames(st) ]
d <- d[apply(d, 1, sd) != 0, ]
st <- st[ colnames(d), ]

## GET THE GENESETS
genesets <- load.gmt.data2(getFileLocation(synGet("syn2321865")))
genesets <- lapply(genesets, function(x){
  x <- x[ x != "" ]
  x <- x[ !is.na(x) ]
  intersect(x, rownames(d))
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
rownames(diffExprPvalues) <- rownames(d)
colnames(diffExprPvalues) <- colnames(st)

diffExprFCs <- sapply(diffExprResults, function(x){
  2^x$coefficients[, "factor(resp)1"]
})
rownames(diffExprFCs) <- rownames(d)
colnames(diffExprFCs) <- colnames(st)

pvalFile <- file.path(tempdir(), paste("diffProtExprPvalues-", group, "-", ds, ".tsv", sep=""))
write.table(diffExprPvalues, file=pvalFile, quote=F, sep="\t", col.names=NA)
pvalSyn <- synStore(File(path=pvalFile, parentId="syn2476109", group=group, dataset=ds, method="eBayes", stat="pvalue", evalDate=as.character(Sys.Date())), 
                    activity=Activity(name="differential protein expression",
                                      used=list(
                                        list(name=basename(code1), url=code1, wasExecuted=F),
                                        list(name=basename(code2), url=code2, wasExecuted=F),
                                        list(name=basename(code3), url=code3, wasExecuted=F),
                                        list(name=basename(code4), url=code4, wasExecuted=F),
                                        list(entity=dd, wasExecuted=F),
                                        list(entity=synGet(grpResId, downloadFile=F), wasExecuted=F),
                                        list(name=basename(thisCode), url=thisCode, wasExecuted=T)
                                      )))

fcFile <- file.path(tempdir(), paste("diffProtExprFCs-", group, "-", ds, ".tsv", sep=""))
write.table(diffExprFCs, file=fcFile, quote=F, sep="\t", col.names=NA)
fcSyn <- synStore(File(path=fcFile, parentId="syn2476109", group=group, dataset=ds, method="eBayes", stat="fc", evalDate=as.character(Sys.Date())), 
                  activity=Activity(name="differential protein expression",
                                    used=list(
                                      list(name=basename(code1), url=code1, wasExecuted=F),
                                      list(name=basename(code2), url=code2, wasExecuted=F),
                                      list(name=basename(code3), url=code3, wasExecuted=F),
                                      list(name=basename(code4), url=code4, wasExecuted=F),
                                      list(entity=dd, wasExecuted=F),
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
  op <- GSA(x=as.matrix(d), y=resp, genesets=genesets, genenames=rownames(d), resp.type="Two class unpaired", nperms=10000, minsize=3)
  op
})

gsaHiResults <- sapply(gsaResults, function(r){
  r$pvalues.hi
})
rownames(gsaHiResults) <- names(genesets)
colnames(gsaHiResults) <- colnames(st)

gsaFile <- file.path(tempdir(), paste("gsaProt-", group, "-", ds, ".tsv", sep=""))
write.table(gsaHiResults, file=gsaFile, quote=F, sep="\t", col.names=NA)
gsaSyn <- synStore(File(path=gsaFile, parentId="syn2476109", group=group, dataset=ds, method="gsa", evalDate=as.character(Sys.Date())), 
                   activity=Activity(name="geneset protein evaluation",
                                     used=list(
                                       list(name=basename(code1), url=code1, wasExecuted=F),
                                       list(name=basename(code2), url=code2, wasExecuted=F),
                                       list(name=basename(code3), url=code3, wasExecuted=F),
                                       list(name=basename(code4), url=code4, wasExecuted=F),
                                       list(entity=dd, wasExecuted=F),
                                       list(entity=synGet("syn2321865", downloadFile=F), wasExecuted=F),
                                       list(entity=synGet(grpResId, downloadFile=F), wasExecuted=F),
                                       list(name=basename(thisCode), url=thisCode, wasExecuted=T)
                                     )))

gsaHiScores <- sapply(gsaResults, function(r){
  r$GSA.scores
})
rownames(gsaHiScores) <- names(genesets)
colnames(gsaHiScores) <- colnames(st)

gsaScores <- file.path(tempdir(), paste("gsaProtScores-", group, "-", ds, ".tsv", sep=""))
write.table(gsaHiScores, file=gsaScores, quote=F, sep="\t", col.names=NA)
gsaSyn <- synStore(File(path=gsaScores, parentId="syn2476109", group=group, dataset=ds, method="gsaScores", evalDate=as.character(Sys.Date())), 
                   activity=Activity(name="geneset protein evaluation",
                                     used=list(
                                       list(name=basename(code1), url=code1, wasExecuted=F),
                                       list(name=basename(code2), url=code2, wasExecuted=F),
                                       list(name=basename(code3), url=code3, wasExecuted=F),
                                       list(name=basename(code4), url=code4, wasExecuted=F),
                                       list(entity=dd, wasExecuted=F),
                                       list(entity=synGet(grpResId, downloadFile=F), wasExecuted=F),
                                       list(entity=synGet("syn2321865", downloadFile=F), wasExecuted=F),
                                       list(name=basename(thisCode), url=thisCode, wasExecuted=T)
                                     )))
