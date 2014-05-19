require(synapseClient)
require(rGithubClient)
require(ggplot2)

## SOURCE IN BACKGROUND FUNCTIONS FROM JG
crcRepo <- getRepo("Sage-Bionetworks/crcsc")
sourceRepoFile(crcRepo, "groups/G/pipeline/JGLibrary.R")
code1 <- getPermlink(crcRepo, "groups/G/pipeline/JGLibrary.R")
sourceRepoFile(crcRepo, "groups/G/pipeline/subtypePipelineFuncs.R")
code2 <- getPermlink(crcRepo, "groups/G/pipeline/subtypePipelineFuncs.R")

## SOURCE CODE TO READ IN DATA
sourceRepoFile(crcRepo, "evals/getDataFuncs.R")
code3 <- getPermlink(crcRepo, "evals/getDataFuncs.R")

## THIS SCRIPT
thisCode <- getPermlink(crcRepo, "evals/summarizeGenesetsConsensus.R")
resFold <- synGet("syn2476109")

## QUERY FOR OUR RESULTS
resQ <- synapseQuery("SELECT id, name, group, dataset, method, stat, evalDate FROM file WHERE parentId=='syn2476109'")
gsQ <- resQ[ resQ$file.method == "gsa", ]

## gse4183 ONLY HAS FOUR SUBTYPES CALLS
gsQ <- gsQ[ gsQ$file.dataset != "gse4183", ]
gp <- "cms4"


## PULL IN ALL OF THE RESULT MATRICES FROM THE GENESET ANALYSES
gsaDat <- lapply(as.list(gsQ$file.id), function(x){
  a <- synGet(x)
  read.delim(getFileLocation(a), as.is=T, header=T, row.names=1)
})
names(gsaDat) <- gsQ$file.dataset


## PLOTS PER GROUP
## FISHER META ANALYSIS P-VALUE
# 1 - pchisq(-2 * sum(log(pvals)),2 * length(pvals))
chsq <- -2*Reduce("+", lapply(gsaDat, log))
pval <- apply(chsq, c(1,2), function(x){
  1-pchisq(x, 2*length(gsaDat))
})
x <- -1*log10(pval)
x <- apply(x, c(1,2), function(y){
  min(y, 20)
})

plotDF <- data.frame(geneset = rep(rownames(x), ncol(x)),
                     subtype = rep(colnames(x), each=nrow(x)),
                     count = as.numeric(x))
p <- ggplot(data=plotDF, aes(x=subtype, y=count, fill=subtype)) +
  geom_bar(stat="identity") + xlab("") + ylab("-log10(fisher meta analysis pval)") + ggtitle(gp) +
  facet_wrap(facets=(~ geneset), ncol=3)

plotPath <- file.path(tempdir(), paste("genesetsMetaAnalysis-", gp, ".png", sep=""))
png(plotPath, width=900, height=600)
show(p)
dev.off()
synPlot <- synStore(File(path=plotPath, parentId="syn2476110"),
                    activity=Activity(name="geneset plots",
                                      used=list(
                                        list(url=thisCode, name=basename(thisCode), wasExecuted=T),
                                        list(entity=resFold, wasExecuted=F))))
