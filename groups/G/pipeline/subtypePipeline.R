library(synapseClient)
library(rGithubClient)
library(affy)
library(hgu133plus2.db)
library(hgu133a2.db)
library(org.Hs.eg.db)

crcRepo <- getRepo("Sage-Bionetworks/crcsc")
sourceRepoFile(crcRepo, "groups/G/pipeline/JGLibrary.R")
code1 <- getPermlink(crcRepo, "groups/G/pipeline/JGLibrary.R")
sourceRepoFile(crcRepo, "groups/G/pipeline/subtypePipelineFuncs.R")
code2 <- getPermlink(crcRepo, "groups/G/pipeline/subtypePipelineFuncs.R")
thisScript <- getPermlink(crcRepo, "groups/G/pipeline/subtypePipeline.R")

# password will be request after calling this
synapseLogin("justin.guinney@sagebase.org")

coreExprList <- list(
  agendia_gse42284=toEntrez("syn2192792",NULL, discoverprint_19742_Map,loader=agendia_data_loader, filename="GSE42284_normalized_data_matrix.txt"),
  agendia_ico208=toEntrez("syn2192796",NULL, discoverprint_19742_Map,loader=agendia_data_loader, filename="ICO208_normalized_data.txt"),
  agendia_vhb70=toEntrez("syn2192799",NULL,discoverprint_32627_Map, loader=agendia_data_loader, filename="VHB70_normalized_data.txt"),
  kfsyscc=toEntrez("syn2169565",NULL,u133plus2Map),
  french=toEntrez("syn2171434",NULL,u133plus2Map),
  amc_ajccii=toEntrez("syn2159423",NULL,u133plus2Map,sep=","),
  nki_az=toEntrez("syn2176657",NULL,u133plus2Map),
  petacc3=toEntrez("syn2175581",NULL,petaccMap),
  tcga_rnaseq=toEntrez("syn2161141",NULL,symbolMap))


publicExprList <- list(
  gse10961=toEntrez("syn2177194","syn2177195"),
  gse13067=toEntrez("syn2177888","syn2177889"),
  gse13294=toEntrez("syn2177894","syn2177895"),
  gse14333=toEntrez("syn2181079","syn2181006"),
  gse15960=toEntrez("syn2177199","syn2177200"),
  gse17537=toEntrez("syn2178128","syn2178129"),
  gse20916=toEntrez("syn2177899","syn2177898"),
  gse2109=toEntrez("syn2177169","syn2177168"),
  gse23878=toEntrez("syn2177902","syn2178063"),
  gse37892=toEntrez("syn2178082","syn2178089"),
  gse4107=toEntrez("syn2177179","syn2177180"),
  gse4183=toEntrez("syn2177187","syn2177188"),
  gse8671=toEntrez("syn2181088","syn2181090"))

allData <- c(coreExprList,publicExprList)

plotCorrOfCorr(coreExprList, publicExprList)

# find common set of entrez genes and subset on those
midxs <- groupMatch(lapply(allData, rownames))
mExprList <- lapply(1:length(midxs), function(i){
  allData[[i]][midxs[[i]],]
})
names(mExprList) <- names(allData)

# use KFSYSCC data as training data.
trainX <- exprs(mExprList[["kfsyscc"]])

# select genes with highest mad (quantile > 60%)
madX <- apply(trainX, 1, mad)
trainX <- trainX[madX > quantile(madX, .6),]

set.seed(2013)
## find clusters in training set (fixed at 4)
clusterX <- findClusters(trainX,k=4)
table(clusterX$groups)

# find top differentially expressed genes between clusters
diffGenes <- findDiffGenes(trainX, clusterX$groups)
trainXs <- trainX[rownames(trainX) %in% diffGenes,]

# principal component plots of clusters
par(mfrow=c(2,2))
pc1 <- svd(trainX - rowMeans(trainX))
plot((pc1$d^2 / sum(pc1$d^2))[1:100],main="Eig distribution",ylab="% var")
plot(pc1$v[,1], pc1$v[,2], col=factor(clusterX$groups), xlab="PC1",ylab="PC2",
     main="PC Plot (Pre Gene Selection)",pch=19,cex=.8)

pc2 <- svd(trainXs - rowMeans(trainXs))
plot((pc2$d^2 / sum(pc2$d^2))[1:100],main="Eig distribution",ylab="% var")
plot(pc2$v[,1], pc2$v[,2], col=factor(clusterX$groups), xlab="PC1",ylab="PC2",
     main="PC Plot (Post Gene Selection)",pch=19,cex=.8)

# build classifier
train <- list(x=trainXs,y=factor(clusterX$groups))
fit <- pamr.train(data=train)
cvfit <- pamr.cv(fit, train)
opt.threshold <- cvfit$threshold[which(cvfit$error == min(cvfit$error))]
if(length(opt.threshold) > 1){
  opt.threshold <- opt.threshold[length(opt.threshold)]
}
pdf("cv.pdf")
pamr.plotcv(cvfit)
dev.off()

Xmean <- rowMeans(trainXs)
Xsd <- apply(trainXs, 1, sd)

# predict on rest of data sets
predList <- lapply(mExprList, function(eset){
  X <- exprs(eset)
  idxs <- match(rownames(trainXs), rownames(X))
  Xtest <- X[idxs,]
  # normalize data to have same mean and sd as training data
  Xntest <- normalize_to_X(Xmean,Xsd, Xtest)
  pred <- pamr.predict(fit, Xntest, opt.threshold, type="posterior")
  return (pred)
})


###################################################
# store all the predictions in synapse as tsv files

synResultDir <- "syn2229267"
groupName <- "GroupG"
trainSynId <- experimentData(allData[["kfsyscc"]])@name

lapply(names(predList), function(datasetName){
  tmp <- round(predList[[datasetName]],digits=4)
  colnames(tmp) <- paste("cluster", colnames(tmp),sep="")
  pred <- data.frame(sampleName=rownames(tmp), tmp)
  synId <- experimentData(allData[[datasetName]])@name
  filePath <- file.path(tempdir(), paste(groupName,"_",synId,"_",datasetName,".tsv",sep=""))
  write.table(pred, file=filePath,sep="\t",quote=FALSE,row.names=FALSE)
  
  # store results in synapse 
  synFile  <- File(path=filePath, parentId=synResultDir)
  synFile  <- synStore(synFile, used=list(list(entity=trainSynId, wasExecuted=F),
                                          list(entity=synId, wasExecuted=F),
                                          list(url=code1, name=basename(code1), wasExecuted=F),
                                          list(url=code2, name=basename(code2), wasExecuted=F),
                                          list(url=thisScript, name=basename(thisScript), wasExecuted=T)))
  unlink(filePath)
})
