library(synapseClient)
library(rGithubClient)
library(survival)
library(ggplot2)
library(gplots)
library(scales)
crcRepo <- getRepo("Sage-Bionetworks/crcsc")
sourceRepoFile(crcRepo, "groups/G/pipeline/JGLibrary.R")


synapseLogin("justin.guinney@sagebase.org")
source("evals/getDataFuncs.R")
source("evals/tcgaEval.R")
source("evals/evalFuncs.R")
source("groups/G/scoring//AssessmentPipelineFunc.R")


###############
# load all group results

groupResults <- lapply(names(groupFolders), getGroupResults)  
names(groupResults) <- names(groupFolders)

####
# plot subtype distribution

for(groupName in names(groupResults)){
  pdf(paste("./evalPlots/GroupSubtypeFreq_",groupName,".pdf",sep=""))
  plotSubtypeDistribution(groupResults[[groupName]])
  dev.off()
}

####
# plot dataset used by group
pdf("./evalPlots/GroupDatasetUse.pdf")
plotDatasetsUsedByGroup(groupResults)
dev.off()

######
## TCGA assessment
tcgaData <- getTCGAMutationStats()
for(groupId in names(groupResults)){
  pMatrix <- groupResults[[groupId]][["tcga_rnaseq"]]
  R <- testCinAssoc(pMatrix, tcgaData$cin)
  pdf(paste("./evalPlots/tcga/", groupId,"_cin.pdf",sep=""))
  plotCinAssoc(R=R)
  dev.off()
  R <- testMutAssoc(pMatrix, tcgaData$mutPatTbl)
  pdf(paste("./evalPlots/tcga/", groupId,"_mut.pdf",sep=""))
  plotMutAssoc(R=R)
  dev.off()
}

#############
## acquire phenodata
phenoObjs <- lapply(names(patientDatasets), getPhenoObjs)
names(phenoObjs) <- names(patientDatasets)
phenoObjs <- phenoObjs[!sapply(phenoObjs, is.na)]

# plot data size distribution
df <- data.frame(size=sapply(phenoObjs, function(x) nrow(x$data)))
f <- factor(rownames(df), levels=rownames(df)[order(df$size)])
df$dataset <- f

sp <- ggplot(df,aes(x=dataset,y=size)) + 
  geom_bar(stat="identity",position="dodge") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
pdf("./evalPlots/DatasetDist.pdf")
print(sp)
dev.off()

#############################
# assesss phenotype associations
L <- lapply(groupResults, assessPhenotypesPerGroup, phenoObjs=phenoObjs,bySubtype=FALSE)
R <- melt(L, measure.vars="value")
cn <- colnames(R)
cn[cn=="L1"] <- "group"
colnames(R) <- cn

# plot by feature
for(feature in unique(R$phenotype)){
  
  data <- R[R$phenotype==feature,]
  sp <- ggplot(data, aes(x=group, y=value)) + 
    geom_crossbar(aes(fill = factor(group),ymin=low,ymax=high),position="dodge") + 
    scale_y_continuous(trans = log10_trans(),
                       breaks = trans_breaks("log10", function(x) 10^x),
                       labels = trans_format("log10", math_format(10^.x)))
  #sp <- sp + scale_x_discrete(limits = levels(data$feature),name="Data set")
  sp <- sp + geom_hline(yintercept=.001, colour="red",linetype="dashed") + ylab("Significance")
  #sp <- sp + geom_vline(xintercept=1:(length(levels(data$ds))-1)+.5, colour="white")
  sp <- sp + theme(panel.grid.minor.x=element_blank(), 
                   panel.grid.major.x=element_blank(),
                   axis.title.x = element_blank(),
                   axis.text.x = element_blank())
  sp <- sp + facet_wrap(~dataset,ncol=4,scales="free")
  sp <- sp + ylab("Significance") + xlab("") + aes(ymin=1)
  
  pdf(paste("evalPlots/byFeature/",as.character(feature),".pdf",sep=""),width=8,height=5)
  print(sp)
  
  uds <- as.character(unique(data$dataset))
  if(feature %in% names(phenoObjs[[uds[1]]]$discreteFields)){
    par(mfrow=c(ceiling(length(uds) / 4),4))
    for(ds in uds){
      phenoObj <- phenoObjs[[ds]]
      field <- phenoObj$discreteFields[[feature]]
      tbl <- table(phenoObj$data[, field])
      textplot(cbind(tbl),cex=.9)
      title(ds)
    }
  }
  dev.off()
}

# plot by group
for(group in unique(df$group)){
  
  data <- df[df$group==group,]
  sp <- ggplot(data, aes(x=feature, y=med)) + 
    geom_crossbar(aes(fill = factor(ds),ymin=low,ymax=high),position="dodge") + 
    scale_y_continuous(trans = log10_trans(),
                       breaks = trans_breaks("log10", function(x) 10^x),
                       labels = trans_format("log10", math_format(10^.x))) + 
    scale_colour_discrete(name = "Data set") + ylab("Pvalue")
  sp <- sp + scale_x_discrete(limits = levels(data$feature),name="Data set")
  sp <- sp + geom_hline(yintercept=.001, colour="red",linetype="dashed")
  sp <- sp + geom_vline(xintercept=1:(length(levels(data$feature))-1)+.5, colour="white")
  sp <- sp + theme(panel.grid.minor.x=element_blank(), 
                   panel.grid.major.x=element_blank(),
                   axis.title.x = element_blank(),
                   axis.text.x = element_text(angle = 90, hjust = 1))
  pdf(paste("evalPlots/byGroup/",as.character(group),".pdf",sep=""),width=8,height=5)
  print(sp)
  dev.off()
}


#sp + facet_wrap(~ds,ncol=2)

############################
## comparison of overlap
globalConcordanceComparison(groupResults)

groupNames <- names(groupResults)
groupSubtypeNames <- lapply(groupNames, function(x){ 
  nsubtypes <- ncol(groupResults[[x]][[1]])
  paste(x, as.character(1:nsubtypes),sep=".")
})
names(groupSubtypeNames) <- groupNames
groupSubtypeSizes <- sapply(groupSubtypeNames, length)
Msize <- sum(groupSubtypeSizes)
A <- array(0, dim=c(Msize,Msize,length(datasets)), dimnames=list(unlist(groupSubtypeNames),unlist(groupSubtypeNames),names(datasets)))
N <- length(groupNames)
R <- list()
for(i in 1:(N-1)){
  group1 <- groupNames[i]
  for(j in (i+1):N){
    group2 <- groupNames[j]
    
    for(ds in names(datasets)){
      pMatrix1 <- groupResults[[group1]][[ds]]
      pMatrix2 <- groupResults[[group2]][[ds]]
      cat(group1, group2, ds, "\n")
      if(!is.null(pMatrix1) & !is.null(pMatrix2)){
        
        idxs <- match(clean.names(rownames(pMatrix1)), 
                      clean.names(rownames(pMatrix2)))
        pMatrix1 <- pMatrix1[!is.na(idxs),]
        pMatrix2 <- pMatrix2[na.omit(idxs),]
        
        
        tmp <- assessCrossGroupConcordance(pMatrix1,pMatrix2)
        #qts <- quantile(tmp$pvals, c(.05, .5, .95))
        #df <- data.frame(lower=qts[1],med=qts[2],high=qts[3],ds=ds, groups=paste(group1,":",group2,sep=""))
        #R[[length(R)+1]] <- df
        
        x <- sum(groupSubtypeSizes[0:(i-1)]) + 1
        y <- sum(groupSubtypeSizes[0:(j-1)]) + 1
        xrange <- x:(x+groupSubtypeSizes[i]-1)
        yrange <- y:(y+groupSubtypeSizes[j]-1)
        A[xrange,yrange, ds] <- tmp$perAtoB
        A[yrange,xrange, ds] <- t(tmp$perBtoA)
        diag(A[,,ds]) <- 1
      }
    }
  }
}
alldf <- do.call("rbind", R)
rm(R)

sp <- ggplot(alldf, aes(x=groups, y=med)) + 
  geom_crossbar(aes(fill = factor(groups),ymin=lower,ymax=high)) + 
  scale_y_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)))
sp <- sp + theme(axis.ticks=element_blank(), 
                 axis.text.x=element_blank(),
                 axis.text.y=element_text(size = 10, colour = "grey50"))
sp <- sp + facet_wrap(~ds,ncol=4,scales="free_y")
sp <- sp + ylab("") + xlab("") + aes(ymax=1)

pdf("evalPlots/cooccurencePlot.pdf",width=10,height=6)
print(sp)
dev.off()

### show percent similarity heatmap
cs <- cumsum(groupSubtypeSizes)

for(ds in dimnames(A)[[3]]){
  pdf(paste("evalPlots/clusterHeatmapByDataset/",ds,".pdf",sep=""))
  makeConcordanceFrequencyHeatmap(A[,,ds],title=ds,breaks=cs)
  dev.off()
}

dssize <- unlist(sapply(dimnames(A)[[3]], function(ds){ nrow(groupResults[[1]][[ds]]) }))
weights <- dssize / sum(dssize)
Alldata <- apply(A, c(1,2), function(x){
  sum(x[names(weights)] * weights)
})
pdf(paste("evalPlots/clusterHeatmapByDataset/ALL.pdf",sep=""))
makeConcordanceFrequencyHeatmap(Alldata,title="ALL",breaks=cs)
dev.off()