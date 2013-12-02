library(synapseClient)
library(rGithubClient)
library(survival)
library(ggplot2)
library(scales)
crcRepo <- getRepo("Sage-Bionetworks/crcsc")
sourceRepoFile(crcRepo, "groups/G/pipeline/JGLibrary.R")
synapseLogin("justin.guinney@sagebase.org")

source("groups/G/scoring//AssessmentPipelineFunc.R")


datasets <- coreDatasets
#TODO add public datasets

groupResultParents <- list(GroupA="syn2274064",GroupB="syn2274065",
                           GroupC="syn2274066",GroupD="syn2274067",
                           GroupE="syn2274069",GroupF="syn2274068",
                           GroupG="syn2274063")


# load all pmatrices
groupResults <- lapply(groupResultParents, function(parentId){
  tmp <- synapseQuery(paste('SELECT id, name FROM entity WHERE parentId=="',parentId,'"',sep=""))
  N <- nrow(tmp)
  pmatrices <- lapply(tmp$entity.id, function(synId){
    pMatrix <- read.table(synGet(synId)@filePath, sep="\t",header=T,as.is=T,row.names=1,check.names=FALSE)
  })
  synIds <- sapply(tmp$entity.name, function(x){ gsub(".*?_(syn.*?)_.*","\\1",x)})
  names(pmatrices) <- sapply(synIds, getDatanameForExprSynId)
  
  return (pmatrices)
})
groupResults <- groupResults[sapply(groupResults, length) > 0]


tmplist <- list()
for(group in names(groupResults)){
  groupResult <- groupResults[[group]]
  for(ds in names(groupResult)){
    if(ds %in% names(corePhenoObjs)){
      phenoObj <- corePhenoObjs[[ds]]
      cat(ds, group,"\n")
      tmp <- assessPhenotypes(phenoObj, groupResult[[ds]])
      tmp2 <- data.frame(group=rep(group, length(tmp)),
                         ds=rep(ds,length(tmp)),
                         feature=names(tmp),
                         t(sapply(tmp,function(x){
        quantile(x, c(.05, .5, .95))
      })),check.names=FALSE)
      tmplist[[length(tmplist) + 1]] <- tmp2
    }
  }
}
df <- do.call("rbind",tmplist)
colnames(df) <- c("group","ds","feature","low","med","high")

# plot by data set
for(ds in unique(df$ds)){
  data <- df[df$ds==ds,]
  sp <- ggplot(data, aes(x=feature, y=med)) + 
    geom_crossbar(aes(fill = factor(group),ymin=low,ymax=high),position="dodge") + 
    scale_y_continuous(trans = log10_trans(),
                       breaks = trans_breaks("log10", function(x) 10^x),
                       labels = trans_format("log10", math_format(10^.x)))
  sp <- sp + geom_hline(yintercept=.001, colour="red",linetype="dashed") + ylab("Significance")
  sp <- sp + geom_vline(xintercept=1:(sum(table(data$feature) > 0)-1)+.5, colour="white")
  sp <- sp + theme(panel.grid.minor.x=element_blank(), 
                   panel.grid.major.x=element_blank(),
                   axis.title.x = element_blank(),
                   axis.text.x = element_text(angle = 90, hjust = 1))
  pdf(paste("evalPlots/byDataset/",as.character(ds),".pdf",sep=""),width=8,height=5)
  print(sp)
  dev.off()
}

# plot by feature
for(feature in unique(df$feature)){
  
  data <- df[df$feature==feature,]
  sp <- ggplot(data, aes(x=group, y=med)) + 
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
  sp <- sp + facet_wrap(~ds,ncol=4,scales="free")
  sp <- sp + ylab("Significance") + xlab("") + aes(ymin=1)
  
  pdf(paste("evalPlots/byFeature/",as.character(feature),".pdf",sep=""),width=8,height=5)
  print(sp)
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


sp + facet_wrap(~ds,ncol=2)

############################
## comparison of overlap


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
        if(nrow(pMatrix1) !=  nrow(pMatrix2)){ 
          warning("Mismatch pmatrices.")
          idxs <- match(rownames(pMatrix1), rownames(pMatrix2))
          pMatrix1 <- pMatrix1[!is.na(idxs),]
          pMatrix2 <- pMatrix2[na.omit(idxs),]
        }
        
        tmp <- assessCrossGroupConcordance(pMatrix1,pMatrix2)
        qts <- quantile(tmp$pvals, c(.05, .5, .95))
        df <- data.frame(lower=qts[1],med=qts[2],high=qts[3],ds=ds, groups=paste(group1,":",group2,sep=""))
        R[[length(R)+1]] <- df
        
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
makeFrequencyHeatmap(Alldata,title="ALL",breaks=cs)
dev.off()