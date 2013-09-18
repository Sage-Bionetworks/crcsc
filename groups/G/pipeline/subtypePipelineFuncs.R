
u133plus2Map <- function(probeIds){mget(probeIds,hgu133plus2ENTREZID,ifnotfound=NA) }
u133a2Map <- function(probeIds){mget(probeIds,hgu133a2ENTREZID,ifnotfound=NA) }
symbolMap <- function(probeIds){mget(probeIds,org.Hs.egSYMBOL2EG,ifnotfound=NA) }
petaccMap <- function(probeIds){
  file <- getFileLocation(synGet("syn2199825"))
  tbl <- read.table(file, sep="\t",as.is=TRUE,header=TRUE,comment="")
  idxs <- match(probeIds, tbl$ProbesetID)
  stopifnot(!any(is.na(idxs)))
  return(tbl$Entrez.GeneID[idxs])
}
discoverprint_19742_Map <- function(probeIds){
  file <- getFileLocation(synGet("syn2192791"))
  tbl <- read.table(gzfile(file), sep="\t",header=T,as.is=TRUE)
  idxs <- match(probeIds, tbl$probe_id)
  stopifnot(!any(is.na(idxs)))
  geneNames <- tbl$symbol[idxs]
  mget(geneNames,org.Hs.egSYMBOL2EG,ifnotfound=NA)
}
discoverprint_32627_Map <- function(probeIds){
  file <- getFileLocation(synGet("syn2192801"))
  tbl <- read.table(file, sep="\t",header=TRUE,as.is=TRUE,comment="",quote="")
  idxs <- match(probeIds, tbl$ProbeName)
  stopifnot(!any(is.na(idxs)))
  geneNames <- tbl$GeneName[idxs]
  mget(geneNames,org.Hs.egSYMBOL2EG,ifnotfound=NA)
}

agendia_data_loader <- function(synId, filename){
  require(impute)

  file <- getFileLocation(synGet(synId))
  con <- unz(file,filename=filename)
  tbl <- read.table(con,sep="\t",header=TRUE,as.is=TRUE)
  duplicated.idxs <- which(duplicated(tbl$probe_id))
  if(length(duplicated.idxs) > 0){
    cat("Found ", length(duplicated.idxs), " duplicated probes. Removing...\n")
    tbl <- tbl[-duplicated.idxs,]
  }
  M <- as.matrix(tbl[,-1])
  rownames(M) <- tbl[,1]
  Mn <- impute.knn(M)$data
  return(Mn)
}

fastLoad <- function(synId,sep="\t",...){
  file <- getFileLocation(synGet(synId))
  
  if(grepl(".*gz$",file)){ con <- gzfile(file)
  }else if(grepl(".*zip$",file)){ con <- unz(file,...)
  }else{ con <- open(file)}
  tmp <- read.table(con, sep=sep,nrows=1)
 
  if(grepl(".*gz$",file)){ con <- gzfile(file)
  }else if(grepl(".*zip$",file)){ con <- unz(file,...)
  }else{ con <- open(file)}
  
  ncolumns <- length(tmp)
  colClasses <- c("character", rep("numeric",ncolumns-1))
  tbl <- as.matrix(read.table(file,sep=sep,header=TRUE,row.names=1,colClasses=colClasses,comment=""))
  close(file)
  return (tbl)
}


toEntrez <- function(synId, synAnnId=NULL, entrezMapper=u133plus2Map,sep="\t", loader=NULL,...){

  getGSMId <- function(ids){ gsub("(GSM\\d*).*","\\1",ids) }
  
  cat(synId,"\n")  
  if(!is.null(loader)){
    expr <- loader(synId, ...)
  }else{
    expr <- fastLoad(synId,sep,...)
  }
  
  entrezIds <- entrezMapper(rownames(expr))
  #browser()
  
  # in case of rare multiple mappings, select 1st
  entrezIds <- unlist(sapply(entrezIds, function(x){ ifelse(length(x) > 1, x[1], x)}))
  
  exprEntrez <- data.matrix(t(apply(expr[!is.na(entrezIds),], 1, as.numeric)))
  colnames(exprEntrez) <- getGSMId(colnames(expr))
  entrezNoNA <- entrezIds[!is.na(entrezIds)]
  
  bestIdxs <- oneToManyOperation(unique(entrezNoNA), entrezNoNA, function(idxs){
    sidxs <- order(apply(exprEntrez[idxs,,drop=FALSE], 1, mad), decreasing=T)
    return(idxs[sidxs[1]])
  })
  
  bestExpr <- exprEntrez[bestIdxs,]
  rownames(bestExpr) <- entrezNoNA[bestIdxs]
  
  experimentData <- new("MIAME", name=synId)
  
  if(!is.null(synAnnId)){
    file <- getFileLocation(synGet(synAnnId))
    pdata <- read.table(file, sep="\t",header=TRUE,row.names=1)
    rownames(pdata) <- getGSMId(rownames(pdata))
    idxs <- match(colnames(bestExpr), rownames(pdata))
   
    eset <- new("ExpressionSet",exprs=bestExpr,experimentData=experimentData,
                phenoData=new("AnnotatedDataFrame",pdata[idxs,]))
  }else{
    eset <- new("ExpressionSet",exprs=bestExpr,experimentData=experimentData)
  }
  
  #if(!is.null(postFunc)){
  #  eset <- postFunc(eset)
  #}
  return (eset) 
}


findClusters <- function(X, method="hclust", k=3, 
                         numbootstraps=10, samplingPercent=.6){
  n <- ncol(X)
  p <- nrow(X)
  

  numSampledFeatures <- round(samplingPercent * p)
 
  M <- array(dim=c(n, n, numbootstraps))
  for(i in 1:numbootstraps){
    Xs <- X[sample(p, numSampledFeatures),]
    dissimilarity <- 1 - cor(Xs, method="spearman")
    #rbf <- rbfdot(sigma = 1)
    #dissimilarity  <- kernelMatrix(rbf, t(X))
    distance <- as.dist(dissimilarity)
    fit <- hclust(distance, method="ward") 
    groups <- cutree(fit, k=k)
    M[,,i] <- matrix(crossprod(t(as.matrix(groups))) %in% c(1,4,9),nrow=length(groups))
  }
  
  # generate new distance matrix based on number of times samples co-cluster (normalized to percentage)
  D <- 1 - apply(M, c(1,2), sum) / numbootstraps
  
  # generate new hierarchical cluster using co-occurence matrix
  distance <- as.dist(D)
  fit <- hclust(distance, method="ward") 
  groups <- cutree(fit, k=k)
  
  scores <- sapply(1:k, function(i){
    mask <- groups==i
    sapply(1:n, function(j){
      mean(D[j, mask]) } )
  })
  nscores <- scores / rowSums(scores)
  
  return (list(groups=groups, scores=nscores))  
}

# find differentially expressed genes between sub-types
findDiffGenes <- function(X, groups, top=100){
  f <- factor(groups)
  topGenes <- lapply(levels(f), function(lvl){
    gf <- factor(lvl == f)
    adjp <- p.adjust(apply(X, 1, function(x){ wilcox.test(x ~ gf)$p.value}), method="BH")
    # just take top genes, rather than filter on adjusted pvalue
    return (rownames(X)[order(adjp)[1:top]])
  })
  unique(do.call("c", topGenes))
}

plotCorrOfCorr <- function(coreEsets, publicEsets,method="spearman", subsample=.2){
  
  midxs <- groupMatch(lapply(c(coreEsets, publicEsets), rownames))
  coreIdxs <- midxs[1:length(coreEsets)]
  pubIdxs <- midxs[-(1:length(coreEsets))]
  
  N <-  length(midxs[[1]])
  sidxs <- sample(N, round(subsample * N))
  coreC <- lapply(1:length(coreEsets), function(i){
    cor(t(exprs(coreEsets[[i]])[coreIdxs[[i]][sidxs],]),method=method)
  })
  pubC <- lapply(1:length(publicEsets), function(i){
    cor(t(exprs(publicEsets[[i]])[pubIdxs[[i]][sidxs],]),method=method)
  })
  
  compute.cvec <- function(c1, c2){
    sapply(1:ncol(c1), function(j){ cor(c1[,j],c2[,j],method=method) })
  }

  for(i in 1:length(coreEsets)){
    pdf(paste("plots/", names(coreEsets)[i],"_densityPlot.pdf",sep=""),width=8,height=4)
    par(oma=c(0,0,2,0))
    C <- coreC[[i]]
    core.cvecs <- lapply(coreC[-i], compute.cvec, C)  
    pub.cvecs <- lapply(pubC, compute.cvec, C)
    df <- data.frame(dens=do.call("c", core.cvecs), groups=rep(names(coreEsets)[-i], each=ncol(C)))
    p1 <- densityplot(~dens, data=df, groups=groups,plot.points=FALSE,ref=TRUE,auto.key=list(space="top"),xlim=c(-.5,1))
    df <- data.frame(dens=do.call("c", pub.cvecs[1:6]), groups=rep(names(publicEsets)[1:6], each=ncol(C)))
    p2 <- densityplot(~dens, data=df, groups=groups,plot.points=FALSE,ref=TRUE,auto.key=list(space="top"),xlim=c(-.5,1))
    df <- data.frame(dens=do.call("c", pub.cvecs[7:13]), groups=rep(names(publicEsets)[7:13], each=ncol(C)))
    p3 <- densityplot(~dens, data=df, groups=groups,plot.points=FALSE,ref=TRUE,auto.key=list(space="top"),xlim=c(-.5,1))
    print(p1, position=c(0, 0, .33, 1),more=TRUE)
    print(p2, position=c(.3, 0, .66, 1),more=TRUE)
    print(p3, position=c(.66, 0, 1, 1))
    dev.off()
  }
}

