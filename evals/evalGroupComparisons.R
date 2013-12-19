
computeMutualInfoPerDataset <- function(groupResults,datasetNames){
  lvls <- lapply(groupResults, lvls4Group)
  MIs <- lapply(datasetNames, function(ds){
    R <- lapply(seq_along(groupResults), function(idx){ 
      tmp <- groupResults[[idx]][[ds]]
      if(!is.null(tmp)){ return(pmatrixToFactor(tmp, lvls[[idx]]))}
      else{ return(NULL)}
    }) 
    n <- length(R)
    MI <- unlist(sapply(1:(n-1), function(i){
      sapply((i+1):n, function(j){
        idxs <- groupMatch(clean.names(names(R[[i]])), clean.names(names(R[[j]])))
        mi.empirical(table(R[[i]][idxs[[1]]],R[[j]][idxs[[2]]]))
      })
    }))
    return(MI)
  })
  names(MIs) <- datasetNames
  df <- melt(MIs)
  medianVal <- median(na.omit(df$value))
  
  sp <- ggplot(df, aes(x=L1, y=value)) + geom_boxplot()
  sp <- sp + ylab("Mutual Information") + 
           theme(panel.grid.minor.x=element_blank(), 
                   panel.grid.major.x=element_blank(),
                   axis.title.x = element_blank(),
                   axis.text.x = element_text(angle = 90, hjust = 1))
  sp <- sp + geom_hline(yintercept=medianVal, colour="red",linetype="dashed")
  
  pdf("./evalPlots/MutualInfo_datasets.pdf",width=10,height=6)
  print(sp)
  dev.off()
}

computeMutualInfoPerGroup <- function(groupResults){
  lvls <- lapply(groupResults, lvls4Group)
  R <- aggregrateResultsPerGroup(groupResults)
  idxs <- groupMatch(lapply(R,rownames))
  
  L <- lapply(seq_along(R), function(i){
    tmp <- R[[i]][idxs[[i]],]
    F1 <- pmatrixToFactor(tmp, lvls[[i]])
    sapply(seq_along(R), function(j){
      if(i == j) { return(NA)}
      tmp <- R[[j]][idxs[[j]],]
      F2 <- pmatrixToFactor(tmp, lvls[[j]])
      mi.empirical(table(F1,F2))
    })
  })  
  names(L) <- names(groupResults)
  df <- melt(L)
  sp <- ggplot(df, aes(x=L1, y=value,fill=L1)) + geom_boxplot()
  sp <- sp + ylab("Mutual Information") + 
    theme(panel.grid.minor.x=element_blank(), 
          panel.grid.major.x=element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1))
  #sp <- sp + geom_hline(yintercept=medianVal, colour="red",linetype="dashed")
  
  pdf("./evalPlots/MutualInfo_group.pdf",width=8,height=6)
  print(sp)
  dev.off()
}


generateJSON4GroupOverlaps <- function(groupResults){
  
  ngroups <- length(groupResults)
  lvls <- lapply(groupResults, lvls4Group)
  R <- aggregrateResultsPerGroup(groupResults)
  idxs <- groupMatch(lapply(R,rownames))
  
  ALL <- list();
  for(i in 1:ngroups){
    groupId1 <- names(groupResults)[i]
    tmp <- R[[i]][idxs[[i]],]
    F1 <- pmatrixToFactor(tmp, lvls[[i]])
    for(j in 1:ngroups){
      groupId2 <- names(groupResults)[j]
      
      tmp <- R[[j]][idxs[[j]],]
      F2 <- pmatrixToFactor(tmp, lvls[[j]])
      tbl <- table(F1,F2)
      ALL[[paste(groupId1,"_",groupId2,sep="")]] = list(rownames=rownames(tbl),colnames=colnames(tbl),array=tbl)
    }
  } 
  
  str = toJSON(ALL)
  writeLines(str, con=paste("./evalPlots/d3/data/ALLgroups.json",sep=""))
  
}