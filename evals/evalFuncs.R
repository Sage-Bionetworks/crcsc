library(abind)
library(ggplot2)
library(reshape)


# converts a probability subtype matrix to a single factor
pmatrixToFactor <- function(pMatrix, lvls){
  idxs <- apply(pMatrix, 1, which.max)
  f <- factor(colnames(pMatrix)[idxs],levels=lvls)
  names(f) <- rownames(pMatrix)
  return (f)
}

lvls4Group <- function(groupResult){
  unique(unlist(lapply(groupResult,function(ds) unique(colnames(ds)))))
}

aggregrateResultsPerGroup <- function(groupResults){
  idxs <- groupMatch(lapply(groupResults, names))
  commonDatasets <- names(groupResults[[1]])[idxs[[1]]]
  commonDatasets <- setdiff(commonDatasets,"unknown")
  
  # reduce to common datasets among groups
  # concatenate results across all datasets for each group
  R <- lapply(groupResults,function(groupR){
    groupR <- groupR[names(groupR) %in% commonDatasets]
    R <- do.call("rbind",lapply(names(groupR), function(dsname){
      dsR <- groupR[[dsname]]
      #cat(length(colnames(dsR)),"\n")
      rownames(dsR) <- paste(dsname,".",clean.names(rownames(dsR)),sep="")
      return (dsR)
    }))
  })
  return(R)
}

globalConcordanceComparison <- function(groupResults){
  
  R <- aggregrateResultsPerGroup(groupResults)
  
  groupNames <- names(groupResults)
  groupSubtypeNames <- lapply(groupNames, function(x){ 
    subtypes <- colnames(groupResults[[x]][[1]])
    paste(x, subtypes,sep=".")
  })
  names(groupSubtypeNames) <- groupNames
  groupSubtypeSizes <- sapply(groupSubtypeNames, length)
  Msize <- sum(groupSubtypeSizes)
  allSubtypes <- unlist(groupSubtypeNames)
  
  A <- matrix(0, nrow=Msize, ncol=Msize,dimnames=list(allSubtypes,allSubtypes))
  diag(A) <- 1
  N <- length(groupNames)
  for(i in 1:(N-1)){
    group1 <- groupNames[i]
    for(j in (i+1):N){
      group2 <- groupNames[j]
      
      pMatrix1 <- R[[group1]]
      pMatrix2 <- R[[group2]]
      cat(group1, group2, "\n")
      if(!is.null(pMatrix1) & !is.null(pMatrix2)){
        
        idxs <- match(clean.names(rownames(pMatrix1)), 
                      clean.names(rownames(pMatrix2)))
        pMatrix1 <- pMatrix1[!is.na(idxs),]
        pMatrix2 <- pMatrix2[na.omit(idxs),]
        
        
        tmp <- assessCrossGroupConcordance(pMatrix1,pMatrix2)
        
        x <- sum(groupSubtypeSizes[0:(i-1)]) + 1
        y <- sum(groupSubtypeSizes[0:(j-1)]) + 1
        xrange <- x:(x+groupSubtypeSizes[i]-1)
        yrange <- y:(y+groupSubtypeSizes[j]-1)
        A[xrange,yrange] <- tmp$perAtoB
        A[yrange,xrange] <- t(tmp$perBtoA)
      }
    }
  }
  
  ### show percent similarity heatmap
  pdf(paste("evalPlots/clusterHeatmapByDataset/ALL_concat.pdf",sep=""),width=8,height=6)
  makeConcordanceFrequencyHeatmap(A,title="ALL concatenated",breaks=cumsum(groupSubtypeSizes))
  dev.off()
}

makeConcordanceFrequencyHeatmap <- function(A, title, breaks){
  m <- melt(A)
  
  p <- ggplot(m, aes(X1, X2)) + 
    geom_tile(aes(fill = value),colour = "white")
  #scale_fill_gradient(low = "white",high = "steelblue")
  p <- p + geom_vline(xintercept=(breaks + .5), lwd=1)
  p <- p + geom_hline(yintercept=(breaks + .5), lwd=1)
  p <- p + xlab("") + 
    ylab("") +
    theme(axis.ticks=element_blank(), 
          axis.text.x=element_text(size = 10, angle = 270, hjust = 0, colour = "grey50"),
          axis.text.y=element_text(size = 10, colour = "grey50")) +
    ggtitle(title)
  print(p)
}

clean.names <- function(x){
  x <- gsub("(.*)\\.CEL.*","\\1",x)
  x <- gsub("^X(\\d.*)","\\1",x)
  x <- gsub("\\.","-", x)
  gsub("^(GSM\\d*)\\D.*","\\1",x)
}

assessPhenotypesPerGroup <- function(groupResult, phenoObjs, bySubtype){

  groupResultDS <- groupResult[names(groupResult) %in% names(phenoObjs)]
  L <- lapply(names(groupResultDS), function(ds){
    cat("-------",ds,"---------\n")
    tmp <- assessPhenotypePerPmatrix(phenoObjs[[ds]], groupResultDS[[ds]], bySubtype)
    if(bySubtype){
      # summarize across phenotype
      R <- lapply(tmp, function(fmatrix){
        # summarize across subtypes
        apply(fmatrix, 2, function(subtypeM){
          # summarize into quantiles
          q <- quantile(subtypeM, c(.05, .5, .95))
          names(q) <- c("low","med","high")
          return (q)
        })
      })
      df <- melt(R)
      colnames(df) <- c("quant","subtype","value","phenotype")
      return (df)
      
    }else{
      #browser()
      df <- data.frame(
        #ds=rep(ds,length(tmp)),
        phenotype=names(tmp),
        t(sapply(tmp,function(x){
          q <- quantile(x, c(.05, .5, .95))
          names(q) <- c("low","med","high")
          q
        })),check.names=FALSE)
      return(df)
    }
  })
 # browser()
  names(L) <- names(groupResultDS)
  df <- melt(L, measure.vars="med")
  cn <- colnames(df)
  cn[length(cn)] <- "dataset"
  colnames(df) <- cn
  
  return (df)  
}

assessPhenotypePerPmatrix <- function(vObj, pMatrix, bySubtype=FALSE){
  #browser()
  phenotypeTbl <- vObj$data
  
  idxs <- groupMatch(clean.names(rownames(phenotypeTbl)),clean.names(rownames(pMatrix)))
  if(length(idxs[[1]]) == 0){
    cat("No match between phenotype and pMatrix\n")
    return()
  }
  phenotypeTbl.m <- phenotypeTbl[idxs[[1]],]
  pMatrix.m <- pMatrix[idxs[[2]],]
  
  #browser()
  
  discreteResults <- lapply(vObj$discreteFields, function(fieldName){
    cat(fieldName,"\n")
    discretePhenotype <- phenotypeTbl.m[, fieldName]
    mask <- !is.na(discretePhenotype)
    test.discrete(discretePhenotype[mask], pMatrix.m[mask,],bySubtype)
  })
  continuousResults <- lapply(vObj$continuousFields, function(fieldName){
    cat(fieldName,"\n")
    continuousPhenotype <- as.numeric(phenotypeTbl.m[, fieldName])
    mask <- !is.na(continuousPhenotype)
    test.continuous(continuousPhenotype[mask], pMatrix.m[mask,],bySubtype)
  })
  censoredResults <- lapply(vObj$censoredFields, function(pair){
    time <- as.numeric(phenotypeTbl.m[, pair[1]])
    status <- phenotypeTbl.m[, pair[2]]
    mask <- !is.na(time) & !is.na(status)
    test.censored(Surv(time[mask],status[mask]), pMatrix.m[mask,],bySubtype)
  })
  return (c(discreteResults, continuousResults,censoredResults))
}

assessCrossGroupConcordance <- function(pmatrix1, pmatrix2, ...){
  Nc1 <- ncol(pmatrix1)
  Nc2 <- ncol(pmatrix2)
  
  smatrix1 <- samplePMatrix(pmatrix1, ...)
  smatrix2 <- samplePMatrix(pmatrix2, ...)
  p <- nrow(smatrix1)
  
  compute_percent_overlap <- function(isAtoB){
    tmp <- do.call("abind", c(lapply(1:p, function(idx){
      tmp2 <- table(factor(smatrix1[idx,],levels=c(1:Nc1)),
                    factor(smatrix2[idx,],levels=c(1:Nc2)))
      if(isAtoB){ tmp3 <- tmp2 %*% diag(1/colSums(tmp2)) }
      else { tmp3 <- t(t(tmp2) %*% diag(1/rowSums(tmp2))) }
      tmp3[is.na(tmp3)] <- 0
      return(tmp3)
    }), along=3))
    apply(tmp, c(1,2), mean)
  }
  
  perAtoB <- compute_percent_overlap(TRUE)
  perBtoA <- compute_percent_overlap(FALSE)
  
  if(any(perAtoB > 1)){
    browser()
  }
  
  pvals <- sapply(1:p, function(idx){
    # use approximation; fisher test too slow
    chisq.test(smatrix1[idx,],smatrix2[idx,])$p.value
  })
  return (list(pvals=pvals,perAtoB=perAtoB,perBtoA=perBtoA))
}

samplePMatrix <- function(pMatrix, nSampling=100, noise=10^-3){
  nSamples <- nrow(pMatrix)
  nCats <- ncol(pMatrix)
  pMatrix <- pMatrix + matrix(abs(rnorm(nSamples * nCats,0,noise)),nrow=nSamples)
  
  subtypeCalls <- apply(pMatrix, 1, function(x){
    csum <- c(0,cumsum(x))
    csum <- csum/csum[length(csum)]
    as.numeric(cut(runif(nSampling), breaks=csum,include.lowest=TRUE))
  })
  return (subtypeCalls)
}

test.censored <- function(survObj, pMatrix,bySubtype=FALSE,...){
  subtypeCalls <- samplePMatrix(pMatrix,...)
  if(bySubtype){
    subtypes <- names(pMatrix)
    R <- sapply(1:length(subtypes),function(subtypeIdx){
      pvals <- apply(subtypeCalls, 1, function(x){
        sdf <- survdiff(survObj ~ factor(x==subtypeIdx))
        p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
        return(p.val)
      })
      return (pvals)
    })
    colnames(R) <- subtypes
    return (R)
  }else{
    pvals <- apply(subtypeCalls, 1, function(x){
      sdf <- survdiff(survObj ~ factor(x))
      p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
      return(p.val)
    })
    return (pvals)
  }
}

test.continuous <- function(continuousPhenotype, pMatrix,bySubtype=FALSE,...){
  subtypeCalls <- samplePMatrix(pMatrix,...)
  if(bySubtype){
    subtypes <- names(pMatrix)
    R <- sapply(1:length(subtypes),function(subtypeIdx){
      pvals <- apply(subtypeCalls, 1, function(x){
        wilcox.test(continuousPhenotype ~ factor(x==subtypeIdx))$p.value
      })
      pvals
    })
    colnames(R) <- subtypes
    return (R)
  }else{
    pvals <- apply(subtypeCalls, 1, function(x){
      kruskal.test(continuousPhenotype, factor(x))$p.value
    })
    return (pvals)
  }
}

test.discrete <- function(discretePhenotype, pMatrix,bySubtype=FALSE,...){
  subtypeCalls <- samplePMatrix(pMatrix,...)
  
  if(bySubtype){
    subtypes <- names(pMatrix)
    R <- sapply(1:length(subtypes),function(subtypeIdx){
      pvals <- apply(subtypeCalls, 1, function(x){
        tryCatch({
          chisq.test(discretePhenotype,factor(x==subtypeIdx))$p.value
        }, error=function(e){
          return(1)
        })
      })
      pvals
    })
    colnames(R) <- subtypes
    return (R)
  }else{
    pvals <- apply(subtypeCalls, 1, function(x){
      #fisher.test(discretePhenotype, factor(x),workspace=2e+07,hybrid=TRUE)$p.value
      chisq.test(discretePhenotype,factor(x))$p.value
    })
    return (pvals)
  }
  
}


# for weighted meta-analysis
Stouffer.test <- function(p, w) { # p is a vector of p-values
  if (missing(w)) {
    w <- rep(1, length(p))/length(p)
  } else {
    if (length(w) != length(p))
      stop("Length of p and w must equal!")
  }
  Zi <- qnorm(1-p) 
  Z  <- sum(w*Zi)/sqrt(sum(w^2))
  p.val <- 1-pnorm(Z)
  return(p.val)
}

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
