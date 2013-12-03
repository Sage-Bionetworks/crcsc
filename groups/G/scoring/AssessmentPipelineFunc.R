library(abind)
library(ggplot2)

makeConcordanceFrequencyHeatmap <- function(A, title, breaks){
  m <- melt(A)
  
  p <- ggplot(m, aes(X1, X2)) + 
    geom_tile(aes(fill = value),colour = "white") + 
    scale_fill_gradient(low = "white",high = "steelblue")
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
  x <- gsub("^X(\\d.*)","\\1",x)
  gsub("\\.","-", x)
}


assessPathwaysKS <- function(exprDatasets, pMatrices){
  
  gs <- load.gmt.data(gmt.file.path="./c2.cp.v4.0.entrez.gmt")
  cat("computing geneset idxs for each data set\n")
  gsetIdxs <- lapply(names(exprDatasets), function(ds){
    genes <- rownames(exprDatasets[[ds]])
    lapply(gs, function(x){
      which(genes %in% x) 
    })
  })
  names(gsetIdxs) <- names(exprDatasets)
  
  cat("compute enrichment per pMatrix\n")
  sapply(names(pMatrices), function(ds){
    expr <- exprDatasets[[ds]]
    pMatrix <- pMatrices[[ds]]
    gsetIdxForDS <- gsetIdxs[[ds]]
    idxs <- groupMatch(clean.names(colnames(expr)),clean.names(rownames(pMatrix)))
    if(length(idxs[[1]]) == 0){
      cat("No match between gsva es and pMatrix\n")
      return()
    }
    expr.m <- expr[,idxs[[1]]]
    pMatrix.m <- pMatrix[idxs[[2]],]
    f <- factor(apply(pMatrix.m, 1, function(x) order(x,decreasing=TRUE)[1]))
    
    pvals <- sapply(levels(f), function(lvl){
      mask <- lvl == f
      ranks <- order(eBayes(lmFit(expr.m, model.matrix(~mask)))$p.value[,2])
      sapply(gsetIdxForDS, function(idxs){
        ks.test(ranks[idxs],ranks[-idxs])$p.value
      })
    })
  })
        
}

assessPathways <- function(gse, pMatrices){  
    
    # match datasets with pmatrices
    idxs <- groupMatch(names(gse), names(pMatrices))
    gse.m <- gse[idxs[[1]]]
    pMatrices.m <- pMatrices[idxs[[2]]]
    

    gs.pvals <- sapply(names(gse.m), function(ds){
      gsvaES <- gse.m[[ds]]
      pMatrix <- pMatrices[[ds]]
      idxs <- groupMatch(clean.names(colnames(gsvaES)),clean.names(rownames(pMatrix)))
      if(length(idxs[[1]]) == 0){
        cat("No match between gsva es and pMatrix\n")
        return()
      }
      gsvaES.m <- gsvaES[,idxs[[1]]]
      pMatrix.m <- pMatrix[idxs[[2]],]
      
      # for speed, we select max prob rather than sample
      f <- factor(apply(pMatrix.m, 1, function(x) order(x,decreasing=TRUE)[1]))    
      fit <- eBayes(lmFit(gsvaES, model.matrix(~f)))
      pvals <- topTableF(fit,sort="none",n=Inf)$adj.P.Val
  
      names(pvals) <- rownames(gsvaES.m)
      #apply(pvals, 1, median)
    })
    
    # weight by data set size
    tmp <- sqrt(sapply(gse.m, ncol))
    weights <- tmp / sum(tmp)
    
    meta_pvals <- apply(gs.pvals, 1, function(p) Stouffer.test(p,weights) )
}

assessPhenotypes <- function(vObj, pMatrix){
  #browser()
  phenotypeTbl <- vObj$data
  
  clean.names <- function(x){
    x <- gsub("^X(\\d.*)","\\1",x)
    gsub("\\.","-", x)
  }
  
  idxs <- groupMatch(clean.names(rownames(phenotypeTbl)),clean.names(rownames(pMatrix)))
  if(length(idxs[[1]]) == 0){
    cat("No match between phenotype and pMatrix\n")
    return()
  }
  phenotypeTbl.m <- phenotypeTbl[idxs[[1]],]
  pMatrix.m <- pMatrix[idxs[[2]],]
  
  discreteResults <- lapply(vObj$discreteFields, function(fieldName){
    cat(fieldName,"\n")
    discretePhenotype <- phenotypeTbl.m[, fieldName]
    mask <- !is.na(discretePhenotype)
    test.discrete(discretePhenotype[mask], pMatrix.m[mask,])
  })
  continuousResults <- lapply(vObj$continuousFields, function(fieldName){
    cat(fieldName,"\n")
    continuousPhenotype <- phenotypeTbl.m[, fieldName]
    mask <- !is.na(continuousPhenotype)
    test.continuous(continuousPhenotype[mask], pMatrix.m[mask,])
  })
  censoredResults <- lapply(vObj$censoredFields, function(pair){
    time <- phenotypeTbl.m[, pair[1]]
    status <- phenotypeTbl.m[, pair[2]]
    mask <- !is.na(time) & !is.na(status)
    test.censored(Surv(time[mask],status[mask]), pMatrix.m[mask,])
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

test.censored <- function(survObj, pMatrix,...){
  subtypeCalls <- samplePMatrix(pMatrix,...)
  
  pvals <- apply(subtypeCalls, 1, function(x){
    sdf <- survdiff(survObj ~ factor(x))
    p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
    return(p.val)
  })
  return (pvals)
}

test.continuous <- function(continuousPhenotype, pMatrix,...){
  subtypeCalls <- samplePMatrix(pMatrix,...)
  
  pvals <- apply(subtypeCalls, 1, function(x){
    kruskal.test(continuousPhenotype, factor(x))$p.value
  })
  return (pvals)
}

test.discrete <- function(discretePhenotype, pMatrix,...){
  subtypeCalls <- samplePMatrix(pMatrix,...)
  
  pvals <- apply(subtypeCalls, 1, function(x){
    #fisher.test(discretePhenotype, factor(x),workspace=2e+07,hybrid=TRUE)$p.value
    chisq.test(discretePhenotype,factor(x))$p.value
  })
  return (pvals)
}

testSangerDrugSensitivity <- function(pMatrix){
  drugResponse <- loadEntity("syn1807986")$objects$sangerAUC@data
  subtypeCalls <- samplePMatrix(pMatrix,...)
}


build.gsva.matrices <- function(datasets){
  require(GSVA)
  require(multicore)
  gsets <- load.gmt.data("./c2.cp.v4.0.entrez.gmt")
  ges <- lapply(datasets, function(ds){
    entrezIds <- featureNames(ds)
    gsetIdxs <- lapply(gsets, function(gs){ 
      na.omit(match(gs, entrezIds))
    })
    gsva(exprs(ds), gsetIdxs,min.sz=10,max.sz=200,parallel.sz=10)$es.obs  
  })
  save(ges,file="./GSVA_coredatasets.rda")
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