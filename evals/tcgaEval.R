

library(plyr)
########### evaluation of TCGA mutation / copy-number / clinical associations
barcode2patid <- function(barcode){ gsub("(TCGA-.*?-.*?)-.*", "\\1", barcode) }


getTCGAMutationStats <- function(){
  #
  mutTbl <- read.table("data//coadread_cleaned.maf",sep="\t",header=TRUE,as.is=TRUE,comment="",quote="")
  mutTbl$patid <- barcode2patid(mutTbl$Tumor_Sample_Barcode)
  patids <- unique(mutTbl$patid)
  
  # only analyze genes present in 10% of patients
  geneCts <- ddply(mutTbl, c("Hugo_Symbol"), function(df) { 
    nrow(df)     
  })
  
  cin <- ddply(mutTbl, c("patid"), function(df){
    nrow(df)
  })
  rownames(cin) <- cin$patid
  cin <- cin[,2,drop=FALSE]
  
  keepGenes <- geneCts$Hugo_Symbol[geneCts$V1 > length(patids) * .05]
  mutTblShort <- mutTbl[mutTbl$Hugo_Symbol %in% keepGenes,]
  mutPatTbl <- ddply(mutTblShort, c("Hugo_Symbol"), function(df) { 
    as.numeric(patids %in% df$patid)     
  })
  rownames(mutPatTbl) <- mutPatTbl[,1]
  mutPatTbl <- mutPatTbl[,-1]
  colnames(mutPatTbl) <- patids
  return(list(mutPatTbl=mutPatTbl,cin=cin))
}

testCinAssoc <- function(pMatrix, cin){
  f <- pmatrixToFactor(pMatrix)
  idxs <- groupMatch(gsub("\\.","-",names(f)), rownames(cin))
  f.m <- f[idxs[[1]]]
  cin.m <- cin[idxs[[2]],]
  R <- subtypeMutAssoc <- sapply(levels(f.m), function(lvl){
    lvl.f <- factor(f.m == lvl)
    wilcox.test(log(cin.m) ~ lvl.f)$p.value
  })
  return(R) 
}

testMutAssoc <- function(pMatrix, mutPatTbl){
  f <- pmatrixToFactor(pMatrix)
  idxs <- groupMatch(gsub("\\.","-",names(f)), colnames(mutPatTbl))
  f.m <- f[idxs[[1]]]
  mutPatTbl.m <- mutPatTbl[, idxs[[2]]]
  
  R <- subtypeMutAssoc <- sapply(levels(f.m), function(lvl){
    lvl.f <- f.m == lvl
    apply(mutPatTbl.m, 1, function(x){
      fisher.test(lvl.f, x,alternative="greater")$p.value
    })
  })
  return(R)
}

plotCinAssoc <- function(R){
  df <- data.frame(subtype=names(R),pvals=-log10(R))
  sp <- ggplot(df, aes(x=subtype, y=pvals,fill=factor(subtype))) + 
    geom_bar(stat="identity") + 
    theme(panel.grid.minor.x=element_blank(), 
          panel.grid.major.x=element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1))
  sp <- sp + geom_hline(yintercept=.001, colour="red",linetype="dashed") + ylab("Significance")
  print(sp)
}

plotMutAssoc <- function(R,maxPerSubtype=10){
  tmp <- apply(R, 2, function(x){
    pvals <- sort(x)[1:maxPerSubtype]
    tmp <- -log10(pvals[1:maxPerSubtype])
    return(data.frame(gene=names(tmp),pvals=tmp))
  })
  
  mx <- max(sapply(tmp, function(x)max(x$pvals)))
  
  sp <- lapply(names(tmp), function(subtype) {
    df <- tmp[[subtype]]
    df$gene <- factor(df$gene, levels=df$gene[order(df$pvals, decreasing=TRUE)])
    sp <- ggplot(df, aes(x=gene, y=pvals)) + 
      geom_bar(stat="identity") + 
      ylim(0,mx) + 
      ggtitle(subtype) + 
      theme(panel.grid.minor.x=element_blank(), 
            panel.grid.major.x=element_blank(),
            axis.title.x = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1))
    return(sp)
  })
  multiplot(plotlist=sp, cols=3)
}
