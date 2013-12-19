
plotByFeature <- function(df, phenotypeName){

  data <- df[df$phenotype==phenotypeName,]
  medData <- data[data$quant=="med",]
  sp <- ggplot(medData, aes(x=subtype, y=-log10(value))) + 
    geom_bar(stat="identity",aes(fill = factor(subtype),position="dodge"))
  sp <- sp + facet_wrap(~dataset,ncol=4)
}

plotSubtypeDistribution <- function(groupResult){
  lvls <- lvls4Group(groupResult)
  R <- sapply(groupResult, function(pmatrix){
    tmp <- table(pmatrixToFactor(pmatrix,lvls))
    return(tmp/sum(tmp))
  })
  df <- melt(R)
  print(ggplot(df, aes(x=X1,y=value, fill=X1)) + geom_bar(stat="identity",position="dodge") +
    facet_wrap(~X2,ncol=8))
}

plotDatasetsUsedByGroup <- function(groupResults){
  df <- melt(lapply(groupResults, function(gr){
    tmp <- melt(sapply(gr, function(ds){ nrow(ds)}))
    tmp$dataset <- rownames(tmp)
    tmp
  }))
  print(ggplot(df, aes(x=L1,y=value, fill=L1)) + geom_bar(stat="identity",position="dodge") +
          facet_wrap(~dataset,ncol=8))
}
