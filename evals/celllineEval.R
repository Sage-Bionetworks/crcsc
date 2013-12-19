

testDrugAssoc <- function(pMatrix, drugTbl, sample=0){
  idxs <- groupMatch(clean.names(rownames(drugTbl)),clean.names(rownames(pMatrix)))
  if(length(idxs[[1]]) == 0){
    cat("No match between phenotype and pMatrix\n")
    return(NULL)
  }
  drugTbl.m <- drugTbl[idxs[[1]],]
  pMatrix.m <- pMatrix[idxs[[2]],]
  f.m <- pmatrixToFactor(pMatrix.m)
  
  # for each drug, test association
  R <- apply(drugTbl.m, 2, function(drugVec){
    na.mask <- !is.na(drugVec)
    subtypeMutAssoc <- sapply(levels(f.m), function(lvl){
      lvl.f <- factor(f.m == lvl)
      wilcox.test(drugVec[na.mask] ~ lvl.f[na.mask])$p.value
    })
  })
  return(R)
}