# TODO: Add comment
# 
# Author: pangelin
###############################################################################

###############################################################################
# Functions for computing MSI signature
###############################################################################

dist.cos = function(x1, x2)
{	
	d = sum(x1 * x2) / sqrt(sum(x1^2) * sum(x2^2))
	
	return (d)
}


## X - a gene expression matrix, genes by columns
## msi.sig - signature, centroids by columns
score.msi.cc = function(X, msi.sig)
{
	s = matrix(0, nrow=nrow(X), ncol=2)
	colnames(s) = colnames(msi.sig)
	
	s[,1] = apply(X, 1, function(z) (dist.cos(z, msi.sig[,1])) )
	s[,2] = apply(X, 1, function(z) (dist.cos(z, msi.sig[,2])) )
	
	return (s)
}

## X - a gene expression matrix, genes by columns
## msi.sig - signature, centroids by columns
score.msi.cor = function(X, msi.sig, method='pearson')
{
	s = cor(t(X), msi.sig, method=method)
	
	return (s)
}

## X - a gene expression matrix, genes by columns
## msi.sig - signature, centroids by columns
score.msi.euclid = function(X, msi.sig)
{
	s = matrix(0, nrow=nrow(X), ncol=2)
	
	s[,1] = apply(X, 1, function(z) (sqrt(sum((z - msi.sig[,1])^2))) )
	s[,2] = apply(X, 1, function(z) (sqrt(sum((z - msi.sig[,2])^2))) )
	
	return (s)
}

## Given distances to two 'centroids' combine them in a
## single score such that higher scores are associated
## with the 1st centroid.
dist.combiner = function(d, method=c('delta'))
{
	if (match(method, c('delta'))) {
		return (d[,1] - d[,2])
	}
	
	return (NULL)
}


## Finds a threshold such that MCC is maximized.
optimize.thr = function(s, y)
{
	z = seq(min(s), max(s), length.out=1000)
	m = NULL
	for (t in z) {
		p = (s >= t) + 0   # predictions
		r = perf.matrix(y, p)
		m = c(m, perf.mcc(r))
	}
	
	return (z[which.max(m)])
}


## Gene symbols to Entrez Gene IDs
symbolMap <- function(ids) {
	sapply(AnnotationDbi::mget(as.character(ids), org.Hs.egSYMBOL2EG,
					ifnotfound = NA),
			function(x) paste(x, collapse = "//"))
}



###############################################################################
###############################################################################

require(org.Hs.eg.db)

if (!exists("work.dir")) {
	work.dir <- "/export/scratch/paolo/SageCRCsubtyping/analysis/work/"
}

# Set work dir
setwd(work.dir)
out.dir <- work.dir

datapath <- "/export/scratch/paolo/SageCRCsubtyping/annotations"

S = read.csv(paste(datapath, "msi64gene.tab", sep="/"), header=TRUE, row.names=2, as.is=TRUE, sep='\t')[,-c(1,2)]

#source(paste(scriptpath, "agendia.R", sep=""))

S.entrezID <- symbolMap(rownames(S))
S <- S[which(S.entrezID!='NA'),]
rownames(S) <- S.entrezID[which(S.entrezID!='NA')]

load("dataset.list.Rdata")

require(foreach)
require(doMC)
registerDoMC(cores=6)

parfunc <- function(n, S1) {
	load(paste0(n,".selected.Rdata"))
	
	idx = intersect(rownames(S1), colnames(Mgenes))
	S1 = S1[idx,]
	Z = Mgenes[,idx]
	
	d.cos.TA = score.msi.cc(Z, S1)
	s.cos.TA = dist.combiner(d.cos.TA)
	
	rm(Mgenes)
	return(s.cos.TA)
}

msi.scores.list <- foreach (n=dataset.list) %dopar% parfunc(n, S)

names(msi.scores.list) <- dataset.list

save(msi.scores.list, file=paste(datapath, "MSIscores.Rdata", sep=""))


## Normalize with ComBat

if (!exists("common.features")) {
	load("common.features.all.Rdata")	
}

m.ALL <- matrix(nrow=0, ncol=length(common.features))
datasets <- vector()
cnt <- 1
for (n in dataset.list) {
	message(paste("Now processing ", n, "..."))	
	load(paste0(n,".selected.Rdata"))
	m.ALL <- rbind(m.ALL, Mgenes[,common.features])
	datasets <- c(datasets, rep(cnt,nrow(Mgenes)))
	cnt <- cnt + 1
}
rm(Mgenes)

require(sva)
require(imputation)
if (any(is.na(m.ALL))) {
	cat("Percentage on NA in m.ALL:", length(which(is.na(m.ALL)))/dim(m.ALL)[1]/dim(m.ALL)[2]*100,"\n")
	X.impute <- meanImpute(m.ALL)
	X.ALL <- X.impute$x
	cat("Percentage on NA in X.ALL:", length(which(is.na(X.ALL)))/dim(X.ALL)[1]/dim(X.ALL)[2]*100,"\n")
}

m.COM <- ComBat(t(X.ALL), datasets, mod=unlist(msi.scores.list), numCovs=1, par.prior=TRUE, prior.plots=TRUE)
m.COM <- t(m.COM)

## Split datasets matrices

cnt <- 1
M.sets <- list()
for (n in dataset.list) {
	message(paste("Now processing ", n, "..."))	
	Mgenes <- m.COM[which(datasets==cnt),]
	M.sets[[n]] <- Mgenes
	save(file=paste0(n,".selected.norm.Rdata"), Mgenes)
	cnt <- cnt + 1	
}
rm(Mgenes)

