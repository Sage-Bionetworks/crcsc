# 
# Author: p angelino
###############################################################################

if (!exists("work.dir")) {
	work.dir <- "/export/scratch/paolo/SageCRCsubtyping/analysis/work/"
}

# Set work dir
setwd(work.dir)
out.dir <- work.dir

require(imputation)

###############################################################################
# matrix reordering
###############################################################################

matrix.sd.reorder <- function(M, probes.ID){
	# reorder the matrix according to gene std
	# and select unique gene symbols
	features <- colnames(M)
	M.stds <- apply(M,2, function(z){
				mad(z,na.rm=TRUE)
			})
	ordering <- order(M.stds,decreasing=TRUE)
	#M.stds <- M.stds[ordering]
	M <- M[,ordering]
	#colnames(M) <- features[ordering]
	
	# extract subset with unique gene names
	M <- M[, !duplicated(features[ordering], fromLast=FALSE)] 
	probes.ID <- probes.ID[ordering]
	probes.ID <- probes.ID[!duplicated(features[ordering], fromLast=FALSE)] 
	attr(M, "probeset") <- probes.ID
	return(M)
}


# Prepare matrix annotation with EntrezIds 
# Prepare reference gene matrix with mad 
for (n in names(exprList)) {
	x = exprList[[n]]
	message(paste("Now processing ", n, "..."))
	# Keep only rows corresponding to Entrez ID
	entrezID <- x$entrezMapper(rownames(x$data.mat))
	keep.idx <- which(!entrezID=="NA")
	x$data.mat <- x$data.mat[keep.idx, ]
	entrezID <- entrezID[keep.idx]
	X <- t(x$data.mat)
	M <- t(x$data.mat)
	colnames(M) <- entrezID
	# probe / gene mapping
	probes.ID <- rownames(x$data.mat)
	names(probes.ID) <- entrezID
	M <- matrix.sd.reorder(M, probes.ID)
	str(M)
	str(X)
	if (any(is.na(M))) {
		cat("Percentage on NA in M:", length(which(is.na(M)))/dim(M)[1]/dim(M)[2]*100,"\n")
		cat("Percentage on NA in X:", length(which(is.na(X)))/dim(X)[1]/dim(X)[2]*100,"\n")
		X.impute <- meanImpute(X)
		X <- X.impute$x
		M <- X
		colnames(M) <- names(probes.ID)
		M <- matrix.sd.reorder(M, probes.ID)
		str(M)
		str(X)
	}
	str(probes.ID)
	save(file=paste0(out.dir,n,".Rdata"), M, X, probes.ID)
	rm(M, X, probes.ID)
}

# Load dataset list
if (!exists("dataset.list")) {
	load("dataset.list.Rdata")
}

# Find common genes
## ---- Find common genes ---------------------
gene.lists <- list()
for (n in names(exprList)) {
	x = exprList[[n]]
	message(paste("Now processing ", n, "..."))
	entrezID <- x$entrezMapper(rownames(x$data.mat))
	keep.idx <- which(!entrezID=="NA")
	entrezID <- entrezID[keep.idx]
	gene.lists[[n]] <- entrezID
	print(paste("genes in",n,"=",length(gene.lists[[n]])))
}

n <- names(exprList)[1]
tmp <- unique(gene.lists[[n]])
print(paste("genes in",n,"=",length(tmp)))
common.features <- tmp
for (n in names(exprList)[2:length(exprList)]) {
	tmp <- unique(gene.lists[[n]])
	print(paste("genes in",n,"=",length(tmp)))
	common.features <- intersect(common.features, tmp)
	print(paste("common features:", length(common.features)))	
}

save(file="common.features.all.Rdata", common.features)