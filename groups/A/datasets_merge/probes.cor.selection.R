# 
# Author: p angelino
###############################################################################

# Set work dir
## ---- EnvSetup ---------------------

if (!exists("work.dir")) {
	work.dir <- "/export/scratch/paolo/SageCRCsubtyping/analysis/work/"
}

# Set work dir
setwd(work.dir)
out.dir <- work.dir

###############################################################################
# Functions declaration
###############################################################################
matrix.sd.reorder <- function(M){
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
	return(M)
}

check_multiple_probes <- function(symbol, annotations){
	probes <- annotations[names(annotations)==symbol];
	if (length(probes)>1) {
		multiple <- TRUE
	} else {
		multiple <- FALSE
	}
	probes.list <- list(multiple = multiple, probes = probes)
	return(probes.list)
}


match.probesets <- function(common.features, ref, test){	
	# Compute correlation between probeset in two datasets
	# Update the gene matrices with the largest correlation probes
	# M: by gene matrix
	# X: by probe matrix
	for (x.gene in common.features) {
		print(x.gene)
		ref.probes.list <- check_multiple_probes(x.gene, ref$probes.ID)
		X.probes.list <- check_multiple_probes(x.gene, test$probes.ID)
		ref.cor <- NULL
		for (x.probe in ref.probes.list$probes) {
			ref.cor <- rbind(ref.cor, cor(ref$X[,x.probe], ref$M))
		}
		X.cor <- NULL
		for (x.probe in X.probes.list$probes) {
			X.cor <- rbind(X.cor, cor(test$X[,x.probe], test$M))
		}		
		print(dim(ref.cor))
		print(dim(X.cor))
		cor.cor <- cor(t(ref.cor),t(X.cor))
		max.ind <- arrayInd(which.max(cor.cor), dim(cor.cor))
		test$M[,x.gene] <- test$X[,X.probes.list$probes[max.ind[2]]];
		ref$M[,x.gene] <- ref$X[,ref.probes.list$probes[max.ind[1]]];
	}
}

matching.probesets   <- function(common.features, ref, test){
	# Compute correlation between probeset in two datasets
	# returns a list of correlation indexes between probes, by gene
	# M: by gene matrix
	# X: by probe matrix
#TODO add probe names to output list	
	gene.list <- list()
	for (x.gene in common.features) {
		print(x.gene)
		ref.probes.list <- check_multiple_probes(x.gene, ref$probes.ID)
		X.probes.list <- check_multiple_probes(x.gene, test$probes.ID)
		ref.cor <- NULL
		for (x.probe in ref.probes.list$probes) {
			ref.cor <- rbind(ref.cor, cor(ref$X[,x.probe], ref$M, use="pairwise.complete.obs"))
		}
		X.cor <- NULL
		for (x.probe in X.probes.list$probes) {
			X.cor <- rbind(X.cor, cor(test$X[,x.probe], test$M, use="pairwise.complete.obs"))
		}
		print(dim(ref.cor))
		print(dim(X.cor))
		#cat("ref: ",which(is.na(ref.cor)),"\n")
		#cat("test: ", which(is.na(X.cor)),"\n")
		gene.list[[x.gene]] <- cor(t(ref.cor),t(X.cor), use="pairwise.complete.obs")
		#print(gene.list[[x.gene]])
	}
	return(gene.list)
}

build.xgene.matrix.test <- function(Mp, Mx, selected.probes, gene.list, probes.id){
	# modified build.xgene.matrix function
	# it catches NA in the data an tries to find alternative probes
	#Mx <- matrix(nrow = nrow(Mp), ncol = length(gene.list))
	probe.names <- attr(Mx, "probeset")
	probe.cor <- vector(mode = "numeric", length = length(probe.names))
	names(probe.cor) <- names(probe.names)
	for (x.gene in gene.list) {
		probes.list <- check_multiple_probes(x.gene, probes.id)
		if (all(selected.probes[[x.gene]] %in% NA)){
			cat("warning: selected probes are NA for gene", x.gene,"\n")
			cat(selected.probes[[x.gene]],"\n")
			cat("number of probes:", length(probes.list$probes),"\n",
					c(probes.list$probes),"\n")
			if(length(c(probes.list$probes))>1){
				probes.sd <- apply(Mp[,c(probes.list$probes)], 2, function(z){sd(z,na.rm=TRUE)})
				print(probes.sd)
				print("removing NAs")
				probes.sd <- probes.sd[which(!is.na(probes.sd))]
				this.probe <- names(probes.sd)[probes.sd == max(probes.sd)]			
			} else {
				this.probe <- probes.list$probes				
			}
			print(this.probe)
		} else {
			this.probe <- probes.list$probes[selected.probes[[x.gene]]]
		}
		if(any(Mp[,this.probe] %in% NA)){
			#cat("NA values found in probe", this.probe, "of gene", x.gene,"\n")
			if(length(c(probes.list$probes))>1){
				probes.NAs <- apply(Mp[,c(probes.list$probes)], 2, function(z){any(z %in% NA)})
				cat(probes.NAs,which(probes.NAs == FALSE),"\n")
				probes.not.NAs <- which(probes.NAs == FALSE) 
				if(length(probes.not.NAs) != 0) {
					if(length(probes.not.NAs) > 1) {
						probes.sd <- apply(Mp[,c(probes.list$probes[probes.not.NAs])], 2, function(z){sd(z,na.rm=TRUE)})
						print(probes.sd)
						this.probe <- names(probes.sd)[probes.sd == max(probes.sd)]			
					} else {
						this.probe <- probes.list$probes[which(probes.NAs == FALSE)]
					}
				}
			} else {
				#cat("No alternative probe available\n")
			}
		}	
#		print(c(i,this.probe))
		Mx[,x.gene] <- Mp[,this.probe]
		probe.names[x.gene] <- this.probe 
		probe.cor[x.gene] <- attr(selected.probes[[x.gene]],"cor")
	}
	
	#colnames(Mx) <- gene.list
	#rownames(Mx) <- rownames(Mp)
	attr(Mx, "probeset") <- probe.names[colnames(Mx)]
	attr(Mx, "cor") <- probe.cor[colnames(Mx)]
	return(Mx)		
}

select.test.probes <- function(x, test, selected.probes){
	# x: list of common genes
	# test: environment for dataset
	# selected.probes: selected probes in the ref dataset
	x.names <- names(x)
	selected.probes.names <- names(selected.probes)
	test$selected.probes <- list()
	for (x.gene in x.names) {
		i <- x[[x.gene]]
		if (x.gene %in% selected.probes.names) {
			if (dim(i)[2] > 1) {
				x.cor = i[selected.probes[[x.gene]],]	
				# Select first maximum if there are more maxima
				test$selected.probes[x.gene] <- seq(along=x.cor)[x.cor == max(x.cor)][1]		
			} else {
				test$selected.probes[x.gene] <- 1 
			}	
		} else {
			max.id <- arrayInd(which.max(i),dim(i))
			test$selected.probes[x.gene] <- max.id[2]
		}
		attr(test$selected.probes[[x.gene]], "cor") <- i[test$selected.probes[[x.gene]]]
	}
}

###############################################################################
###############################################################################

# environment
library(nclust)
library(GEOquery)
library(org.Hs.eg.db)
library(ggplot2)
library(grema)


## ---- ProbesCor ---------------------

if (!exists("common.features")) {
	load("common.features.all.Rdata")	
}

load("dataset.list.Rdata")
# store dataset size
dataset.size <- list()
for (n in dataset.list) {
	message(paste("Now processing ", n, "..."))	
	test <- new.env()
	load(paste0(n,".Rdata"), envir=test)
	dataset.size[[n]] <- nrow(test$M)
}

reference.dataset <- "french"
dataset.list <- dataset.list[which(dataset.list!=reference.dataset)]
sub.selection <- common.features[1:5005]

# load reference dataset
ref <- new.env()
load(paste0(reference.dataset,".Rdata"), envir=ref)
ref$MM <- ref$M
		
require(foreach)
require(doMC)
registerDoMC(cores=6)

doparfunc <- function(n, ref, sub.selection){
	test <- new.env()
	load(paste0(n,".Rdata"), envir=test)
	common.features <- intersect(colnames(ref$MM),colnames(test$M))
	test$M[which(is.na(test$M))] <- 0
	test$X[which(is.na(test$M))] <- 0
	common.features <- common.features[!(common.features %in% "")]
	ref$M <- ref$MM[,intersect(sub.selection,common.features)]
	test$M <- test$M[,intersect(sub.selection,common.features)]
	probes.cor <- matching.probesets(common.features = common.features,ref,test)
	rm(test)
	return(probes.cor)
}

probes.cor.results <- foreach (n=dataset.list) %dopar% doparfunc(n, ref, sub.selection =  sub.selection)
names(probes.cor.results) <- dataset.list

rm(ref)

features <- list()
for (n in dataset.list) {
	message(paste("Now processing ", n, "..."))	
	features[[n]] <- length(probes.cor.results[[n]]) 
}

save(file="probes.cor.Rdata", probes.cor.results, features, dataset.size)

## ---- SelectFrenchProbes ---------------------

selected.probes <- list()
for (x.gene in common.features) {
	ref.probes.nbr <- dim(probes.cor.results[[1]][[x.gene]])[1]
	if (ref.probes.nbr > 1) {
		probe.ranking <- 0
		for (n in dataset.list) {
			max.id <- apply(probes.cor.results[[n]][[x.gene]],1, function(x){ max(x) })
			# each dataset vote for probe selection, vote is weighted by dataset size
			probe.ranking <- probe.ranking + rank(max.id)*dataset.size[[n]] 
		}
		tmp <- seq(along=probe.ranking)[probe.ranking == max(probe.ranking)]
		selected.probes[[x.gene]] <- tmp[1]
	} else {
		selected.probes[[x.gene]] <- 1 		
	}
	attr(selected.probes[[x.gene]], "cor") <- 1
}

#save(file="selected.probes2.Rdata", selected.probes)
save(file="selected.probes.Rdata", selected.probes)

ref <- new.env()
load(paste0(reference.dataset,".Rdata"), envir=ref)
Mgenes <- build.xgene.matrix.test(ref$X, ref$M, selected.probes, common.features, ref$probes.ID)
save(file=paste0(reference.dataset,".selected.Rdata"), Mgenes)
rm(ref, Mgenes)


## ---- SelectOtherProbes ---------------------
# in other datasets
#
#
foreach (n=dataset.list) %dopar% {
#for (n in dataset.list) {
	message(paste("Now processing ", n, "..."))	
	test <- new.env()
	load(paste0(n,".Rdata"), envir=test)
	select.test.probes(probes.cor.results[[n]], test, selected.probes)
	Mgenes <- build.xgene.matrix.test(test$X, test$M, test$selected.probes, 
			intersect(common.features,colnames(test$M)), test$probes.ID)
	save(file=paste0(n,".selected.Rdata"), Mgenes)
	rm(test, Mgenes)
}


