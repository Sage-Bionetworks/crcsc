compute.cor <- FALSE

## ---- GeneSubselectQTL ---------------------

gene.cor.cor <- function(ref.X, ref.M, test.X, test.M){
	# Compute correlation between genes in two datasets
	# returns a vector of correlation indexes by gene
	# M: by reference matrix
	# X: by gene matrix
	gene.cor <- vector()
	cnt <- 1
	for (x.gene in colnames(test.X)) {
		cat("[",cnt,"] ", x.gene, "\n")
		ref.cor <- cor(ref.X[,x.gene], ref.M, use="pairwise.complete.obs")
		X.cor <- cor(test.X[,x.gene], test.M, use="pairwise.complete.obs")
		gene.cor[x.gene] <- cor(t(ref.cor),t(X.cor), use="pairwise.complete.obs")
		cnt <- cnt + 1
	}
	return(gene.cor)
}

require(foreach)
require(doMC)
registerDoMC(cores=6)

if (!exists("work.dir")) {
	work.dir <- "/export/scratch/paolo/SageCRCsubtyping/analysis/work/"
}

# Set work dir
setwd(work.dir)
out.dir <- work.dir

# quantiles
qtl.low <- 0.05
qtl.high <- 0.95

# Load dataset list
if (!exists("dataset.list")) {
	load("dataset.list.Rdata")
}

# Load datasets
if (!exists("M.sets")) {
	message(paste("Loading datasets for gene filtering..."))	
	M.sets <- list()
	for (n in dataset.list) {
		message(paste("Now processing ", n, "..."))	
		load(paste0(n,".selected.norm.Rdata"))
		M.sets[[n]] <- Mgenes
	}
}


# Load commmon genes
if (!exists("common.features.all")) {
	load("common.features.all.Rdata")
}

# compute IQR and MEDIAN
iqr.f <- function(z){
	dd <- quantile(z,probs=c(qtl.high, qtl.low, 0.5), na.rm=TRUE)
	out <- vector()
	out[1] <- (dd[1]-dd[2])/dd[3] #IQR/median
	out[2] <- dd[1]-dd[2]	#IQR
	out[3] <- dd[3]			#MEDIAN
	return(out)
}

M.sd <- foreach (n=dataset.list) %dopar% apply(M.sets[[n]],2, iqr.f)
names(M.sd) <- dataset.list


#interquartile range
IQR <- vector("list",length(M.sd))
names(IQR) <- names(M.sd)
#interquartile range / median
IQRn <- vector("list",length(M.sd))
names(IQRn) <- names(M.sd)
#median
MEDIAN <- vector("list",length(M.sd))
names(MEDIAN) <- names(M.sd)
#intequartile range / median of(IRQ)
IQRrn <- vector("list",length(M.sd))
names(IQRrn) <- names(M.sd)
for (dset in names(M.sd)) {
	IQRn[[dset]] <- M.sd[[dset]][1,]
	IQR[[dset]] <- M.sd[[dset]][2,]
	MEDIAN[[dset]] <- M.sd[[dset]][3,]
}
iqr.med <- lapply(IQR, median, na.rm=TRUE)
for (dset in names(M.sd)) {
	IQRrn[[dset]] <- IQR[[dset]]/iqr.med[[dset]]
}	

# Correlations
if (!file.exists("cor_iqr.Rdata")){
	compute.cor <- TRUE
}
if (compute.cor) {	
	ref <- new.env()

	reference.dataset <- "french"
	dataset.rlist <- dataset.list[which(dataset.list!=reference.dataset)]
	sub.selection <- common.features[1:5005]
	load(paste0(reference.dataset,".Rdata"), envir=ref)

	doparfunc <- function(n, ref){
		test <- new.env()
		load(paste0(n,".Rdata"), envir=test)
		corcor <- gene.cor.cor(M.sets[[reference.dataset]],ref$M[,sub.selection],M.sets[[n]],test$M[,sub.selection])
		rm(test)
		return(corcor)
	}

	COR <- foreach (n=dataset.rlist) %dopar% doparfunc(n, ref)
	names(COR) <- dataset.rlist

	save(file="cor_iqr.Rdata", COR, M.sd)
} else {
	load("cor_iqr.Rdata")
}


# Compute dataset weight
dset.w <- vector("numeric",length(names(M.sets)))
names(dset.w) <- names(M.sets)
for (n in names(M.sets)) {
	dset.w[n] <- nrow(M.sets[[n]])	 
}
tot.size <- sum(dset.w)
dset.w <- dset.w/tot.size

# Thresholds
iqr.threshold.list <- c(0.5, 0.75, 1)
iqr.threshold.names <- c("05","075","1")
cor.threshold.list <- c(0.3, 0.4, 0.5, 0.6)
#iqr.threshold.list <- c(0.5)
#iqr.threshold.names <- c("05")
#cor.threshold.list <- c(0.4)

# defining variable to be used in the loops
cor.output <- list()
coriqr.output <- list()
output.names <- vector("numeric", length(cor.threshold.list*length(iqr.threshold.list)))
summary.output <- list()

# gene filtering loop
cnt <- 1
pdf(file="voting.hist.pdf")
iqr.done <- NULL
cor.done <- NULL
for (j in seq(along = iqr.threshold.list)) {
	# Loop on iqr thresholds
	iqr.threshold <- iqr.threshold.list[j]
	prefix = iqr.threshold.names[j]
	desc2 = paste0("iqr threshold ",iqr.threshold)	
	fid <- file(paste0("iqr",prefix,".summary.out"),"w") 
	
	for (i in seq(along = cor.threshold.list)) {
		# Loop on cor thresholds
		cor.threshold <- cor.threshold.list[i]
		cat("IQR threshold:", iqr.threshold, "\n", file=fid)
		# apply IQR threshold
		iqr.vote <- data.frame(matrix(0,length(common.features),length(names(COR))+1),row.names=common.features)
		names(iqr.vote) <- c(names(COR),"french")
		for (dset in c(names(COR),"french")) {
			iqr.vote[[dset]][IQRrn[[dset]] >= iqr.threshold] <- dset.w[dset]
			cat("IQR. Dataset: ",dset,", number of matching genes: ", length(which(iqr.vote[[dset]]==dset.w[dset])),"\n", file=fid)
		}
		iqr.vote$sum <- rowSums(iqr.vote)
		if (!(iqr.threshold %in% iqr.done)) {
			hist(iqr.vote$sum, main=paste("iqr",iqr.threshold))
			iqr.done <- c(iqr.done, iqr.threshold)
		}
		cat("Correlation threshold:", cor.threshold, "\n", file=fid)
# apply correlation threshold
		cor.vote <- data.frame(matrix(0,length(common.features),length(names(COR))),row.names=common.features)
		names(cor.vote) <- names(COR)
		for (dset in names(COR)) {
			cor.vote[[dset]][COR[[dset]] >= cor.threshold] <- dset.w[dset]
			cat("COR. Dataset: ",dset,", number of matching genes: ", length(which(cor.vote[[dset]]==dset.w[dset])),"\n", file=fid)
		}
		cor.vote$sum <- rowSums(cor.vote)
		if (!(cor.threshold %in% cor.done)) {
			hist(cor.vote$sum,main=paste("cor", cor.threshold))
			cor.done <- c(cor.done, cor.threshold)
		}
# select gene across datasets
		iqr.selection <- rownames(iqr.vote)[iqr.vote$sum>=0.9]
		cor.selection <- rownames(cor.vote)[cor.vote$sum>=0.5]
		sub.selection <- intersect(iqr.selection, cor.selection)
		gene.list <- sub.selection
		save(file=paste0("gene.filtered.",iqr.threshold,".",cor.threshold,".Rdata"),gene.list)
		cat("Number of common selected gene (COR):", length(cor.selection),"\n", file=fid)
		cat("Number of common selected gene (IQR):", length(iqr.selection),"\n", file=fid)
		cat("Number of common selected gene (IQR&&COR):", length(sub.selection),"\n", file=fid)
		summary.output[[cnt]] <- c(length(cor.selection),length(iqr.selection),length(sub.selection))
		cat("********************************************\n\n", file=fid)		
		cnt <- cnt + 1			
	}
	close(fid)
}

X11()
hist(cor.vote$sum)
dev.off()

## ---- GeneSubselectPlot ---------------------

require(ggplot2)
require(gridExtra)
pdf(file="gene.filter.cor.iqr.pdf")
for (dset in names(COR)) {
	tmp.cor <- unlist(COR[dset])
	tmp.iqr <- unlist(IQR[dset])
	hist_top <- ggplot()+geom_histogram(aes(tmp.cor)) + xlab("COR")
	empty <- ggplot()+geom_point(aes(1,1), colour="white")+
			theme(axis.ticks=element_blank(), 
					panel.background=element_blank(), 
					axis.text.x=element_blank(), axis.text.y=element_blank(),           
					axis.title.x=element_blank(), axis.title.y=element_blank())
	
	scatter <- ggplot()+geom_point(aes(tmp.cor, tmp.iqr)) +
			xlab("COR") + ylab("IQR/median") + ggtitle(dset)
	hist_right <- ggplot()+geom_histogram(aes(tmp.iqr))+coord_flip() + xlab("IQR/median")
	
	grid.arrange(hist_top, empty, scatter, hist_right, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))
	#readline()
}
dev.off()

## ---- GeneSubselectPlot2 ---------------------
pdf(file="gene.filter.cor.iqrn.pdf")
for (dset in names(COR)) {
	tmp.cor <- unlist(COR[dset])
	tmp.iqr <- unlist(IQRrn[dset])
	hist_top <- ggplot()+geom_histogram(aes(tmp.cor)) + xlab("COR")
	empty <- ggplot()+geom_point(aes(1,1), colour="white")+
			theme(axis.ticks=element_blank(), 
					panel.background=element_blank(), 
					axis.text.x=element_blank(), axis.text.y=element_blank(),           
					axis.title.x=element_blank(), axis.title.y=element_blank())
	
	scatter <- ggplot()+geom_point(aes(tmp.cor, tmp.iqr)) +
			xlab("COR") + ylab("IQR") + ggtitle(dset)
	hist_right <- ggplot()+geom_histogram(aes(tmp.iqr))+coord_flip() + xlab("IQR")
	
	grid.arrange(hist_top, empty, scatter, hist_right, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))
	#readline()
}

for (dset in names(COR)) {
	tmp.cor <- unlist(MEDIAN[dset])
	tmp.iqr <- unlist(IQR[dset])		
	median.iqr.plot <- ggplot()+geom_point(aes(tmp.cor, tmp.iqr)) +
			xlab("MEDIAN") + ylab("IQR") + ggtitle(dset)
	print(median.iqr.plot)
}
dev.off()
