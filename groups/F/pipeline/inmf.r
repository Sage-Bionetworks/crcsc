# Functions for performing the iNMF subtyping
# 
# Author: Andreas Schlicker (a.schlicker@nki.nl)
###############################################################################

##' Map iNMF signatures to identifiers used in another data set.
##' @param which one of "step1", "step2.c1", "step2.c2"
##' @param type one of "ps", "symbol". "ps" uses Affy HGU133+2 signatures as 
##' starting point. "symbol" uses signatures that were already mapped to 
##' gene symbols. 
##' @param mapping a matrix mapping data set-specific IDs (first column) to either Affy 
##' HGU133+2 probe sets or gene symbols (second column). If NULL (default), the 
##' signature IDs will be kept.
##' @return a named list with the signatures
##' @author Andreas Schlicker
signatureList = function(which=c("step1", "step2.c1", "step2.c2"), type=c("ps", "symbol", "entrez"), mapping=NULL) {
	which = match.arg(which)
	type = match.arg(type)
	
	res = list()
	if (type == "ps") {
		if (which == "step1") {
			res = list("1"=inmf.clust1.ps, "2"=inmf.clust2.ps)
		} else if (which == "step2.c1") {
			res = list("1.1"=inmf.clust1.1.ps, "1.2"=inmf.clust1.2.ps, "1.3"=inmf.clust1.3.ps)
		} else if (which == "step2.c2") {
			res = list("2.1"=inmf.clust2.1.ps, "2.2"=inmf.clust2.2.ps)
		}
	} else if (type == "symbol") {
		if (which == "step1") {
			res = list("1"=inmf.clust1.sym, "2"=inmf.clust2.sym)
		} else if (which == "step2.c1") {
			res = list("1.1"=inmf.clust1.1.sym, "1.2"=inmf.clust1.2.sym, "1.3"=inmf.clust1.3.sym)
		} else if (which == "step2.c2") {
			res = list("2.1"=inmf.clust2.1.sym, "2.2"=inmf.clust2.2.sym)
		}
	} else if (type == "entrez") {
		if (which == "step1") {
			res = list("1"=inmf.clust1.entrez, "2"=inmf.clust2.entrez)
		} else if (which == "step2.c1") {
			res = list("1.1"=inmf.clust1.1.entrez, "1.2"=inmf.clust1.2.entrez, "1.3"=inmf.clust1.3.entrez)
		} else if (which == "step2.c2") {
			res = list("2.1"=inmf.clust2.1.entrez, "2.2"=inmf.clust2.2.entrez)
		}
	}
	
	if (!is.null(mapping)) {
		return(lapply(res, function(x) { unique(mapping[which(mapping[, 2] %in% x), 1])}))
	}
	
	res
}

##' Removes identifiers which are not contained in the allIds vector.
##' @param signatures list of signatures that are to be filtered
##' @param allIds vector with IDs that should be retained
##' @return names list after filtering IDs
##' @author Andreas Schlicker
filterSignatures = function(signatures, allIds) {
	lapply(signatures, function(x) { intersect(x, allIds)})
}

##' Prepare a list containing all signatures for the iNMF subtyping. 
##' The return value can directly be used as input for the iNMF function.
##' @param allIds vector containing all IDs for which expression data is available
##' @param type one of "ps", "symbol". "ps" uses Affy HGU133+2 signatures as 
##' starting point. "symbol" uses signatures that were already mapped to 
##' gene symbols.
##' @param mapping a matrix mapping data set-specific IDs (first column) to either Affy 
##' HGU133+2 probe sets or gene symbols (second column). If NULL (default), the 
##' signature IDs will be kept.
##' @return named list with all three signatures
##' @author Andreas Schlicker
signaturesList = function(allIds, type=c("ps", "symbol"), mapping=NULL) {
	list(step1= filterSignatures(signatureList("step1", type, mapping), allIds),
		 step2.c1 = filterSignatures(signatureList("step2.c1", type, mapping), allIds),
		 step2.c2 = filterSignatures(signatureList("step2.c2", type, mapping), allIds))
}

##' Map sample clusters to the correct subtype.  
##' Computes the average expression value for each sample cluster and each signature.
##' Then iteratively assigns subtype IDs to the clusters starting with the cluster
##' with the highest certainty for a given subtype. 
##' The function assumes that all samples and features are contained in the expression 
##' matrix. Missing values will be ignored.
##' @param samples named list with sample clusters
##' @param named list of signatures
##' @param exprs expresssion matrix with samples in columns and features in rows
##' @return the mapping between cluster ID and subtype label
##' @author Andreas Schlicker
assignClusterId = function(samples, signatures, exprs) {
	# Calculated mean expression for each sample cluster and each signature
	means = sapply(samples, function(samp) { sapply(signatures, function(x) { mean(exprs[x, samp], na.rm=TRUE) }) }, simplify=TRUE)
	# To each subtype, assign the cluster with the largest mean expression
	mapping = apply(means, 1, function(x) { names(x)[which(x == max(x))] })
	
	# Check whether one cluster was assigned more than once
	repeated = names(table(mapping))[table(mapping) > 1]
	if (length(repeated) > 0) {
		# Get the missing cluster ID
		missing = setdiff(names(samples), unique(mapping))
		# Get the submatrix with mean values
		sub.mat = means[names(mapping)[which(mapping == repeated)], c(repeated, missing)]
		# Clear the mapping with the repeated cluster
		mapping[names(mapping)[which(mapping == repeated)]] = NA
		# For each subtype without assignment, calculate the difference in mean expression
		diff = apply(sub.mat, 1, function(x) { abs(x[1] - x[2]) })
		# The subtype with the largest difference in mean expression
		temp = names(diff)[diff==max(diff)]
		# Assign the cluster with the largest mean expression to this subtype
		mapping[temp] = names(means[temp, ])[means[temp, ] == max(means[temp, ])]
		# Assign the missing cluster to the subtype without mapping
		mapping[is.na(mapping)] = setdiff(names(samples), mapping)
	}
	# Reverse the mapping between subtype and cluster ID
	res = names(mapping)
	names(res) = mapping
	
	res
	
#	# Calculate mean expression for each cluster and subtype signature
#	means = lapply(samples, function(samp) { sapply(signatures, function(x) { mean(exprs[x, samp], na.rm=TRUE) }) })
#	
#	# Calculate the margin for each cluster. The margin is defined as the difference between 
#	# highest average signature expression and second highest expression. The larger the margin
#	# is, the more certain are we that this is the right subtype.
#	margin = sapply(means, function(x) { sort(x, decreasing=TRUE)[1] - sort(x, decreasing=TRUE)[2]})
#	names(margin) = names(means)
#	
#	# Build the cluster to subtype ID mapping
#	mapping = rep(NA, times=length(signatures))
#	names(mapping) = names(signatures)
#	# Go through all clusters starting with the one with highest margin
#	for (n in names(sort(margin, decreasing=TRUE))) {
#		# Cycle through all subtypes in the order of the mean expression
#		for (subtype in names(sort(means[[n]], decreasing=TRUE))) {
#			# If we didn't assign that subtype label to any cluster yet, do it.
#			# If there is already a cluster ID, some other cluster had a higher certainty.
#			if (is.na(mapping[subtype])) {
#				mapping[subtype] = n
#				break
#			}
#		}
#	}
#	
#	# Have to reverse mapping of subtype label to cluster ID
#	res = names(mapping)
#	names(res) = mapping
#	
#	# And we return the mapping between cluster ID and subtype label
#	res
}

##' Performs one step of the iNMF clustering. Essentially, this is a 
##' hierarchical clustering of the given samples using the given signatures.
##' Cluster IDs are assigned using function "assignClusterId" but should be 
##' manually verified using a heatmap or similar. The number of clusters is
##' equal to the number of signatures. Also computes silhouette values for the 
##' clustering. 
##' @param exprs expression matrix with samples in columns and features in rows
##' @param signatures list of signatures to use.
##' @param samples vector with sample names. If NULL (default), all samples will be
##' clustered.
##' @param silhouette boolean value indicating whether the silhouette is to be 
##' computeted
##' @return Named list with two elements. The first element "clustering" contains a
##' named list with the assignment of samples to clusters. The names correspond to the 
##' names used for the signatures. The second element "silhouette" contains an object 
##' of class silhouette, if silhouette==TRUE.
##' @author Andreas Schlicker
subtype = function(exprs, signatures, samples=NULL, silhouette=TRUE) {
	if (silhouette) {
		require(cluster) || stop("I need package \"cluster\" to calculate the silhouette!")
	}
	
	if (is.null(samples)) {
		samples = colnames(exprs)
	}
	dist.mat = 1 - cor(exprs[unlist(signatures), samples], use="pairwise.complete.obs")
	
	clustering = cutree(hclust(as.dist(dist.mat), method="complete"), k=length(signatures))
	
	clust.list = list()
	for (i in unique(clustering)) {
		clust.list[[as.character(i)]] = names(clustering[clustering == i])	
	}
	
	idMap = assignClusterId(clust.list, signatures, exprs)
	
#	res = list()
#	idMap = list()
#	for (i in unique(clustering)) {
#		samps = names(clustering[clustering == i])
#		cId = assignClusterId(samps, signatures, exprs)
#		if (length(intersect(names(res), cId)) > 0) {
#			warning(paste("Cluster ID ", cId, " was assigned more than once!"))
#		}
#		res[[cId]] = samps
#		idMap[[as.character(i)]] = cId
#	}
	
	res = list()
	for (n in names(clust.list)) {
		res[[idMap[n]]] = clust.list[[n]]
	}
	
	sil = NULL
	if (silhouette) {
		sil = silhouette(clustering, dist.mat)
		for (i in names(idMap)) {
			sil[which(sil[, "cluster"] == i), "cluster"] = idMap[[i]]
			sil[which(sil[, "neighbor"] == i), "neighbor"] = idMap[[i]]
		}
		rownames(sil) = samples
	}
	
	if (is.null(sil)) {
		list(clustering=res)
	} else {
		list(clustering=res, silhouette=sil)
	}
}

##' Create the quality control heatmap. Column and row colors indicate cluster
##' membership for samples and features, respectively. 
##' @param exprs expression matrix with samples in columns and features in rows
##' @param clustering named list with sample to cluster assignment
##' @param signatures named list with feature to signature assignment
##' @return heatmap.2 for plotting
##' @author Andreas Schlicker
createHeatmap = function(exprs, clustering, signatures) {
	require(gplots) || stop("I need package \"gplots\" for doing this.")
	
	samples = c()
	for (n in sort(names(clustering))) {
		samples = c(samples, clustering[[n]])
	}
	
	features = c()
	for (n in sort(names(signatures))) {
		features = c(features, signatures[[n]])
	}
	
	bounds = quantile(exprs[unlist(signatures), unlist(clustering)], probs=c(0.05, 0.95), na.rm=TRUE)
	csd = c()
	for (i in 1:length(clustering)) {
		csd = c(csd, rep(palette()[i %% 8], times=length(clustering[[sort(names(clustering))[i]]])))
	}
	rsd = c()
	for (i in 1:length(signatures)) {
		rsd = c(rsd, rep(palette()[i %% 8], times=length(signatures[[sort(names(signatures))[i]]])))
	}
	heatmap.2(exprs[features, samples],
		trace="none",
		scale="none",
		col=colorpanel(49, low="blue", high="yellow"),
		breaks=seq(bounds[1], bounds[2], length.out=50),
		ColSideColors=csd,
		RowSideColors=rsd,
		Colv=NA,
		Rowv=NA,
		dendrogram="none")
}

##' Run iNMF subtyping of the specified data set.
##' @param exprs expression matrix with samples in columns and features in rows
##' @param signatures named list with subtyping signatures mapped to feature ids
##' used in the expression matrix. The following names have to be used "step1", 
##' "step2.c1" and "step2.c2". The function "signaturesList" can be used to 
##' obtain this list.
##' @param bootstrap boolean indicating whether running in bootstrap mode. If TRUE,
##' samples will be randomly selected
##' @param silhouette boolean indicating whether the silhouette should be calculated
##' @param plotHeatmaps boolean indicating whether clustering heatmaps should be 
##' saved in files
##' @param directory path to the directory where heatmaps will be saved
##' @param filePrefix file name prefix for the heatmap files
##' @return a named list with the clustering results
##' @author Andreas Schlicker
iNMF = function(exprs, signatures, bootstrap=FALSE, silhouette=TRUE, plotHeatmaps=TRUE, directory=".", filePrefix="") {
	if (bootstrap) {
		samples = sample(colnames(exprs), ncol(exprs), TRUE)
	} else {
		samples = colnames(exprs)
	}
	
	clust1 = subtype(exprs, signatures$step1, samples, silhouette)
	clust2.c1 = subtype(exprs, signatures$step2.c1, clust1$clustering$'1', silhouette)
	clust2.c2 = subtype(exprs, signatures$step2.c2, clust1$clustering$'2', silhouette)
	
	if (plotHeatmaps) {
		png(paste(directory, "/", filePrefix, "step1.png", sep=""), width=4000, height=3000, res=300)
		createHeatmap(exprs, clust1$clustering, signatures$step1)
		dev.off()
		
		png(paste(directory, "/", filePrefix, "step2_c1.png", sep=""), width=4000, height=3000, res=300)
		createHeatmap(exprs, clust2.c1$clustering, signatures$step2.c1)
		dev.off()
		
		png(paste(directory, "/", filePrefix, "step2_c2.png", sep=""), width=4000, height=3000, res=300)
		createHeatmap(exprs, clust2.c2$clustering, signatures$step2.c2)
		dev.off()
	}
	
	list(step1=clust1, step2.c1=clust2.c1, step2.c2=clust2.c2)
}

##' Creates a subtyping matrix. If a samples was assigned 
##' to a subtypes, the corresponding entry is set to 1 else to 0.
##' As the first level subtypes can be deduced from the second level
##' subtypes, they are not taken into account for creating the matrix.
##' @param solution subtyping solution as returned by iNMF
##' @param samples vector with all samples contained in the data set.
##' If set to NULL (default), the function assumes that the step1
##' subtyping contains all samples. This parameter allows to account for
##' samples missing in the solution, e.g. when performing bootstrapping. 
##' @return subtyping matrix with subtypes in columns and samples in rows
##' @author Andreas Schlicker
subtypingMatrix = function(solution, samples=NULL) {
	if (class(solution) != "list") {
		print(solution)
		print(class(solution))
		return(matrix(0, ncol=5, nrow=length(samples)))
	}
	combined = c(solution$step2.c1$clustering,
				 solution$step2.c2$clustering)
	if (is.null(samples)) {
		samples = unlist(combined)
	}
	mat = matrix(0, ncol=5, nrow=length(samples))
	colnames(mat) = c("1.1", "1.2", "1.3", "2.1", "2.2")
	rownames(mat) = samples
	
	for (n in names(combined)) {
		tab = table(combined[[n]])
		mat[names(tab), n] = mat[names(tab), n] + tab
	}
	
	mat
}

coClustering = function(solution, samples=NULL) {
	if (class(solution) != "list") {
		print(solution)
		print(class(solution))
		return(matrix(0, ncol=length(samples), nrow=length(samples)))
	}
	combined = c(solution$step2.c1$clustering,
				 solution$step2.c2$clustering)
	if (is.null(samples)) {
		samples = unlist(combined)
	}
	mat = matrix(0, ncol=length(samples), nrow=length(samples))
	colnames(mat) = samples
	rownames(mat) = samples
	
	for (n in names(combined)) {
		temp.samp = unique(combined[[n]])
		for (i in temp.samp) {
			mat[i, temp.samp] = mat[i, temp.samp] + 1
		}
	}
	
	mat
}

##' Runs iNMF subtyping in bootstrapping mode. For each run, the samples are sampled from
##' all samples in the data set. The result is a samples by subtype matrix giving the 
##' probability for each sample-subtype combination.
##' @param exprs expression matrix with samples in columns and features in rows
##' @param signatures named list with subtyping signatures mapped to feature ids
##' used in the expression matrix. The following names have to be used "step1", 
##' "step2.c1" and "step2.c2". The function "signaturesList" can be used to 
##' obtain this list.
##' @param runs number of bootstrapping runs; default: 1000
##' @param procCores number of processing cores to use; default: 1. The function uses the 
##' packages "foreach" and "doMC" for parallel computations. If these packages can not be
##' loaded, procCores is ignored and only one core is used.
##' @param seed random seed; default: NULL. If computation is performed in parallel, this
##' seed is not set for each thread. Otherwise all would produce the same random number sequence
##' and thus use the same sample selection.
##' @return a probability matrix with subtypes in columns and samples in rows. Each value
##' is computed as number of times the sample was assigned to a given subtype divided by
##' the total number of subtype assignments over all runs. 
##' @author Andreas Schlicker
bootstrappediNMF = function(exprs, signatures, runs=1000, procCores=1, seed=NULL) {
	if (!is.null(seed)) {
		set.seed(seed)
	}
	
	if (procCores > 1 && (!require(foreach) || !(require(doMC)))) {
		warning("Could not load package \"foreach\" or package \"doMC\". Using only one processing core.")
		procCores = 1
	}
	
	if (procCores > 1) {
		registerDoMC(cores=procCores)
		subs = mclapply(1:runs, 
				iNMF, 
				exprs=exprs, 
				signatures=signatures, 
				bootstrap=TRUE, 
				silhouette=FALSE, 
				plotHeatmaps=FALSE,
				mc.cores=procCores)	
	} else {
		subs = lapply(1:runs, 
				iNMF, 
				exprs=exprs, 
				signatures=signatures, 
				bootstrap=TRUE, 
				silhouette=FALSE, 
				plotHeatmaps=FALSE)
	}

	
	res = subtypingMatrix(subs[[1]], colnames(exprs))
	coclust = coClustering(subs[[1]], colnames(exprs))
	for (i in 2:runs) {
		res = res + subtypingMatrix(subs[[i]], colnames(exprs))
		coclust = coclust + coClustering(subs[[i]], colnames(exprs))
	}
	
	list(subtype.p=res / apply(res, 1, sum), cocluster=coclust / diag(coclust))
}
