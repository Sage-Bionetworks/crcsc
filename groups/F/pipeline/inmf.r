# Functions for performing the iNMF subtyping
# 
# Author: Andreas Schlicker (a.schlicker@nki.nl)
###############################################################################

##' Extract the given signature from the data frame containing all iNMF signatures.
##' @param type ID type to be used for the signature, either "ps", "symbol" or "entrez"
##' @param subtype either "1", "1.1", "1.2", "1.3", "2", "2.1", "2.2"
##' @param signatures data frame with all signatures
##' @return vector with all IDs of given type belonging to the requested signature
##' @author Andreas Schlicker
getSig = function(type, subtype, signatures) {
	signatures[which(signatures[, "type"]==type & signatures[, "subtype"]==subtype), "acc"]
}

##' Map iNMF signatures to identifiers used in another data set.
##' @param step one of "step1", "step2.c1", "step2.c2"
##' @param type one of "ps", "symbol" or "entrez". "ps" uses Affy HGU133+2 signatures as 
##' starting point. "symbol" uses signatures that were already mapped to 
##' gene symbols. "entrez" uses signatures that were previously mapped to Entrez gene IDs.
##' @param signatures data frame with all signatures
##' @param mapping a matrix mapping data set-specific IDs (first column) to either Affy 
##' HGU133+2 probe sets or gene symbols (second column). If NULL (default), the 
##' signature IDs will be kept.
##' @return a named list with the signatures
##' @author Andreas Schlicker
signatureList = function(step=c("step1", "step2.c1", "step2.c2"), type=c("ps", "symbol", "entrez"), signatures, mapping=NULL) {
	step = match.arg(step)
	type = match.arg(type)
	
	res = list()
	if (step == "step1") {
		res[["1"]] = getSig(type, "1", signatures)
		res[["2"]] = getSig(type, "2", signatures)
	} else if (step == "step2.c1") {
		res[["1.1"]] = getSig(type, "1.1", signatures)
		res[["1.2"]] = getSig(type, "1.2", signatures)
		res[["1.3"]] = getSig(type, "1.3", signatures)
	} else if (step == "step2.c2") {
		res[["2.1"]] = getSig(type, "2.1", signatures)
		res[["2.2"]] = getSig(type, "2.2", signatures)
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
##' @param type one of "ps", "symbol" or "entrez". "ps" uses Affy HGU133+2 signatures as 
##' starting point. "symbol" uses signatures that were already mapped to 
##' gene symbols. "entrez" uses signatures that were previously mapped to Entrez gene IDs.
##' @param signatures data frame with all signatures
##' @param mapping a matrix mapping data set-specific IDs (first column) to either Affy 
##' HGU133+2 probe sets or gene symbols (second column). If NULL (default), the 
##' signature IDs will be kept.
##' @return named list with all three signatures
##' @author Andreas Schlicker
signaturesList = function(allIds, type=c("ps", "symbol", "entrez"), signatures, mapping=NULL) {
	type = match.arg(type)
	
	list(step1= filterSignatures(signatureList("step1", type, signatures, mapping), allIds),
		 step2.c1 = filterSignatures(signatureList("step2.c1", type, signatures, mapping), allIds),
		 step2.c2 = filterSignatures(signatureList("step2.c2", type, signatures, mapping), allIds))
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
	# Calculate the margin for each cluster assignment
	margin = apply(means, 1, function(x) { sort(x, decreasing=TRUE)[1] - sort(x, decreasing=TRUE)[2] })
	
	# Mapping between subtype label and cluster ID
	mapping = rep("", times=length(names(signatures)))
	names(mapping) = names(signatures)
	# All cluster IDs
	clusts = colnames(means)
	# Cycle through all subtypes in order of decreasing margin
	for (n in names(sort(margin, decreasing=TRUE))) {
		if (length(clusts) > 1) {
			# There are more than one cluster left
			# Take one of the subtypes with the largest means
			mapping[n] = names(which(means[n, clusts] == max(means[n, clusts])))[1]
			clusts = setdiff(clusts, mapping[n])
		} else {
			# Assign the remaining subtype
			mapping[n] = clusts
		}
	}	
	# Reverse the mapping between subtype and cluster ID
	res = names(mapping)
	names(res) = mapping
	
	res
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
	
	if (length(samples) == 1) {
		means = unlist(lapply(signatures, function(x) { mean(exprs[x, samples]) }))
		res = list()
		res[names(means)[which(means == max(means))][1]] = c(samples)
		return(list(clustering=res))
	}
	
	dist.mat = 1 - cor(exprs[unlist(signatures), samples], use="pairwise.complete.obs")
	
	clustering = cutree(hclust(as.dist(dist.mat), method="complete"), k=min(length(signatures), length(samples)))
	
	clust.list = list()
	for (i in unique(clustering)) {
		clust.list[[as.character(i)]] = names(clustering[clustering == i])	
	}
	
	idMap = assignClusterId(clust.list, signatures, exprs)

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
