# Functions to help with data handling using Synapse
# 
# Author: Andreas Schlicker
###############################################################################

##' Load a matrix from Synapse. Can be used for 
##' data matrices and for annotation matrices.
##' @param synId the Synapse file ID
##' @param file a filename to be used for zip files; default = ""
##' @param sep value separator used in the file; default = "\t"
##' @return a matrix with all the read values
##' @author Andreas Schlicker; based on a function from Justin Guinney (justin.guinney@sagebase.org)
loadMatrix = function(synId, file="", sep="\t", quote=""){
	# Get the file name
	filename = getFileLocation(synGet(synId))
	
	# Unzip zip files
	if (length(grep("zip$", filename)) > 0) {
		temp = unz(filename, file)
	} else {
		temp = filename
	}
	
	# Read the data
	raw = read.table(temp, 
					 sep=sep,
					 header=TRUE,
			 		 quote=quote,
					 stringsAsFactors=FALSE,
					 as.is=TRUE,
					 fill=TRUE,
					 comment.char="",
					 check.names=FALSE)
	
	# Remove duplicates
	dupls = which(duplicated(raw[, 1]))
	if (length(dupls) > 0) {
		raw = raw[-dupls, ]
	}
	
	# Create the matrix
	# Some annotation matrices only have two columns. Keep both of them
	if (ncol(raw) > 2 && is.character(raw[, 1])) {
		mat = as.matrix(raw[, -1])
	} else {
		mat = as.matrix(raw)
	}
	if (is.character(raw[, 1])) {
		rownames(mat) = raw[, 1]
	}
	
	mat
}

##' Create mapping of data set-specific IDs to signature IDs.
##' @param anno an annotation matrix; rownames have to be equal to the IDs
##' used in the expression matrix
##' @param colName name of the column that contains the ID from the signature
##' @return a mapping matrix
##' @author Andreas Schlicker
getMapping = function(anno, colName) {
	map = cbind(ids=rownames(anno), anno[, colName])
	rownames(map) = rownames(anno)
	
	map
}

##' Average replicates in the matrix. Colnames are split at "_" and 
##' the element at "loc" in the resulting vector is used to identify replicates.
##' The name of the first replicate is used as column name in the output.
##' @param mat data matrix, samples in columns and features in rows
##' @param loc "first" or "last", element of the column name (separated by "_")
##' that identifies replicates; default: "first"
##' @return data matrix with averaged replicates
##' @author Andreas Schlicker
averageReplicates = function(mat, loc=c("first", "last")) {
	require(stringr) || stop("Can't load required package \"stringr\"!")
	# Get the cell line name
	loc = match.arg(loc)
	clName = unlist(lapply(str_split(colnames(mat), "_"), function(x) { if (loc == "last") el=length(x) else el = 1; x[el] }))
	# Which samples map to which cell line
	replicates = lapply(unique(clName), function(x) { which(clName == x) })
	names(replicates) = unique(clName)
	
	resMat = matrix(NA, nrow=nrow(mat), ncol=length(replicates))
	rownames(resMat) = rownames(mat)
	colnames(resMat) = unlist(lapply(replicates, function(x) { colnames(mat)[x[1]] }))
	for (i in 1:length(replicates)) {
		resMat[, i] = apply(mat[, replicates[[i]], drop=FALSE], 1, mean)
	}
	
	resMat
}
