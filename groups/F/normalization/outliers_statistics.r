# Count outliers per category for each dataset 
# 
# Author: Andreas Schlicker
######################################################################################


library(synapseClient)
library(rGithubClient)
library(stringr)

# GitHub repository
crcRepo = getRepo("andreas-schlicker/crcsc")
synapseHelper = getPermlink(crcRepo, "groups/F/pipeline/synapse_helper.r")
sourceRepoFile(crcRepo, "groups/F/pipeline/synapse_helper.r")
thisScript = getPermlink(crcRepo, "groups/F/normalization/outliers_statistics.r")

synapseLogin()

# Folder with all public datasets
pubFolder = "syn2176663"
privateFolder = "syn2157405"

outlierSummary = list()
outlierMats = list()

# Get result files
allData = rbind(synapseQuery(paste('SELECT id, name FROM entity WHERE parentId=="', pubFolder, '"', sep="")),
				synapseQuery(paste('SELECT id, name FROM entity WHERE parentId=="', privateFolder, '"', sep="")))


# Strict outlier calling.
# Any sample with outlier call in any category; ignoring MA plots 
strictOutliers = function(x) { 
	names(which(x[, "Sum"]-x[, "MAplot"] > 0))
}
		
# Relaxed outlier calling.
# Outlier needs to be call in Boxplot and either ArrayDist or SVD.
relaxedOutliers = function(x) {
	names(which((x[, "ArrayDist"] == 1 | x[, "SVD"] == 1) & x[, "Boxplot"] == 1))
}
		
for (i in 1:nrow(allData)) {
	# Get the files 
	files = synapseQuery(paste('SELECT id, name FROM entity WHERE parentId=="', allData[i, 2], '"', sep=""))
	outlierFile = which(str_detect(files[, 1], "outlier_summary"))
	
	if (length(outlierFile) > 0) {
		for (n in outlierFile) {
			# Get fRMA normalized gene expression
			failed = TRUE
			tries = 0
			while (failed && (tries < 5)) {
				outlierSum = tryCatch(loadMatrix(files[n, 2]),
						error=function(e) NA)
				if (class(outlierSum) == "matrix") {
					failed=FALSE
				}
				tries = tries + 1
			}
		
			outlierSummary[[files[n, 1]]] = apply(outlierSum, 2, function(x) { sum(x > 0) })
			outlierMats[[files[n, 1]]] = outlierSum
			
			# Call outliers
			outliers = list(strict=strictOutliers(outlierSum),
							relaxed=relaxedOutliers(outlierSum))
			
			baseFileName = gsub("summary.tsv", "", files[n, 1])
			for (k in names(outliers)) {
				# Write temporary outlier file
				filePath = file.path(tempdir(), paste(baseFileName, k, ".txt", sep=""))
				write.table(outliers[[k]], file=filePath, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
				
				# List with used resources
				resources = list(list(entity=files[outlierFile, 2], wasExecuted=F),
								 list(url=synapseHelper, name=basename(synapseHelper), wasExecuted=F),
								 list(url=thisScript, name=basename(thisScript), wasExecuted=T))
				
				# Store results in synapse and forget about the temporary file 
				synFile = File(path=filePath, parentId=allData[i, 2])
				failed = TRUE
				tries = 0
				while (failed && (tries < 5)) {
					res = tryCatch(synStore(synFile, used=resources),
							error=function(e) NA)
					if (!is.na(res)) {
						failed=FALSE
					}
					tries = tries + 1
				}
				unlink(filePath)
			}
		}
	}
	
	#summaryMat = matrix(0, ncol=4, nrow=length(outlierSummary))
	#colnames(summaryMat) = names(outlierSummary[[1]])[1:4]
	#rownames(summaryMat) = gsub("_outlier_summary.tsv", "", names(outlierSummary))
	#for (n in names(outlierSummary)) {
	#	summaryMat[gsub("_outlier_summary.tsv", "", n), ] = signif(outlierSummary[[n]][1:4] / nrow(outlierMats[[n]]), digits=2)
	#}
	
	#numberOutliers = cbind(Strict=sapply(names(outlierMats), function(x) { length(strictOutliers(outlierMats[[x]])) }),
	#		Relaxed=sapply(names(outlierMats), function(x) { length(relaxedOutliers(outlierMats[[x]])) }))
	#rownames(numberOutliers) = gsub("_outlier_summary.tsv", "", rownames(numberOutliers))
}

synapseLogout()
