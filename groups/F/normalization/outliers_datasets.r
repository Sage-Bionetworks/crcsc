# Normalize all public datasets using fRMA and store data in Synapse 
# 
# Author: Andreas Schlicker
# Method outlierDetection was provided by Justin Guinney (justin.guinney@sagebase.org)
######################################################################################

# Flag samples as outliers based on SVD.
outlierDetection = function(x, numPC=2, nsd=2.5, 
		makePlot=TRUE, fileName=NULL){
	pc = svd(x - rowMeans(x))
	
	v = scale(pc$v[,1:numPC])
	d = apply(v,1,function(x){ sqrt(sum(x^2))})
	isOutlier = as.integer(d > nsd)
	names(isOutlier) = colnames(x)
	
	if(makePlot){
		require(plotrix)
		png(fileName, width=4000, height=3000, res=300)
		par(mfrow=c(1,2))
		lim = max(c(d,nsd))
		plot((pc$d^2 / sum(pc$d^2))[1:20],ylab="% var explained",xlab="PCs")
		plot(v[,1],v[,2],xlim=c(-lim,lim),ylim=c(-lim,lim),asp=1,
				pch=19,cex=.7,col=c("black","red")[factor(isOutlier)],
				xlab="PC1",ylab="PC2")
		draw.circle( 0, 0, nsd,border="blue",lty=2)
		dev.off()
	}
	
	return(isOutlier)
}

library(synapseClient)
library(rGithubClient)
library(stringr)
library(Biobase)
library(arrayQualityMetrics)

# GitHub repository
crcRepo = getRepo("andreas-schlicker/crcsc")
synapseHelper = getPermlink(crcRepo, "groups/F/pipeline/synapse_helper.r")
sourceRepoFile(crcRepo, "groups/F/pipeline/synapse_helper.r")
thisScript = getPermlink(crcRepo, "groups/F/normalization/outliers_datasets.r")

synapseLogin()

allData = list(amc=list(synId="syn2019118", exprs="syn2363559", prefix="GSE33113"),
	 		   french=list(synId="syn2019116", exprs="syn2363561", prefix="GSE39582"),
			   kfsyscc=list(synId="syn2019114", exprs="syn2363564", prefix="KFSYSCC"),
			   nki=list(synId="syn2176651", exprs="syn2363566", prefix="GSE35896"),
			   gse42284=list(synId="syn2019117", exprs="syn2192792", prefix="GSE42284", file="GSE42284_normalized_data_matrix.txt"),
			   ico208=list(synId="syn2019117", exprs="syn2192796", prefix="ICO208", file="ICO208_normalized_data.txt"),
			   vhb70=list(synId="syn2019117", exprs="syn2192799", prefix="VHB70", file="VHB70_normalized_data.txt"),
			   mdanderson=list(synId="syn2026954", exprs="syn2233387", prefix="MDA219"),
			   petacc3=list(synId="syn2172604", exprs="syn2175581", prefix="PETACC3"),
			   tcga=list(synId="syn2023932", exprs="syn2325328", prefix="TCGACRC_expression_merged"))

for (i in 1:length(allData)) {
	# Create the directory for the result plots
	outlierDir = file.path("OUTLIERS", allData[[i]]$prefix)
	dir.create(outlierDir, recursive=TRUE)
	
	# Get normalized gene expression
	exprsMat = loadMatrix(allData[[i]]$exprs, file=ifelse(!is.null(allData[[i]]$file), allData[[i]]$file, ""))
	exprsSet = ExpressionSet(assayData=exprsMat)
		
	# Perform SVD outlier detection 
	isOutlier = outlierDetection(exprsMat, fileName=file.path(outlierDir, "outliers_svd.png"))
		
	# Run arrayQualityMetrics
	aqm = arrayQualityMetrics(exprsSet, outdir=file.path(outlierDir, "AQM"), spatial=FALSE)
		
	outlierTable = aqm$arrayTable[, c(3:5)]
	rownames(outlierTable) = aqm$arrayTable[, 2]
	colnames(outlierTable)[1:3] = c("ArrayDist", "Boxplot", "MAplot")
	outlierTable[outlierTable == ""] = 0
	outlierTable[outlierTable == "x"] = 1
	outlierTable = cbind(outlierTable, SVD=isOutlier[rownames(outlierTable)])
	outlierTable = cbind(outlierTable, Sum=rowSums(data.matrix(outlierTable)))
		
	# Write temporary file with expression data
	filePath = file.path(tempdir(), paste(allData[[i]]$prefix, "_outlier_summary.tsv", sep=""))
	write.table(outlierTable, file=filePath, sep="\t", quote=FALSE)
		
	# List with used resources
	resources = list(list(entity=allData[[i]]$exprs, wasExecuted=F),
					 list(url=synapseHelper, name=basename(synapseHelper), wasExecuted=F),
					 list(url=thisScript, name=basename(thisScript), wasExecuted=T))
				
	# Store results in synapse and forget about the temporary file 
	synFile = File(path=filePath, parentId=allData[[i]]$synId)
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

synapseLogout()
