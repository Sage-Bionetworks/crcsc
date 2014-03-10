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
library(arrayQualityMetrics)

# GitHub repository
crcRepo = getRepo("andreas-schlicker/crcsc")
synapseHelper = getPermlink(crcRepo, "groups/F/pipeline/synapse_helper.r")
sourceRepoFile(crcRepo, "groups/F/pipeline/synapse_helper.r")
thisScript = getPermlink(crcRepo, "groups/F/normalization/outliers_publicdata.r")

synapseLogin()

# Folder with all public datasets
pubFolder = "syn2176663"

# Get result files
allData = synapseQuery(paste('SELECT id, name FROM entity WHERE parentId=="', pubFolder, '"', sep=""))
# Remove possible confidence files
allData = allData[str_detect(allData[, "entity.name"], "GSE"), ]

for (i in 1:nrow(allData)) {
	# Get the files 
	files = synapseQuery(paste('SELECT id, name FROM entity WHERE parentId=="', allData[i, 2], '"', sep=""))
	frmaFile = which(str_detect(files[, 1], "frma"))
	geoId = unlist(str_split(allData[i, 1], "_"))[1]
	outlierDir = file.path("OUTLIERS", geoId)
	dir.create(outlierDir, recursive=TRUE)
	
	if (length(frmaFile) > 0) {
		# Get fRMA normalized gene expression
		exprsMat = loadMatrix(files[frmaFile, 2])
		# Dummy sample annotation
		annMat = data.frame(ID=colnames(exprsMat))
		rownames(annMat) = colnames(exprsMat)
		
		exprsSet = new("ExpressionSet", exprs=exprsMat, 
					   phenoData=new("AnnotatedDataFrame", data=annMat, 
							   varMetadata=data.frame(labelDescription=c("Sample"), row.names=c("ID"))), 
				annotation="hgu133plus2")
		
		# Perform SVD outlier detection 
		isOutlier = outlierDetection(exprsMat, 
				fileName=file.path(outlierDir, "outliers_svd.png"))
		
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
		filePath = file.path(tempdir(), paste(geoId, "_outlier_summary.tsv", sep=""))
		write.table(outlierTable, file=filePath, sep="\t", quote=FALSE)
		
		# List with used resources
		resources = list(list(entity=files[rawFile, 2], wasExecuted=F),
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

synapseLogout()
