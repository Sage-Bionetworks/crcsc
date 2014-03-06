# Normalize all public datasets using fRMA and store data in Synapse 
# 
# Author: Andreas Schlicker
###############################################################################

library(synapseClient)
library(rGithubClient)
library(frma)
library(stringr)
library(affy)

# GitHib repository
crcRepo = getRepo("andreas-schlicker/crcsc")
thisScript = getPermlink(crcRepo, "groups/F/normalization/frma_datasets.r")

synapseLogin()

# Get result files
allData = list(amc=list(synId="syn2019118", cel="syn2027094", prefix="GSE33113", filename="GSE33113_RAW.tar"),
			   french=list(synId="syn2019116", cel="syn2026929", prefix="GSE39582", filename="GSE39582_RAW.tar"),
			   kfsyscc=list(synId="syn2019114", cel="syn2025141", prefix="KFSYSCC", filename="kfsyscc.colon.tar.gz"),
			   nki=list(synId="syn2176651", cel="syn2176652", prefix="GSE35896", filename="GSE35896_RAW.tar"))

for (i in names(allData)) {
	# Get the archive
	#filename = getFileLocation(synGet(allData[[i]]$cel))
	tD = file.path("TMP", allData[[i]]$filename)
	dir.create(tD)
	
	command = paste("tar xf ", allData[[i]]$filename, " -C ", tD, sep="")
	if (str_detect(allData[[i]]$filename, "gz")) {
		command = paste("gunzip -c ", allData[[i]]$filename, " | tar xf - -C ", tD, sep="")
	}
	system(command)
	
	es = frma(ReadAffy(celfile.path=tD, compress=TRUE), summarize="robust_weighted_average")
	colnames(es) = gsub(".gz", "", colnames(es))
	unlink(tD, recursive=TRUE)
		
	# Write temporary file with expression data
	filePath = file.path(tempdir(), paste(allData[[i]]$prefix, "_frma_expression.tsv", sep=""))
	write.table(exprs(es), file=filePath, sep="\t", quote=FALSE)
		
	# List with used resources
	resources = list(list(entity=files[rawFile, 2], wasExecuted=F),
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
