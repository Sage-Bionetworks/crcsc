# Apply the iNMF subtyping to all data sets in the CRC subtyping consortium.
# 
# Author: Andreas Schlicker
###############################################################################

library(synapseClient)
library(rGithubClient)
library(stringr)

# GitHib repository
crcRepo = getRepo("andreas-schlicker/crcsc")
# iNMF subtyping functions
sourceRepoFile(crcRepo, "groups/F/pipeline/inmf.r")
iNMFFunctions = getPermlink(crcRepo, "groups/F/pipeline/inmf.r")
# Functions for interacting with Synapse 
sourceRepoFile(crcRepo, "groups/F/pipeline/synapse_helper.r")
helperFunctions = getPermlink(crcRepo, "groups/F/pipeline/synapse_helper.r")
# Putting everything together
thisScript = getPermlink(crcRepo, "groups/F/pipeline/subtyping_celllines.r")

synapseLogin()

# File with iNMF signatures (Affy probe sets, Entrez gene IDs and gene symbols) stored in Synapse
iNMFSignaturesSynId = "syn2245383"
iNMFSignatures = read.table(getFileLocation(synGet(iNMFSignaturesSynId)), header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)

# Parameters needed to subtype each data set
# synID: Synapse ID for the expression data
# sigID: which type of iNMF signature do we need? Possible values: ps, symbol, entrez; default: ps
# mapSynId: File mapping data set-sepcific ids to iNMF signature IDs
# mapId: name of column in the mapSynId file
# file: file to extract from zipped source files
# is.logr: is the data already log-ratio?
# sep: separator used in the data file; default: \t
# quote: quote symbol used in the data file; default: none
# annSynId: Synapse ID of the file containing sample annotation
allData = list(gse8332=list(synId="syn2181082", is.logr=TRUE, loc="last"),
			   gsk=list(synId="syn2181084", is.logr=FALSE, loc="first"),
			   sanger=list(synId="syn2181097", is.logr=FALSE),
			   ccle=list(synId="syn2292137", is.logr=FALSE, sigId="entrez"))

# Subtyping
allResults = list()
for (n in names(allData)) {
	x = allData[[n]]
	
	sep = "\t"
	if (!is.null(x$sep)) {
		sep = x$sep
	}
	
	quote = ""
	if (!is.null(x$quote)) {
		quote = x$quote
	}
	
	# Get the expression data
	file = ""
	if (!is.null(x$file)) {
		file = x$file
	}
	data.mat = loadMatrix(x$synId, file=file, sep=sep, quote=quote)
	
	if (n == "ccle") {
		rownames(data.mat) = unlist(lapply(str_split(rownames(data.mat), "_"), function(x) { x[1] }))
	}
	
	if (!x$is.logr) {
		data.mat = data.mat - apply(data.mat, 1, mean)
	}
	
	if (!is.null(x$loc)) {
		data.mat = averageReplicates(data.mat, x$loc)
	}
	
	sig.id = "ps"
	if (!is.null(x$sigId)) {
		sig.id = x$sigId
	}
	
	mapping = NULL
	if (!is.null(x$mapSynId)) {
		# Get the mapping from data set IDs to signature IDs
		map.mat = loadMatrix(x$mapSynId)
		mapping = getMapping(map.mat, x$mapId)
	}
	# Get the iNMF signatures
	signatures = signaturesList(rownames(data.mat), sig.id, iNMFSignatures, mapping)
	
	# Run bootstrapped iNMF subtyping
	allResults[[n]] = bootstrappediNMF(data.mat, signatures, runs=10000, procCores=20)
}

# Save the results in Synapse
synResultDir = "syn2274068"
groupName = "GroupF"

lapply(names(allResults), 
	   function(datasetName){
			# Round probability values
		    iNMFRes = round(allResults[[datasetName]][[1]], digits=4)
			colnames(iNMFRes) = paste("subtype", colnames(iNMFRes), sep="")
			# And create data frame
			iNMF.df = data.frame(sampleName=rownames(iNMFRes), iNMFRes)
			# Synapse ID of the expression data
			synId = allData[[datasetName]]$synId
			# That's where we got the ID mapping got from. If any was used
			mapSynId = allData[[datasetName]]$mapSynId
			
			filePath = file.path(tempdir(), paste(groupName,"_",synId,"_",datasetName,".tsv",sep=""))
			write.table(iNMF.df, file=filePath, sep="\t", quote=FALSE, row.names=FALSE)
			
			# List with used resources
			resources = list(list(entity=iNMFSignaturesSynId, wasExecuted=F),
							 list(entity=synId, wasExecuted=F),
							 list(url=iNMFFunctions, name=basename(iNMFFunctions), wasExecuted=F),
					 		 list(url=helperFunctions, name=basename(helperFunctions), wasExecuted=F),
							 list(url=thisScript, name=basename(thisScript), wasExecuted=T))
			if (!is.null(mapSynId)) {
				resources = c(resources, list(list(entity=mapSynId, wasExecuted=F)))
			}
			
			# Store results in synapse and forget about the temporary file 
			synFile = File(path=filePath, parentId=synResultDir)
			synFile = synStore(synFile, used=resources)
			unlink(filePath)
		})

synapseLogout()
