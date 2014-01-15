# Apply the iNMF subtyping to all data sets in the CRC subtyping consortium.
# 
# Author: Andreas Schlicker
###############################################################################

library(synapseClient)
library(rGithubClient)

# GitHib repository
crcRepo = getRepo("andreas-schlicker/crcsc")
# iNMF subtyping functions
sourceRepoFile(crcRepo, "groups/F/pipeline/inmf.r")
iNMFFunctions = getPermlink(crcRepo, "groups/F/pipeline/inmf.r")
# Functions for interacting with Synapse 
sourceRepoFile(crcRepo, "groups/F/pipeline/synapse_helper.r")
helperFunctions = getPermlink(crcRepo, "groups/F/pipeline/synapse_helper.r")
# Putting everything together
thisScript = getPermlink(crcRepo, "groups/F/pipeline/subtyping.r")

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
coreExprList = list(
		agendia_gse42284=list(synId="syn2192792", sigId="symbol", mapSynId="syn2192791", mapId="symbol", file="GSE42284_normalized_data_matrix.txt", is.logr=TRUE),
		agendia_ico208=list(synId="syn2192796", sigId="symbol", mapSynId="syn2192791", mapId="symbol", file="ICO208_normalized_data.txt", is.logr=TRUE),
		agendia_vhb70=list(synId="syn2192799", sigId="symbol", mapSynId="syn2192801", mapId="GeneName", file="VHB70_normalized_data.txt", is.logr=TRUE),
		kfsyscc=list(synId="syn2169565", is.logr=FALSE),
		french=list(synId="syn2171434", is.logr=TRUE),
		amc_ajccii=list(synId="syn2159423", sep=",", quote="\"", is.logr=FALSE),
		nki_az=list(synId="syn2176657", is.logr=FALSE),
		petacc3=list(synId="syn2175581", sigId="entrez", mapSynId="syn2199825", mapId="Entrez.GeneID", is.logr=FALSE),
		tcga_rnaseq=list(synId="syn2161141", sigId="symbol", is.logr=FALSE),
		mdanderson=list(synId="syn2233387", sigId="symbol", mapSynId="syn2233216", mapId="GeneName", is.logr=TRUE),
		tcga_rnaseq_ga=list(synId="syn2326094", sigId="symbol", is.logr=FALSE),
		tcga_rnaseq_hi=list(synId="syn2326100", sigId="symbol", is.logr=FALSE),
		tcga_rnaseq_merged=list(synId="syn2325328", sigId="symbol", is.logr=FALSE),
		tcga_microarray=list(synId="syn2316354", quote="\"", sigId="symbol", mapSynId="syn2316355", mapId="Gene.Symbol", is.logr=FALSE))

publicExprList = list(
		gse10961=list(synId="syn2177194", is.logr=FALSE, annSynId="syn2177195"),
		gse13067=list(synId="syn2177888", is.logr=FALSE, annSynId="syn2177889"),
		gse13294=list(synId="syn2177894", is.logr=FALSE, annSynId="syn2177895"),
		gse14333=list(synId="syn2181079", is.logr=FALSE, annSynId="syn2181006"),
		gse15960=list(synId="syn2177199", is.logr=FALSE, annSynId="syn2177200"),
		gse17536=list(synId="syn2178137", is.logr=FALSE, annSynId="syn2178136"),
		gse17537=list(synId="syn2178128", is.logr=FALSE, annSynId="syn2178129"),
		gse20916=list(synId="syn2177899", is.logr=FALSE, annSynId="syn2177898"),
		gse2109=list(synId="syn2177169", is.logr=FALSE, annSynId="syn2177168"),
		gse23878=list(synId="syn2177902", is.logr=FALSE, annSynId="syn2178063"),
		gse37892=list(synId="syn2178082", is.logr=FALSE, annSynId="syn2178089"),
		gse4107=list(synId="syn2177179", is.logr=FALSE, annSynId="syn2177180"),
		gse4183=list(synId="syn2177187", is.logr=FALSE, annSynId="syn2177188"),
		gse8671=list(synId="syn2181088", is.logr=FALSE, annSynId="syn2181090"))

allData = c(coreExprList, publicExprList)

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
	
	if (!x$is.logr) {
		data.mat = data.mat - apply(data.mat, 1, mean)
	}
	
	sig.id = "ps"
	if (!is.null(x$sigId)) {
		sig.id = x$sigId
	}
	
	mapping = NULL
	if (!is.null(x$mapSynId)) {
		# Get the mapping from data set IDs to signature IDs
		map.mat = loadMatrix(x$mapSynId, quote=quote)
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
