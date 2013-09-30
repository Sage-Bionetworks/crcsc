# Apply the iNMF subtyping to all data sets in the CRC subtyping consortium.
# 
# Author: Andreas Schlicker
###############################################################################


library(synapseClient)
library(rGithubClient)

crcRepo = getRepo("andreas-schlicker/crcsc")

sourceRepoFile(crcRepo, "groups/F/pipeline/inmf.r")
iNMFFunctions = getPermlink(crcRepo, "groups/F/pipeline/inmf.r")
sourceRepoFile(crcRepo, "groups/F/pipeline/synapse_helper.r")
helperFunctions = getPermlink(crcRepo, "groups/F/pipeline/synapse_helper.r")
sourceRepoFile(crcRepo, "groups/F/pipeline/plotting.r")
plottingFunctions = getPermlink(crcRepo, "groups/F/pipeline/plotting.r")

thisScript = getPermlink(crcRepo, "groups/F/pipeline/subtyping.r")

# password will be request after calling this
synapseLogin("a.schlicker@nki.nl")

iNMFSignaturesSynId = "syn2245383"
iNMFSignatures = read.table(getFileLocation(synGet(iNMFSignaturesSynId)), header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)

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
		mdanderson=list(synId="syn2233387", sigId="symbol", mapSynId="syn2233216", mapId="GeneName", is.logr=TRUE))

publicExprList = list(
		gse10961=list(synId="syn2177194", is.logr=FALSE, annSynId="syn2177195"),
		gse13067=list(synId="syn2177888", is.logr=FALSE, annSynId="syn2177889"),
		gse13294=list(synId="syn2177894", is.logr=FALSE, annSynId="syn2177895"),
		gse14333=list(synId="syn2181079", is.logr=FALSE, annSynId="syn2181006"),
		gse15960=list(synId="syn2177199", is.logr=FALSE, annSynId="syn2177200"),
		gse17537=list(synId="syn2178128", is.logr=FALSE, annSynId="syn2178129"),
		gse20916=list(synId="syn2177899", is.logr=FALSE, annSynId="syn2177898"),
		gse2109=list(synId="syn2177169", is.logr=FALSE, annSynId="syn2177168"),
		gse23878=list(synId="syn2177902", is.logr=FALSE, annSynId="syn2178063"),
		gse37892=list(synId="syn2178082", is.logr=FALSE, annSynId="syn2178089"),
		gse4107=list(synId="syn2177179", is.logr=FALSE, annSynId="syn2177180"),
		gse4183=list(synId="syn2177187", is.logr=FALSE, annSynId="syn2177188"),
		gse8671=list(synId="syn2181088", is.logr=FALSE, annSynId="syn2181090"))


allResults = list()
for (n in c(names(coreExprList), names(publicExprList()))) {
	x = coreExprList[[n]]
	
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
		map.mat = loadMatrix(x$mapSynId)
		mapping = getMapping(map.mat, x$mapId)
	}
	# Get the iNMF signatures
	signatures = signaturesList(rownames(data.mat), sig.id, iNMFSignatures, mapping)
	
	# Run bootstrapped iNMF subtyping
	coreResults[[n]] = bootstrappediNMF(data.mat, signatures, runs=1000, procCores=15)
}

plotting.cols = c("gray10", "gray40", "gray70", "springgreen4", "springgreen2")
plotting.shapes = 15:19
for (x in names(allResults)) {
	png(paste(x, "probability.png", sep="_"), width=6000, height=3000, res=300)
	print(probabilityPlot(allResults[[x]][[1]], colors=plotting.cols))
	dev.off()
	
	png(paste(x, "probability_faceted.png", sep="_"), width=6000, height=3000, res=300)
	print(facetedProbabilityPlot(allResults[[x]][[1]], colors=plotting.cols))
	dev.off()
	
	png(paste(x, "margins.png", sep="_"), width=6000, height=3000, res=300)
	print(marginPlot(allResults[[x]][[1]], colors=plotting.cols, shapes=plotting.shapes))
	dev.off()
	
	png(paste(x, "cooccurrence.png", sep="_"), width=3000, height=3000, res=300)
	coclusteringPlot(allResults[[x]][[2]])
	dev.off()
}

