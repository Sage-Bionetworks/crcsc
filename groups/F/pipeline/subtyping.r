# Apply the iNMF subtyping to all data sets in the CRC subtyping consortium.
# 
# Author: Andreas Schlicker
###############################################################################


library(synapseClient)
library(rGithubClient)

crcRepo = getRepo("andreas-schlicker/crcsc")
iNMFSignatures = getPermlink(crcRepo, "groups/F/pipeline/inmf_signatures.rdata")
load(iNMFSignatures)
sourceRepoFile(crcRepo, "groups/F/pipeline/inmf.r")
iNMFFunctions = getPermlink(crcRepo, "groups/F/pipeline/inmf.r")
sourceRepoFile(crcRepo, "groups/F/pipeline/synapse_helper.r")
helperFunctions = getPermlink(crcRepo, "groups/F/pipeline/synapse_helper.r")
thisScript = getPermlink(crcRepo, "groups/F/pipeline/subtyping.r")

# password will be request after calling this
synapseLogin("a.schlicker@nki.nl")

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

coreResults = lapply(coreExprList, function(x) {
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
			if (!is.null(x$file))
				file = x$file
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
			signatures = signaturesList(rownames(data.mat), sig.id, mapping)
			
			# Run bootstrapped iNMF subtyping
			bootstrappediNMF(data.mat, signatures, runs=10, procCores=10)
		})


publicExprList = list(
		gse10961=toEntrez("syn2177194","syn2177195"),
		gse13067=toEntrez("syn2177888","syn2177889"),
		gse13294=toEntrez("syn2177894","syn2177895"),
		gse14333=toEntrez("syn2181079","syn2181006"),
		gse15960=toEntrez("syn2177199","syn2177200"),
		gse17537=toEntrez("syn2178128","syn2178129"),
		gse20916=toEntrez("syn2177899","syn2177898"),
		gse2109=toEntrez("syn2177169","syn2177168"),
		gse23878=toEntrez("syn2177902","syn2178063"),
		gse37892=toEntrez("syn2178082","syn2178089"),
		gse4107=toEntrez("syn2177179","syn2177180"),
		gse4183=toEntrez("syn2177187","syn2177188"),
		gse8671=toEntrez("syn2181088","syn2181090"))

allData = c(coreExprList,publicExprList)

