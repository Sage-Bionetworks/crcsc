# Assign confidence call to all subtype assignments and upload that to Synapse
# 
# Author: schlandi
###############################################################################

library(synapseClient)
library(rGithubClient)

# GitHib repository
crcRepo = getRepo("andreas-schlicker/crcsc")
# Plotting functions for calculating the margin
sourceRepoFile(crcRepo, "groups/F/pipeline/plotting.r")
plottingFunctions = getPermlink(crcRepo, "groups/F/pipeline/plotting.r")
# This script
thisScript = getPermlink(crcRepo, "groups/F/pipeline/confidence_call.r")

synapseLogin()

# Load results from phase 1 subtyping
load("subtyping_inmf_results_phase1.rdata")

allData = list(agendia_gse42284="syn2192792", agendia_ico208="syn2192796", agendia_vhb70="syn2192799", 
			   kfsyscc="syn2169565", french="syn2171434", amc_ajccii="syn2159423", nki_az="syn2176657",
			   petacc3="syn2175581", tcga_rnaseq="syn2161141", mdanderson="syn2233387", tcga_rnaseq_ga="syn2326094",
			   tcga_rnaseq_hi="syn2326100", tcga_rnaseq_merged="syn2325328", tcga_microarray="syn2316354",
			   gse10961="syn2177194", gse13067="syn2177888", gse13294="syn2177894", gse14333="syn2181079",
			   gse15960="syn2177199", gse17536="syn2178137", gse17537="syn2178128", gse20916="syn2177899",
			   gse2109="syn2177169", gse23878="syn2177902", gse37892="syn2178082", gse4107="syn2177179",
			   gse4183="syn2177187", gse8671="syn2181088")

# Save the results in Synapse
synResultDir = "syn2274068"
groupName = "GroupF"

lapply(names(allResults), 
		function(datasetName){
			# Round probability values
			iNMFRes = apply(allResults[[datasetName]][[1]], 1, margin)
			# And create data frame
			iNMF.df = data.frame(sampleName=names(iNMFRes), confidence=ifelse(iNMFRes > 0.1, "HIGH", "LOW"))
			
			# Synapse ID of the expression data
			synId = allData[[datasetName]]
			
			filePath = file.path(tempdir(), paste(groupName, "_", synId, "_", datasetName, "_conf.tsv",sep=""))
			write.table(iNMF.df, file=filePath, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
			
			# List with used resources
			resources = list(list(url=plottingFunctions, name=basename(plottingFunctions), wasExecuted=F),
							 list(url=thisScript, name=basename(thisScript), wasExecuted=T))
			
			# Store results in synapse and forget about the temporary file 
			synFile = File(path=filePath, parentId=synResultDir)
			synFile = synStore(synFile, used=resources)
			unlink(filePath)
		})

synapseLogout()
