# Assign confidence call to all subtype assignments and upload that to Synapse
# 
# Author: schlandi
###############################################################################

library(synapseClient)
library(rGithubClient)
library(stringr)

# GitHib repository
crcRepo = getRepo("andreas-schlicker/crcsc")
# Plotting functions for calculating the margin
sourceRepoFile(crcRepo, "groups/F/pipeline/plotting.r")
plottingFunctions = getPermlink(crcRepo, "groups/F/pipeline/plotting.r")
# This script
thisScript = getPermlink(crcRepo, "groups/F/pipeline/confidence_call.r")

synapseLogin()

# Save the results in Synapse
synResultDir = "syn2274068"
groupName = "GroupF"

# Get result files
allData = synapseQuery(paste('SELECT id, name FROM entity WHERE parentId=="', synResultDir, '"', sep=""))
# Remove possible confidence files
allData = allData[!str_detect(allData[, "entity.name"], "_conf"), ]

apply(allData, 1,  
	  function(dataset){
		  	# Get the subtyping from Synapse
		  	tempRes = synGet(dataset[, 2])
		  	# and calculate the margin
		  	iNMFRes = apply(read.delim(tempRes@filePath, header=TRUE, row.names=1, sep="\t", quote="", as.is=TRUE), 1, margin)
			# create data frame
			iNMF.df = data.frame(sampleName=names(iNMFRes), confidence=ifelse(iNMFRes > 0.1, "HIGH", "LOW"))
			
			filePath = file.path(tempdir(), paste(str_replace(dataset[, 1], ".tsv", ""), "_conf.tsv", sep=""))
			write.table(iNMF.df, file=filePath, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
			
			# List with used resources
			resources = list(list(entity=dataset[, 2], wasExecuted=F),
							 list(url=plottingFunctions, name=basename(plottingFunctions), wasExecuted=F),
							 list(url=thisScript, name=basename(thisScript), wasExecuted=T))
			
			# Store results in synapse and forget about the temporary file 
			synFile = File(path=filePath, parentId=synResultDir)
			synFile = synStore(synFile, used=resources)
			unlink(filePath)
		})

synapseLogout()
