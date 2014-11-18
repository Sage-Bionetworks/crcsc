# TODO: Add comment
# 
# Author: pangelin
###############################################################################



if (!exists("work.dir")) {
	work.dir <- "/export/scratch/paolo/SageCRCsubtyping/analysis/work/"
}

# Set work dir
setwd(work.dir)
out.dir <- work.dir


# Load list of filtered genes 
load("gene.filtered.0.75.0.5.Rdata")

dataset.list <- c("agendia_gse42284", "agendia_ico208",  "agendia_vhb70",    "amc_ajccii",      
	"french",           "kfsyscc",          "mdanderson",       "nki_az",          
    "petacc3",          "tcgacrc_merged",   "gse10961",		"gse13067",         "gse13294",         
	"gse14333",        
	"gse15960",         "gse17537",         "gse20916",         "gse2109" ,        
    "gse23878",         "gse37892",         "gse4107",          "gse4183" ,        
	"gse8671",          "gse17536"        
)

# Load dataset list
if (!exists("dataset.list")) {
	load("dataset.list.Rdata")
}


# Load datasets
if (!exists("M.sets")) {
	message(paste("Loading datasets for merging..."))	
	M.sets <- list()
	for (n in dataset.list) {
		message(paste("Now processing ", n, "..."))	
		load(paste0(n,".selected.norm.Rdata"))
		M.sets[[n]] <- Mgenes		
	}
}

m.ALL <- NULL
rename.sample <- FALSE
for (n in names(M.sets)) {
	M.sets[[n]] <- M.sets[[n]][,gene.list]
	if (rename.sample) {
		rownames(M.sets[[n]]) <- paste0(n,'-',rownames(M.sets[[n]]))
	}
	m.ALL <- rbind(m.ALL, M.sets[[n]])
}

datasets <- list()
for (n in names(M.sets)) {
	datasets[[n]] <- rep(n,nrow(M.sets[[n]]))
}

datasets <- unlist(datasets)
names(datasets) <- NULL
if (! rename.sample){
	attr(m.ALL, "dataset") <- unlist(datasets)
	save(file="merged.datasets.Rdata", m.ALL)	
} else {
	save(file="merged.datasets.renamed.Rdata", m.ALL)	
}





library(synapseClient)
library(rGithubClient)

crcRepo <- getRepo("sib-bcf/crcsc")

sourceRepoFile(crcRepo, "groups/A/pipeline/synapseHelper.R")
helperFunctions <- getPermlink(crcRepo, "groups/A/pipeline/synapseHelper.R")

synapseLogin("pangelino@gmail.com", apiKey = "")

projectPath <- synGet("syn2417811")

# Uncomment to create results folder on synapse
# myFolder <- Folder(name = "Merged_Datasets", parentId=projectPath$properties$id)
# myFolder <- synStore(myFolder)

myFolder <- synGet("syn2417826")
message(paste("Uploading datasets"))	

if (! rename.sample){
	upFile <- File(path="./merged.datasets.Rdata", parentId=myFolder$properties$id)
	upFile <- synStore(upFile)
} else {
	upFile <- File(path="./merged.datasets.renamed.Rdata", parentId=myFolder$properties$id)
	upFile <- synStore(upFile)
}
for (n in dataset.list) {
	message(paste("Now processing ", n, "..."))	
	upFile <- File(path=paste0(n,".selected.norm.Rdata"), parentId=myFolder$properties$id)
	upFile <- synStore(upFile)
}

#upFile <- File(path="README", parentId=myFolder$properties$id)
#upFile <- synStore(upFile)



