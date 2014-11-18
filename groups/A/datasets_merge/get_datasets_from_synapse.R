# 
# Author: p angelino
###############################################################################

if (!exists("work.dir")) {
	work.dir <- "/export/scratch/paolo/SageCRCsubtyping/analysis/work/"
}

setwd(work.dir)

library(synapseClient)
library(rGithubClient)
library(hgu133plus2.db)
library(hgu133a2.db)
library(org.Hs.eg.db)

crcRepo <- getRepo("sib-bcf/crcsc")

#sourceRepoFile(crcRepo, "groups/A/pipeline/synapseHelper.R")
source("/export/scratch/paolo/workspace/Rworks/crcsc/groups/A/pipeline/synapseHelper.R")
helperFunctions <- getPermlink(crcRepo, "groups/A/pipeline/synapseHelper.R")

synapseLogin("pangelino@gmail.com", apiKey = "")


## Define the datasets to work with
exprList <- list(agendia_gse42284 = list(synId = "syn2192792",
				entrezMapper = discoverprint_19742Map, 
				file = "GSE42284_normalized_data_matrix.txt"),
		agendia_ico208 = list(synId = "syn2192796",
				entrezMapper = discoverprint_19742Map,
				file = "ICO208_normalized_data.txt"),
		agendia_vhb70 = list(synId = "syn2192799",
				entrezMapper = discoverprint_32627Map,
				file = "VHB70_normalized_data.txt"),
		amc_ajccii = list(synId = "syn2363559",
				entrezMapper = u133plus2Map),
		french = list(synId = "syn2363561",
				entrezMapper = u133plus2Map),
		kfsyscc = list(synId = "syn2363564",
				entrezMapper = u133plus2Map),
		mdanderson = list(synId = "syn2233387",
				entrezMapper = discoverprint_32627Map),
		nki_az = list(synId = "syn2363566",
				entrezMapper = u133plus2Map),
		petacc3 = list(synId = "syn2175581",
				entrezMapper = petaccMap),
		tcgacrc_ga = list(synId = "syn2326094",
				entrezMapper = symbolMap),
		tcgacrc_hiseq = list(synId = "syn2326100",
				entrezMapper = symbolMap),
		tcgacrc_merged = list(synId = "syn2325328",
				entrezMapper = symbolMap),
		tcgacrc_microarray = list(synId = "syn2316354",
				quote = '""', 
				entrezMapper = tcga_agilentMap),		
		tcga_rnaseq = list(synId = "syn2161141",
				entrezMapper = symbolMap),
		gse10961 = list(synId = "syn2361860",
				entrezMapper = u133plus2Map),		
		gse13067 = list(synId = "syn2361866",
				entrezMapper = u133plus2Map),
		gse13294 = list(synId = "syn2361870",
				entrezMapper = u133plus2Map),
		gse14333 = list(synId = "syn2361878",
				entrezMapper = u133plus2Map),
		gse15960 = list(synId = "syn2361882",
				entrezMapper = u133plus2Map),
		gse17537 = list(synId = "syn2361889",
				entrezMapper = u133plus2Map),		
		gse20916 = list(synId = "syn2362291",
				entrezMapper = u133plus2Map),
		gse2109 = list(synId = "syn2362299",
				entrezMapper = u133plus2Map),
		gse23878 = list(synId = "syn2362301",
				entrezMapper = u133plus2Map),
		gse37892 = list(synId = "syn2362305",
				entrezMapper = u133plus2Map),
		gse4107 = list(synId = "syn2362307",
				entrezMapper = u133plus2Map),
		gse4183 = list(synId = "syn2362311",
				entrezMapper = u133plus2Map),
		gse8671 = list(synId = "syn2362317",
				entrezMapper = u133plus2Map),
		gse17536 = list(synId = "syn2361887",
				entrezMapper = u133plus2Map))

		
		## Reading relaxed outliers
		outlierFileList <- list(agendia_gse42284 = "syn2374731",
				agendia_ico208 = "syn2374735",
				agendia_vhb70 = "syn2374739",
				amc_ajccii = "syn2374700",
				french = "syn2374704",
				kfsyscc = "syn2374708",
				mdanderson = "syn2374712",
				nki_az = "",
				petacc3 = "syn2374720",
				tcgacrc_ga = "",
				tcgacrc_hiseq = "",
				tcgacrc_merged = "syn2374724",
				tcgacrc_microarray = "",		
				tcga_rnaseq = "",
				gse10961 = "",		
				gse13067 = "syn2373932",
				gse13294 = "",
				gse14333 = "",
				gse15960 = "syn2374643",
				gse17537 = "syn2374653",		
				gse20916 = "",
				gse2109 = "syn2374663",
				gse23878 = "",
				gse37892 = "syn2374677",
				gse4107 = "syn2374683",
				gse4183 = "",
				gse8671 = "",
				gse17536 = "syn2374649")

		
# datasets selected for merging		
dataset.list <- c("agendia_gse42284", "agendia_ico208",  "agendia_vhb70",    "amc_ajccii",      
				"french",           "kfsyscc",          "mdanderson",       "nki_az",          
				"petacc3",          "tcgacrc_merged",   "gse10961",		"gse13067",         "gse13294",         
				"gse14333",        
				"gse15960",         "gse17537",         "gse20916",         "gse2109" ,        
				"gse23878",         "gse37892",         "gse4107",          "gse4183" ,        
				"gse8671",          "gse17536")
exprList <- exprList[dataset.list]		

# save dataset list
save(file=paste0(work.dir,"dataset.list.Rdata"), dataset.list)

# loading datasets		
for (n in names(exprList)) {
	tryCatch({
				x = exprList[[n]]
				message(paste("Now processing ", n, "..."))
				
				sep <- ifelse (is.null(x$sep), "\t", x$sep)
				quote <- ifelse (is.null(x$quote), "", x$quote)
				file <- ifelse (is.null(x$file), "", x$file)
				
				## Load data -> returns a data matrix (features x samples)
				exprList[[n]][["data.mat"]] <- loadMatrix(x$synId, file = file, sep = sep, quote = quote)
				if (!outlierFileList[[n]]=="") {
					outlierFile <- synGet(outlierFileList[[n]])
					outliers <- read.delim(file=outlierFile@filePath, header=F, colClasses = c("character"))
					outliers.id <- which(colnames(exprList[[n]][["data.mat"]]) %in% outliers$V1)
					exprList[[n]][["data.mat"]]
					exprList[[n]][["data.mat"]] <- exprList[[n]][["data.mat"]][,-outliers.id] 
				}				
				
			})
	
}

#synapseLogout()

