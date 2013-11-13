##
## Apply the module-based LDA to all data sets in the CRC subtyping consortium
##
## Author: Charlotte Soneson
## ############################################################################

library(synapseClient)
library(rGithubClient)
library(hgu133plus2.db)
library(org.Hs.eg.db)

crcRepo <- getRepo("sib-bcf/crcsc")

sourceRepoFile(crcRepo, "groups/A/pipeline/moduleScores.R")
ModuleFunctions <- getPermlink(crcRepo, "groups/A/pipeline/moduleScores.R")
sourceRepoFile(crcRepo, "groups/A/pipeline/synapseHelper.R")
helperFunctions <- getPermlink(crcRepo, "groups/A/pipeline/synapseHelper.R")
sourceRepoFile(crcRepo, "groups/A/pipeline/MetaGenesSubtyping.R")
subtypingFunctions <- getPermlink(crcRepo, "groups/A/pipeline/MetaGenesSubtyping.R")

thisScript <- getPermlink(crcRepo, "groups/A/pipeline/groupA_subtyping_pipeline.R")

synapseLogin("charlottesoneson@gmail.com")

## Set the Synapse result directory (phaseI) and the group name
synResultDir <- "syn2229267"
groupName <- "GroupA"

## Load the module definitions and the LDA classifier
modulesSynId <- "syn2293486"
modules <- read.table(getFileLocation(synGet(modulesSynId)), header = TRUE,
                      sep = "\t", stringsAsFactors = FALSE)
classifierSynId <- "syn2293487"
load(getFileLocation(synGet(classifierSynId)))

## Define the data sets to be assigned subtypes
exprList <- list(agendia_gse42284 = list(synId = "syn2192792",
                     entrezMapper = discoverprint_19742Map, 
                     file = "GSE42284_normalized_data_matrix.txt"),
                 agendia_ico208 = list(synId = "syn2192796",
                     entrezMapper = discoverprint_19742Map,
                     file = "ICO208_normalized_data.txt"),
                 agendia_vhb70 = list(synId = "syn2192799",
                     entrezMapper = discoverprint_32627Map,
                     file = "VHB70_normalized_data.txt"),
                 kfsyscc = list(synId = "syn2169565",
                     entrezMapper = u133plus2Map),
                 french = list(synId = "syn2171434",
                     entrezMapper = u133plus2Map),
                 amc_ajccii = list(synId = "syn2159423",
                     entrezMapper = u133plus2Map,
                     sep = ",", quote = "\""),
                 nki_az = list(synId = "syn2176657",
                     entrezMapper = u133plus2Map),
                 petacc3 = list(synId = "syn2175581",
                     entrezMapper = petaccMap),
                 tcga_rnaseq = list(synId = "syn2161141",
                     entrezMapper = symbolMap),
                 mdanderson = list(synId = "syn2233387",
                     entrezMapper = discoverprint_32627Map),
                 gse10961 = list(synId = "syn2177194",
                     entrezMapper = u133plus2Map),
                 gse13067 = list(synId = "syn2177888",
                     entrezMapper = u133plus2Map),
                 gse13294 = list(synId = "syn2177894",
                     entrezMapper = u133plus2Map),
                 gse14333 = list(synId = "syn2181079",
                     entrezMapper = u133plus2Map),
                 gse15960 = list(synId = "syn2177199",
                     entrezMapper = u133plus2Map),
                 gse17537 = list(synId = "syn2178128",
                     entrezMapper = u133plus2Map),
                 gse20916 = list(synId = "syn2177899",
                     entrezMapper = u133plus2Map),
                 gse2109 = list(synId = "syn2177169",
                     entrezMapper = u133plus2Map),
                 gse23878 = list(synId = "syn2177902",
                     entrezMapper = u133plus2Map),
                 gse37892 = list(synId = "syn2178082",
                     entrezMapper = u133plus2Map),
                 gse4107 = list(synId = "syn2177179",
                     entrezMapper = u133plus2Map),
                 gse4183 = list(synId = "syn2177187",
                     entrezMapper = u133plus2Map),
                 gse8671 = list(synId = "syn2181088",
                     entrezMapper = u133plus2Map))

for (n in names(exprList)) {
    x = exprList[[n]]
    message(paste("Now processing ", n, "..."))

    sep <- ifelse (is.null(x$sep), "\t", x$sep)
    quote <- ifelse (is.null(x$quote), "", x$quote)
    file <- ifelse (is.null(x$file), "", x$file)

    ## Load data -> returns a data matrix (features x samples)
    data.mat <- loadMatrix(x$synId, file = file, sep = sep, quote = quote)

    if (n != "petacc3") {
        ## Get the Entrez Gene ID for each row of the matrix
        entrezID <- x$entrezMapper(rownames(data.mat))
        
        ## Keep only rows corresponding to Entrez ID
        keep.idx <- which(!is.na(entrezID))
        data.mat <- data.frame(data.mat[keep.idx, ])
        entrezID <- entrezID[keep.idx]
        
        ## Summarize on Entrez Gene ID level (select feature with maximal MAD)
        data.mat <- do.call(rbind, lapply(split(data.mat, f = entrezID),
                                          function(x) {
                                              x[which.max(apply(x, 1, mad,
                                                                na.rm = TRUE)), , drop = FALSE]
                                          }))
    }

    ## Apply classifier
    results <- MetaGenesSubtyping(data.mat, scaling = "zscore",
                                  run.mode = "classify", gene.modules = modules,
                                  subtypes = NULL,
                                  lda_classifier = lda_classifier,
                                  meta.median = FALSE,
                                  probesetlevel = ifelse(n == "petacc3", TRUE, FALSE))

    tmp <- data.frame(sampleName = rownames(results$subtypes_posterior),
                      round(results$subtypes_posterior, digits = 4))

    ## Save
    filePath <- file.path(tempdir(), paste(groupName, "_", x$synId, "_", n, ".tsv", sep = ""))
    write.table(tmp, file = filePath, sep = "\t", quote = FALSE, row.names = FALSE)

    ## Store in Synapse
    synFile <- File(path = filePath, parentId = synResultDir)
    synFile <- synStore(synFile, used = list(
                                     list(entity = "syn2175581", wasExecuted = FALSE),
                                     list(entity = x$synId, wasExecuted = FALSE),
                                     list(url = ModuleFunctions, wasExecuted = FALSE),
                                     list(url = helperFunctions, wasExecuted = FALSE),
                                     list(url = thisScript, wasExecuted = TRUE),
                                     list(entity = modulesSynId, wasExecuted = FALSE),
                                     list(entity = classifierSynId, wasExecuted = FALSE)))
    unlink(filePath)
}

synapseLogout()
