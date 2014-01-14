library(synapseClient)
library(rGithubClient)

crcRepo <- getRepo("sib-bcf/crcsc")

thisScript <- getPermlink(crcRepo, "groups/A/pipeline/createConfidenceFiles.R")

synapseLogin("charlottesoneson@gmail.com")

## Set the Synapse result directory (phaseI) and the group name
synResultDir <- "syn2274064"
groupName <- "GroupA"

result.files <- list(amc_ajccii = list(id = "syn2299128"),
                     tcga_rnaseq = list(id = "syn2299134"),
                     kfsyscc = list(id = "syn2299124"),
                     french = list(id = "syn2299126"),
                     petacc3 = list(id = "syn2299132"),
                     nki_az = list(id = "syn2299130"),
                     gse2109 = list(id = "syn2299252"),
                     gse4107 = list(id = "syn2299262"),
                     gse4183 = list(id = "syn2299264"),
                     gse10961 = list(id = "syn2299138"),
                     gse15960 = list(id = "syn2299150"),
                     gse13067 = list(id = "syn2299140"),
                     gse13294 = list(id = "syn2299142"),
                     gse20916 = list(id = "syn2299158"),
                     gse23878 = list(id = "syn2299255"),
                     gse37892 = list(id = "syn2299259"),
                     gse17537 = list(id = "syn2299152"),
                     gse17536 = list(id = "syn2313723"),
                     gse14333 = list(id = "syn2299144"),
                     gse8332 = list(id = "syn2313725"),
                     gsk = list(id = "syn2313727"),
                     gse8671 = list(id = "syn2299305"),
                     sanger = list(id = "syn2313729"),
                     agendia_gse42284 = list(id = "syn2299118"),
                     agendia_ico208 = list(id = "syn2299120"),
                     agendia_vhb70 = list(id = "syn2299122"),
                     mdanderson = list(id = "syn2299136"),
                     ccle = list(id = "syn2313731"),
                     tcgacrc_microarray = list(id = "syn2327813"),
                     tcgacrc_merged = list(id = "syn2326831"),
                     tcgacrc_ga = list(id = "syn2326800"),
                     tcgacrc_hiseq = list(id = "syn2326810"))

for (n in names(result.files)) {
    x <- result.files[[n]]

    filename <- getFileLocation(synGet(x$id))
    fln <- strsplit(filename, split = "/")[[1]]
    flname <- fln[length(fln)]

    ## Read the data
    tmp <- read.table(filename, sep = "\t",
                      header = TRUE,
                      as.is = TRUE,
                      fill = TRUE,
                      comment.char = "",
                      check.names = FALSE)

    tmp$maxprobs <- apply(tmp[, c("A", "B", "C", "D", "E")], 1, max)
    tmp$confidence <- c("LOW", "HIGH")[1 + (tmp$maxprobs >= 0.85)]
    tmp$confidence[is.na(tmp$maxprobs)] <- "LOW"

    ## Store in Synapse
    filePath <- file.path(tempdir(), sub(".tsv", "_conf.tsv", flname))
    write.table(tmp[, c("sampleName", "confidence")],
                file = filePath, sep = "\t", quote = FALSE,
                row.names = FALSE, col.names = FALSE)
    synFile <- File(path = filePath, parentId = synResultDir)
    synFile <- synStore(synFile, used = list(
                                     list(entity = x$id, wasExecuted = FALSE),
                                     list(url = thisScript, wasExecuted = TRUE)))
    unlink(filePath)
}

synapseLogout()
