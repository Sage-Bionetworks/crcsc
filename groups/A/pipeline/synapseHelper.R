##
## Script to help with data handling via Synapse
##
## Author: Charlotte Soneson, based on scripts by Justin Guinney and Andreas Schlicker
######################################################################################

## Probeset IDs to Entrez Gene IDs, hgu133plus2
u133plus2Map <- function(ids) {
    sapply(AnnotationDbi::mget(as.character(ids), hgu133plus2ENTREZID,
                               ifnotfound = NA),
           function(x) paste(x, collapse = "//"))
}

## Gene symbols to Entrez Gene IDs
symbolMap <- function(ids) {
    sapply(AnnotationDbi::mget(as.character(ids), org.Hs.egSYMBOL2EG,
                               ifnotfound = NA),
           function(x) paste(x, collapse = "//"))
}

## PETACC-3 probeset IDs to Entrez Gene IDs
petaccMap <- function(ids) {
    file <- getFileLocation(synGet("syn2199825"))
    tbl <- read.table(file, sep = "\t", as.is = TRUE,
                      header = TRUE, comment = "")
    idxs <- match(ids, tbl$ProbesetID)
    stopifnot(!any(is.na(idxs)))
    return(tbl$Entrez.GeneID[idxs])
}

## DiscoverPrint_19742 probeset IDs to Entrez Gene IDs
discoverprint_19742Map <- function(ids) {
    file <- getFileLocation(synGet("syn2192791"))
    tbl <- read.table(file, sep = "\t", as.is = TRUE,
                      header = TRUE)
    idxs <- match(ids, tbl$probe_id)
    stopifnot(!any(is.na(idxs)))
    refseq <- tbl$systematic_name[idxs]
    return(sapply(mget(as.character(refseq), org.Hs.egREFSEQ2EG, ifnotfound = NA),
           function(x) paste(x, collapse = "//")))
}

## DiscoverPrint_32627 probeset IDs to Entrez Gene IDs
discoverprint_32627Map <- function(ids) {
    file <- getFileLocation(synGet("syn2192801"))
    tbl <- read.table(file, sep = "\t", as.is = TRUE,
                      header = TRUE, comment = "", quote = "")
    idxs <- match(ids, tbl$ProbeName)
    stopifnot(!any(is.na(idxs)))
    refseq <- tbl$SystematicName[idxs]
    return(sapply(mget(as.character(refseq), org.Hs.egREFSEQ2EG, ifnotfound = NA),
           function(x) paste(x, collapse = "//")))
}

## #####################################################
## Function to load data
## #####################################################

loadMatrix <- function(synId, file = "", sep = "\t", quote = "", ...) {
    filename <- getFileLocation(synGet(synId))

    ## Unzip if necessary
    if (length(grep("zip$", filename)) > 0) {
        temp <- unz(filename, file)
    } else {
        temp <- filename
    }

    ## Read the data
    raw <- read.table(temp, sep = sep,
                      header = TRUE, quote = quote,
                      as.is = TRUE,
                      fill = TRUE,
                      comment.char = "",
                      check.names = FALSE)

    ## Remove duplicates
    dupls <- which(duplicated(raw[, 1]))
    if (length(dupls) > 0) {
        raw <- raw[-dupls, ]
    }

    ## Create the matrix
    if (is.character(raw[, 1])) {
        M <- as.matrix(raw[, -1])
        rownames(M) <- raw[, 1]
    } else {
        M <- as.matrix(raw)
    }

    return(M)
}
