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

## Probeset IDs to Entrez Gene IDs, hgu133av2
u133a2Map <- function(ids) {
    sapply(AnnotationDbi::mget(as.character(ids), hgu133a2ENTREZID,
                               ifnotfound = NA),
           function(x) paste(x, collapse = "//"))
}

## Entrez Gene IDs + _mt to Entrez Gene IDs (CCLE)
entrezmtMap <- function(ids) {
    sapply(ids, function(x) {strsplit(x, "_mt")[[1]][1]})
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
	out <- tbl$Entrez.GeneID[idxs]
	names(out) <- tbl$ProbesetID[idxs]
	return(out)
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

## TCGA microarray probeset IDs to Entrez Gene IDs
tcga_agilentMap <- function(ids) {
    file <- getFileLocation(synGet("syn2316355"))
    tbl <- read.delim(file, sep = "\t", as.is = TRUE,
                      header = TRUE)
    idxs <- match(ids, tbl$CLID)
	stopifnot(!any(is.na(idxs)))
    genesymbols <- tbl$Gene.Symbol[idxs]	
    entrezID <- rep(NA, length(genesymbols))
    temp <- which(genesymbols != "")
    gstemp <- genesymbols[temp]
    temp2 <- sapply(AnnotationDbi::mget(as.character(gstemp), org.Hs.egSYMBOL2EG,
                                        ifnotfound = NA),
                    function(x) paste(x, collapse = "//"))
    entrezID[temp] <- temp2
    entrezID[entrezID == "NA"] <- NA
	names(entrezID) <- tbl$CLID[idxs]
	return(entrezID)
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

    ## Create the matrix
    if (is.character(raw[, 1])) {
        M <- as.matrix(raw[, -1])
        rownames(M) <- raw[, 1]
    } else {
        M <- as.matrix(raw)
    }

    ## Remove duplicates
    dupls <- which(duplicated(rownames(M)))
    if (length(dupls) > 0) {
        M <- M[-dupls, ]
    }

    return(M)
}

## #####################################################
## Function to average replicates (by Andreas Schlicker)
## #####################################################

averageReplicates = function(mat, loc = c("first", "last")) {
        require(stringr) || stop("Can't load required package \"stringr\"!")
        # Get the cell line name
        loc = match.arg(loc)
        clName = unlist(lapply(str_split(colnames(mat), "_"),
            function(x) { if (loc == "last") el = length(x) else el = 1; x[el] }))
        # Which samples map to which cell line
        replicates = lapply(unique(clName), function(x) { which(clName == x) })
        names(replicates) = unique(clName)
        
        resMat = matrix(NA, nrow = nrow(mat), ncol = length(replicates))
        rownames(resMat) = rownames(mat)
        colnames(resMat) = unlist(lapply(replicates, function(x) { colnames(mat)[x[1]] }))
        for (i in 1:length(replicates)) {
                resMat[, i] = apply(mat[, replicates[[i]], drop = FALSE], 1, mean)
        }
        resMat
}
