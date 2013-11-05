## File : MetaGenesSubtyping.R
##
## Desc : compute meta gene values, train or apply classifier
## 
## Auth : Sarah Gerster & Charlotte Soneson
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## ---- metagenes-subtyping
##' Assign subtype based on probe sets (Almac DSA) or genes
##' 
##' This function can be used to either 
##' train an LDA classifier or to apply an already trained classifier to assign 
##' subtypes to CRC samples. With data from the Almac DSA, one can work
##' on probe set level. Otherwise the
##' classification is done based on the gene  (Entrez IDs) assignments.
##' 
##' @param X expression matrix, either
##'   \itemize{
##'     \item probe sets x samples: rownames are probe set IDs
##'     \item genes x samples: row names are Entrez Gene IDs. --> set 
##'       \code{probesetlevel} to \code{FALSE}}
##' @param scaling The method to be used for (gene-wise) scaling of 
##'   the expression matrix. One of 
##'   \code{'zscore'}, \code{'mean'}, \code{'median'} or \code{'none'}.
##' @param run.mode The mode to run the function in. One of \code{'train'}
##'   or \code{'classify'}.
##' @param gene.modules A data frame, giving the composition of each module.
##'   Should contain at least two columns, named 'Group' and either
##'   'GeneSymbol' or 'ProbesetID', depending on whether or not the classifier
##'   should be built on probe set-level data (i.e., whether \code{probesetlevel}
##'   is \code{FALSE} or \code{TRUE}.
##' @param subtypes If \code{run.mode} = \code{'train'}, a vector of subtypes, 
##'   of the same length as the number of columns in \code{X}.
##' @param lda_classifier If \code{run.mode} = \code{'classify'}, a classifier
##'   object (such as the one returned from the function in training mode).
##' @param meta.median Logical, if the module scores (meta-genes) should be 
##'   median-centered before training/applying the classifier.
##' @param probesetlevel If TRUE, select original probe sets for each gene used in the classifier.
##' @return If \code{run.mode} = \code{'train'}, a list with
##'   \itemize{
##'     \item \code{XmetaC}, an expression matrix of the normalized meta genes
##'     \item \code{lda_classifier}, the trained classifier}
##' If \code{run.mode} = \code{'classify’}, a list with
##'   \itemize{
##'     \item \code{XmetaC}, an expression matrix of the normalized meta genes
##'     \item \code{subtypes_lda} the subtype assignment for each sample
##'     \item \code{subtypes_posterior}, samples x subtypes -- posterior
##'       probability for each subtype}
##' @export
##' @author Sarah Gerster, Charlotte Soneson
MetaGenesSubtyping <- function(X, scaling, run.mode,
                               gene.modules, subtypes = NULL, 
                               lda_classifier = NULL, 
                               meta.median = FALSE, probesetlevel = TRUE) {
    source("moduleScores.R")

    gene.modules$GeneSymbol <- as.character(gene.modules$GeneSymbol)
    if ("ProbesetID" %in% colnames(gene.modules)) {
        gene.modules$ProbesetID <- as.character(gene.modules$ProbesetID)
    }
    
    ## perform some checks of the input arguments
    if (!(run.mode %in% c("train", "classify"))) {
        stop("run.mode must be one of 'train' or 'classify'")
    }
    if (run.mode == "classify" && is.null(lda_classifier)) {
        stop("You must provide a classifier!")
    }
    if (run.mode == "train" && is.null(subtypes)) {
        stop("You must provide subtypes!")
    }
    
    ## prepare expression matrix
    ## --> genes x samples with Entrez Gene ID as row names
    if (probesetlevel) {
        ## Extract the probesets that are assigned to a module and available in the data set
        SelectedProbesets <- as.character(gene.modules$ProbesetID[gene.modules$ProbesetID %in% rownames(X)])
        ## subset X to only include probesets with a matched Entrez Gene ID
        Xsubtyp <- X[SelectedProbesets, ]
        ## set rownames to Entrez Gene ID
        rownames(Xsubtyp) <- gene.modules$GeneSymbol[match(SelectedProbesets, gene.modules$ProbesetID)]
    } else {
        ## select Entrez Gene IDs that are assigned to a module and available in the data set
        SelectedGenes <- as.character(gene.modules$GeneSymbol[gene.modules$GeneSymbol %in% rownames(X)])
        Xsubtyp <- X[SelectedGenes, ]
    }
    
    ## Scale the variables
    if (scaling == "zscore") {
        Xsubtyp <- t(scale(t(Xsubtyp), center = TRUE, scale = TRUE))
    } else if (scaling == "mean") {
        Xsubtyp <- t(scale(t(Xsubtyp), center = TRUE, scale = FALSE))
    } else if (scaling == "median") {
        Xsubtyp <- sweep(Xsubtyp, FUN = "-", MARGIN = 1, 
                         STATS = apply(Xsubtyp, 1, median))
    } else if (scaling == "none") {
        Xsubtyp <- Xsubtyp
    } else {
        stop("Scaling not defined")
    }
    
    ## calculating meta genes
    ## --> yields a sample x metagene matrix
    Xmeta <- moduleScores(Xsubtyp, 
                          groups = gene.modules[, c(1, 2)], 
                          "median")
    
    ## optionally median-center the meta-genes
    if (meta.median) {
        XmetaC <- apply(Xmeta, 2, function(x) x - median(x))
        stopifnot(all(apply(XmetaC, 2, median) < 1e-10))
    } else {
        XmetaC <- Xmeta
    }
    
    ## train classifier or assign subtypes
    if (run.mode == "train") {
        lda_classifier <- lda(XmetaC, grouping = subtypes)
    } else if (run.mode == "classify") {  
        ldaPredicted <- predict(lda_classifier, newdata = XmetaC)
        subtypes.lda <- ldaPredicted$class
        names(subtypes.lda) <- rownames(ldaPredicted$posterior)
    } 
    
    ## return computed objects as list
    if (run.mode == "train") {
        return(list(XmetaC = XmetaC,
                    lda_classifier = lda_classifier))
    } else if (run.mode == "classify") {
        return(list(XmetaC = XmetaC, 
                    subtypes_lda = subtypes.lda,
                    subtypes_posterior = ldaPredicted$posterior))
    }
}
