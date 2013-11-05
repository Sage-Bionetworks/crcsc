require(MASS)

#' Compute module/meta-gene scores
#'
#' Compute a score for each module and each sample by averaging the expression values of all
#' genes that are part of the module.
#'
#' @param x expression matrix (genes in rows, samples in columns)
#' @param groups data frame with columns GeneSymbol and Groups, containing
#' the module assignment of the genes
#' @param mode how to summarize the expression values across the genes in
#' a module.
#'
#' @return a matrix with module scores (samples in rows, modules in columns)
moduleScores <- function(x, groups, mode = c("mean", "median"))
{
    moduleIDs = setdiff(unique(as.character(groups$Groups)), "0")
    
    temp_matrix = matrix(NA, nrow = ncol(x), 
        ncol = length(moduleIDs))
    
    k = 0
    for (i in moduleIDs) {
        k = k + 1
        selection = match(as.character(groups$GeneSymbol[as.character(groups$Groups) == i]),
            rownames(x))
        if (mode == "median"){
            temp_matrix[, k] = apply(x[selection, ], 2, 
                           FUN = function(x) median(x, na.rm = TRUE))
        }
        if (mode == "mean"){
            temp_matrix[, k] = apply(x[selection, ], 2, 
                           FUN = function(x) mean(x, na.rm = TRUE))
        }
    }
    colnames(temp_matrix) = paste("meta", moduleIDs, sep = "")
    rownames(temp_matrix) = colnames(x)
    
    return(temp_matrix)
}
