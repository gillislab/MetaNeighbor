#' Runs the neighbor voting algorithm.
#'
#' The function performs cell type identity prediction based on 'guilt by
#' association' using cross validation. Performance is evaluated by calculating
#' the AUROC for each cell type.
#'
#' @param exp_labels numerical vector that indicates the dataset source of
#' each sample
#' @param cell_labels sample by cell type matrix that indicates the cell type
#' of each sample (0-absent; 1-present)
#' @param network sample by sample adjacency matrix, ranked and standardized
#' between 0-1
#' @param means default \code{TRUE}, determines output formatting
#' @param node_degree_normalization default \code{TRUE}, should predictions
#' be divided by node degree?
#'
#' @return If \code{means = TRUE} (default) a vector containing the mean of
#' AUROC values across cross-validation folds will be returned. If FALSE a list
#' is returned containing a cell type by dataset matrix of AUROC scores, for
#' each fold of cross-validation. Default is over-ridden when more than one cell
#' type is assessed.
#'
#' @examples
#' data("mn_data")
#' data("GOmouse")
#' library(SummarizedExperiment)
#' AUROC_scores = MetaNeighbor(dat = mn_data,
#'                             experiment_labels = as.numeric(factor(mn_data$study_id)),
#'                             celltype_labels = metadata(colData(mn_data))[["cell_labels"]],
#'                             genesets = GOmouse,
#'                             bplot = TRUE)
#' AUROC_scores
#' @seealso \code{\link{MetaNeighbor}}
#' @export
#'

neighborVoting <- function (exp_labels,
                            cell_labels,
                            network,
                            means=TRUE,
                            node_degree_normalization=TRUE){

    # cell_labels : needs to be in 1s and 0s
    x1 <- dim(cell_labels)[2]
    x2 <- dim(cell_labels)[1]
    e <- unique(exp_labels)


    #print("Make genes label CV matrix")
    test_cell_labels <- matrix(cell_labels, nrow=x2, ncol = length(e)*x1)
    exp_cols <- rep(e, each = x1)

    for (i in seq_along(e)){
        d <- which(exp_labels == i)
        a <- which(exp_cols == i)
        test_cell_labels[d,a] <- 0
    }

    #print("Get sums - mat. mul.")
    sum_in  <- (network %*% test_cell_labels)

    #print("Get sums - calc predicts")
    if (node_degree_normalization) {
        #print("Get sums - calc sumall")
        sum_all <- matrix(apply(network, MARGIN = 2, FUN = sum),
                          ncol = dim(sum_in)[2],
                          nrow = dim(sum_in)[1])
        predicts <- sum_in/sum_all
    } else {
        predicts <- sum_in
    }

    #print("Hide training data")
    nans <- which(test_cell_labels == 1, arr.ind = TRUE)
    predicts[nans] <- NA

    #Hide other experiment data
    for (i in seq_along(e)){
        d <- which(exp_labels != i)
        a <- which(exp_cols == i)
        predicts[d,a] <- NA
    }

    #print("Rank test data")
    predicts <- apply(abs(predicts),
                        MARGIN = 2,
                        FUN = rank,
                        na.last = "keep",
                        ties.method = "average")
    filter <- matrix(cell_labels, nrow = x2, ncol = length(e)*x1)

    for (i in seq_along(e)){
        d <- which(exp_labels != i)
        a <- which(exp_cols ==i )
        filter[d,a] <- NA
    }

    negatives <- which(filter == 0, arr.ind = TRUE)
    positives <- which(filter == 1, arr.ind = TRUE)
    predicts[negatives] <- 0

    #print("Calculate ROC - np")
    np <- colSums(filter, na.rm = TRUE) # Postives

    #print("Calculate ROC - nn")
    nn <- apply(filter,
                MARGIN = 2,
                FUN = function(x) sum(x == 0, na.rm = TRUE))     # Negatives

    #print("Calculate ROC - p")
    p  <- apply(predicts, MARGIN = 2, FUN = sum, na.rm = TRUE)

    #print("Calculate ROC - rocN")
    rocNV           <- (p/np - (np+1)/2)/nn
    rocNV           <- matrix(rocNV, ncol = length(e), nrow = x1)
    colnames(rocNV) <- e
    rownames(rocNV) <- colnames(cell_labels)

    if(means==TRUE){
        scores = list(rowMeans(rocNV ,na.rm = TRUE))
    } else {
        scores = list(rocNV)
    }
}
