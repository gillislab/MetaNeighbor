#' Runs MetaNeighbor
#'
#' For each gene set of interest, the function builds a network of rank
#' correlations between all cells. Next,It builds a network of rank correlations
#' between all cells for a gene set. Next, the neighbor voting predictor
#' produces a weighted matrix of predicted labels by performing matrix
#' multiplication between the network and the binary vector indicating cell type
#' membership, then dividing each element by the null predictor (i.e., node
#' degree). That is, each cell is given a score equal to the fraction of its
#' neighbors (including itself), which are part of a given cell type. For
#' cross-validation, we permute through all possible combinations of
#' leave-one-dataset-out cross-validation, and we report how well we can recover
#' cells of the same type as area under the receiver operator characteristic
#' curve (AUROC). This is repeated for all folds of cross-validation, and the
#' mean AUROC across folds is reported. Calls
#' \code{\link{neighborVoting}}.
#'
#' @param dat A SummarizedExperiment object containing gene-by-sample
#' expression matrix.
#' @param i default value 1; non-zero index value of assay containing the matrix
#' data
#' @param experiment_labels A vector that indicates the source/dataset of
#' each sample.
#' @param celltype_labels A character vector or one-hot encoded matrix
#' (cells x cell type) that indicates the cell type of each sample.
#' @param genesets Gene sets of interest provided as a list of vectors.
#' @param bplot default true, beanplot is generated
#' @param fast_version default value FALSE; a boolean flag indicating whether
#' to use the fast and low memory version of MetaNeighbor
#' @param node_degree_normalization default value TRUE; a boolean flag
#' indicating whether to normalize votes by dividing through total node
#' degree.
#' @param batch_size Optimization parameter. Gene sets are processed in groups
#' of size batch_size. The count matrix is first subset to all genes from
#' these groups, then to each gene set individually.
#' @param detailed_results Should the function return the average AUROC across
#' all test datasets (default) or a detailed table with the AUROC for each
#' test dataset?
#' @return A matrix of AUROC scores representing the mean for each gene set
#' tested for each celltype is returned directly
#' (see \code{\link{neighborVoting}}). If detailed_results is set to TRUE,
#' the function returns a table of AUROC scores in each test dataset for each
#' gene set.
#'
#' @seealso \code{\link{neighborVoting}}
#' @examples
#' data("mn_data")
#' data("GOmouse")
#' library(SummarizedExperiment)
#' AUROC_scores = MetaNeighbor(dat = mn_data,
#'                             experiment_labels = as.numeric(factor(mn_data$study_id)),
#'                             celltype_labels = metadata(colData(mn_data))[["cell_labels"]],
#'                             genesets = GOmouse,
#'                             bplot = TRUE)
#' @export
#'

MetaNeighbor <-function(dat, i = 1, experiment_labels, celltype_labels,
                        genesets, bplot = TRUE, fast_version = FALSE,
                        node_degree_normalization = TRUE, batch_size = 10,
                        detailed_results = FALSE) {

    dat <- SummarizedExperiment::assay(dat, i = i)
    
    if (is.vector(celltype_labels)) {
        celltype_labels <- design_matrix(as.character(celltype_labels))
    }

    #check length of experiment_labels equal # of samples
    if(length(experiment_labels) != length(colnames(dat))){
        stop('experiment_labels length does not match number of samples')
    }

    #check length of celltype_labels equal # of samples
    if(length(rownames(celltype_labels)) != length(colnames(dat))){
        stop('celltype_labels length does not match number of samples')
    }

    #check obj contains more than 1 unique study_id
    if(length(unique(experiment_labels)) < 2){
        stop('Found only 1 unique experiment_label. Please use data from more than 1 study!')
    }

    #check genesets matches more than 1 genes in gene_matrix
    genes_in_geneset <- as.character(unlist(genesets))
    genes_in_matrix <- rownames(dat)
    if(length(intersect(genes_in_geneset,genes_in_matrix)) < 1)
        stop('No matching genes between genesets and gene_matrix')

    ROCs              <- vector(mode = "list", length = length(genesets))
    names(ROCs)       <- names(genesets)
    nv_mat            <- matrix(data = 0,
                                ncol = dim(celltype_labels)[2],
                                nrow = length(genesets))
    rownames(nv_mat)  <- names(genesets)
    colnames(nv_mat)  <- colnames(celltype_labels)

    dat <- dat[rownames(dat) %in% unlist(genesets),]
    for (i in 1:ceiling(length(genesets)/batch_size)) {
        batch_start <- batch_size*(i-1)+1
        batch_end <- min(batch_start+batch_size-1, length(genesets))
        subsets <- genesets[batch_start:batch_end]
        subdat <- dat[rownames(dat) %in% unlist(subsets),, drop=FALSE]
        for(l in names(subsets)){
            print(l)
            geneset     <- subsets[[l]]
            m           <- match(rownames(subdat), geneset)
            dat_sub     <- subdat[!is.na(m),, drop=FALSE]
            if (fast_version) {
              ROCs[[l]] <- score_low_mem(
                  dat_sub, experiment_labels, celltype_labels,
                  node_degree_normalization
              )
            } else {
              ROCs[[l]] <- score_default(
                  dat_sub, experiment_labels, celltype_labels,
                  node_degree_normalization
              )
            }
        }
    }
    
    for(i in seq_along(ROCs)){
        nv_mat[i,] <- round(rowMeans(ROCs[[i]], na.rm = TRUE),3)
    }

    if(bplot) {
      plotBPlot(nv_mat)
    }

    if (detailed_results) {
        result <- lapply(ROCs, function(aurocs) {
            aurocs <- tibble::as_tibble(aurocs, rownames = "cell_type")
            tidyr::pivot_longer(aurocs,
                                cols = -cell_type,
                                names_to = "test_dataset",
                                values_to = "auroc",
                                values_drop_na = TRUE)
        })
        result <- dplyr::bind_rows(result, .id = "gene_set")
        return(result)
    } else {
        return(nv_mat)
    }
}

# Compute ROCs according to the default procedure
score_default <- function(dat_sub, experiment_labels, celltype_labels,
                          node_degree_normalization = TRUE) {
  dat_sub     <- stats::cor(as.matrix(dat_sub), method = "s")
  dat_sub     <- as.matrix(dat_sub)
  rank_dat    <- dat_sub
  rank_dat[]  <- rank(dat_sub, ties.method = "average", na.last = "keep")
  rank_dat[is.na(rank_dat)] <- 0
  rank_dat    <- rank_dat/max(rank_dat)
  return(neighborVoting(experiment_labels,
                        celltype_labels,
                        rank_dat,
                        means = FALSE,
                        node_degree_normalization))
}

# Compute ROCs using the approximate low memory version
# For detailed description of vectorized equations, see MetaNeighborUS.R
score_low_mem <- function(dat_sub, study_id, celltype_labels,
                          node_degree_normalization = TRUE) {
  if (nrow(dat_sub) < 2) {
     return(na_matrix(study_id, celltype_labels))  
  }
  dat_sub <- normalize_cols(dat_sub)
  
  # remove cells with no variability
  nonzero_cells <- !matrixStats::colAnyNAs(dat_sub)
  if (sum(nonzero_cells) == 0) {
      return(na_matrix(study_id, celltype_labels))
  }
  dat_sub <- dat_sub[, nonzero_cells, drop=FALSE]
  study_id <- study_id[nonzero_cells]
  celltype_labels <- as.matrix(celltype_labels[nonzero_cells,,drop=FALSE])
  
  unique_study_ids <- unique(study_id)
  aurocs <- c()
  for (study in unique_study_ids) {
    votes <- compute_votes(candidates = dat_sub[, study_id == study, drop=FALSE],
                           voters = dat_sub[, study_id != study, drop=FALSE],
                           voter_id = celltype_labels[study_id != study,, drop=FALSE],
                           node_degree_normalization)
    all_aurocs <- compute_aurocs(
      votes, candidate_id = celltype_labels[study_id == study,, drop = FALSE]
    )
    aurocs <- cbind(aurocs, diag(all_aurocs))
  }
  colnames(aurocs) <- unique_study_ids
  return(aurocs)
}
                                  
# Return AUROC matrix containing only NAs
na_matrix <- function(study_id, celltype_labels) {
  unique_study_ids <- unique(study_id)
  result <- matrix(NA, ncol(celltype_labels), length(unique_study_ids),
                   dimnames = list(colnames(celltype_labels), unique_study_ids))
  return(result)
}

# Compute neighbor voting for a given set of candidates and voters                              
compute_votes <- function(candidates, voters, voter_id,
                          node_degree_normalization = TRUE) {
    votes <- crossprod(candidates, voters %*% voter_id)
    if (node_degree_normalization) {
      vote_names <- dimnames(votes)
      # shift to positive values and normalize node degree
      votes <- matrixStats::t_tx_OP_y(votes, colSums(voter_id), "+")
      node_degree <- colSums(candidates*rowSums(voters)) + ncol(voters)
      votes <- votes / node_degree
      dimnames(votes) <- vote_names
    }
    return(votes)
}
