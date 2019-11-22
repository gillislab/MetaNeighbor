# Contains a collection of utility functions
#

# Scale matrix such that all colums sum to 0 and have l2-norm of 1
normalize_cols <- function(M, ranked = TRUE) {
  M <- as.matrix(M)
  if (ranked) {
    M <- matrixStats::colRanks(M, ties.method = "average", preserveShape = TRUE)
  }
  return(normalize_cols_cpp(M))
}

# Compute AUROCs based on neighbor voting and candidate identities
#
# candidate_id is a binary matrix indicating the cell type of candidates
# If left empty, cell types are assumed to be the row names of the vote matrix.
compute_aurocs <- function(votes, candidate_id = NULL) {
  if (is.null(candidate_id)) {
    positives <- design_matrix(rownames(votes))
  } else {
    positives <- as.matrix(candidate_id)
  }
  n_positives <- colSums(positives)
  n_negatives <- nrow(positives) - n_positives
  sum_of_positive_ranks <- crossprod(
    positives,
    matrixStats::colRanks(votes, ties.method = "average", preserveShape = TRUE)
  )
  colnames(sum_of_positive_ranks) <- colnames(votes)
  result <- (sum_of_positive_ranks / n_positives - (n_positives+1)/2) / n_negatives
  return(result)
}

# Transform a vector with cell_type labels into a binary matrix
design_matrix <- function(cell_type) {
  cell_type <- as.factor(cell_type)
  if (length(levels(cell_type)) > 1) {
    result <- model.matrix(~cell_type-1)
  } else {
    result <- matrix(1, nrow = length(cell_type), ncol = 1)
  }
  colnames(result) <- levels(cell_type)
  return(result)
}

# Return study id from a label in format 'study_id|cell_type'
get_study_id <- function(cluster_name) {
  return(sapply(strsplit(cluster_name, "|", fixed = TRUE), "[", 1))
}

# Return cell type from a label in format 'study_id|cell_type'
get_cell_type <- function(cluster_name) {
  return(sapply(strsplit(cluster_name, "|", fixed = TRUE), "[", 1))
}
