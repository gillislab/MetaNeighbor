# Contains a collection of utility functions
#

#' Make cluster names in format 'study_id|cell_type'
#'
#' @param study_id Character vector containing study ids.
#' @param cell_type Character vector containing cell type names
#  (must be same length as study_id).
#'
#' @return Character vector containing cluster names in the format
#'  study_id|cell_type.
#'
#' @export
makeClusterName <- function(study_id, cell_type) {
  return(paste(standardizeLabel(study_id),
               standardizeLabel(cell_type),
               sep = "|"))
}

#' Remove special characters ("|") from labels to avoid later conflicts
#'
#' @param labels Character vector containing study ids or cell type names.
#' @param replace Special character to replace
#' @param with Character to use instead of special character
#'
#' @return Character vector with replaced special characters.
#'
#' @export
standardizeLabel <- function(labels, replace = "|", with = ".") {
  return(gsub(replace, with, labels, fixed = TRUE))
}

#' Return study ID from a label in format 'study_id|cell_type'
#'
#' @param cluster_name Character vector containing cluster names in the
#' format study_id|cell_type.
#'
#' @return Character vector containing all study ids.
#'
#' @export
getStudyId <- function(cluster_name) {
  return(sapply(strsplit(cluster_name, "|", fixed = TRUE), "[", 1))
}

#' Return cell type from a label in format 'study_id|cell_type'
#'
#' @param cluster_name Character vector containing cluster names in the
#' format study_id|cell_type.
#'
#' @return Character vector containing all cell type names.
#'
#' @export
getCellType <- function(cluster_name) {
  return(sapply(strsplit(cluster_name, "|", fixed = TRUE), "[", 2))
}

# Scale matrix such that all colums sum to 0 and have l2-norm of 1
normalize_cols <- function(M, ranked = TRUE) {
  result <- as.matrix(M)
  if (ranked) {
    result <- matrixStats::colRanks(result, ties.method = "average",
                                    preserveShape = TRUE)
  }
  result <- scale_cols(result)
  dimnames(result) <- dimnames(M)
  return(result)
}

# This is not exactly equivalent to scale (by a factor of sqrt(nrow(M-1)))
# and is a little faster.
# The point is that the dot product of vectors scaled this way is the
# correlation coefficient (which is not true using conventional scaling)
scale_cols <- function(M) {
    cm <- colMeans(M)
    cnorm <- 1 / sqrt(colSums(M**2) - nrow(M) * cm**2)
    matrixStats::t_tx_OP_y(matrixStats::t_tx_OP_y(M, cm, "-"), cnorm, "*")
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

# Compute AUROCs based on a one-vs-best setting
# Inputs are one-vs-rest AUROCs and votes
compute_1v1_aurocs = function(votes, aurocs) {
  result <- matrix(NA, nrow(aurocs), ncol(aurocs), dimnames = dimnames(aurocs))
  for (i in seq_len(ncol(aurocs))) {
    if (all(is.na(aurocs[,i]))) { next }
    top_candidates <- find_top_candidate(votes[,i], aurocs[,i])
    result[top_candidates$best, i] <- top_candidates$score
    result[top_candidates$second, i] <- 1-top_candidates$score
  }
  return(result)
}

# Find best and second best matching clusters
# CAREFUL: the best match is not necessarily the cluster that had highest
# one-vs-rest score!
# We need to consider top 5 candidates and match them up to find best match
find_top_candidate <- function(votes, aurocs) {
  candidates <- extract_top_candidates(aurocs, 5)
  best <- candidates[1]
  votes_best <- votes[names(votes) == best]
  score <- 1
  second_best <- candidates[2]
  for (i in seq(2, length(candidates))) {
    contender <- candidates[i]
    votes_contender <- votes[names(votes) == contender]
    auroc <- c(compute_aurocs(
      as.matrix(c(votes_best, votes_contender)),
      as.matrix(rep(c(1,0), c(length(votes_best), length(votes_contender))))
    ))
    if (auroc < 0.5) {
      second_best <- best
      best <- contender
      score <- 1 - auroc
      votes_best <- votes_contender
    } else if (auroc < score) {
      score <- auroc
      second_best <- contender
    }
  }
  return(list(score = score, best = best, second = second_best))
}

extract_top_candidates <- function(aurocs, n = 10) {
  return(names(head(sort(aurocs, decreasing=TRUE), n = n)))
}

