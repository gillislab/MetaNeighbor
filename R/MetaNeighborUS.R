#' Runs unsupervised version of MetaNeighbor
#'
#' When it is difficult to know how cell type labels compare across datasets this
#' function helps users to make an educated guess about the overlaps without
#' requiring in-depth knowledge of marker genes
#'
#' @param var_genes vector of high variance genes.
#' @param dat SummarizedExperiment object containing gene-by-sample
#' expression matrix.
#' @param i default value 1; non-zero index value of assay containing the matrix
#' data
#' @param study_id a vector that lists the Study (dataset) ID for each sample
#' @param cell_type a vector that lists the cell type of each sample
#' @param trained_model default value NULL; a matrix containing a trained model
#' generated from MetaNeighbor::trainModel. If not NULL, the trained model is
#' treated as training data and dat is treated as testing data. If a trained model
#' is provided, fast_version will automatically be set to TRUE and var_genes will
#' be overridden with genes used to generate the trained_model
#' @param fast_version default value FALSE; a boolean flag indicating whether
#' to use the fast and low memory version of MetaNeighbor
#' @param node_degree_normalization default value TRUE; a boolean flag indicating
#' whether to use normalize votes by dividing through total node degree.
#' @param one_vs_best default value FALSE; a boolean flag indicating whether
#' to compute AUROCs based on a best match against second best match setting
#' (default version is one-vs-rest). This option is currently only relevant
#' when fast_version = TRUE.
#' @param symmetric_output default value TRUE; a boolean flag indicating whether
#' to average AUROCs in the output matrix.
#'
#' @return The output is a cell type-by-cell type mean AUROC matrix, which is
#' built by treating each pair of cell types as testing and training data for
#' MetaNeighbor, then taking the average AUROC for each pair (NB scores will not
#' be identical because each test cell type is scored out of its own dataset,
#' and the differential heterogeneity of datasets will influence scores).
#' If symmetric_output is set to FALSE, the training cell types are displayed
#' as columns and the test cell types are displayed as rows.
#' If trained_model was provided, the output will be a cell type-by-cell
#' type AUROC matrix with training cell types as columns and test cell types
#' as rows (no swapping of test and train, no averaging).
#'
#' @examples
#' data(mn_data)
#' var_genes = variableGenes(dat = mn_data, exp_labels = mn_data$study_id)
#' celltype_NV = MetaNeighborUS(var_genes = var_genes,
#'                              dat = mn_data,
#'                              study_id = mn_data$study_id,
#'                              cell_type = mn_data$cell_type)
#' celltype_NV
#'

#' @export
MetaNeighborUS <- function(var_genes = c(), dat, i = 1, study_id, cell_type,
                           trained_model = NULL, fast_version = FALSE,
                           node_degree_normalization = TRUE, one_vs_best = FALSE,
                           symmetric_output = TRUE) {

    dat    <- SummarizedExperiment::assay(dat, i = i)
    samples <- colnames(dat)
    if (!is.null(trained_model)) {
        trained_model <- as.matrix(trained_model)
        var_genes <- rownames(trained_model)[-1]
    }
    
    #check obj contains study_id
    if(length(study_id)!=length(samples)){
        stop('study_id length does not match number of samples')
    }

    #check obj contains cell_type
    if(length(cell_type)!=length(samples)){
        stop('cell_type length does not match number of samples')
    }

    matching_vargenes <- match(rownames(dat), var_genes)
    matching_vargenes_count   <- sum(!is.na(matching_vargenes))

    if(matching_vargenes_count < 2){
        stop("matching_vargenes should have more than 1 matching genes!",
             call. = TRUE)
    } else if(matching_vargenes_count < 5) {
        warning("matching_vargenes should have more matching genes!",
                immediate. = TRUE)
    }
    dat <- dat[!is.na(matching_vargenes),]
    if (!is.null(trained_model)) {
        trained_model <- trained_model[c("n_cells", rownames(dat)),]
    }

    study_id <- as.character(study_id)
    cell_type <- as.character(cell_type)

    if (is.null(trained_model)) {
        if (fast_version) {
          cell_NV <- MetaNeighborUSLowMem(dat, study_id, cell_type,
                                          node_degree_normalization, one_vs_best)
        } else {
          cell_NV <- MetaNeighborUSDefault(dat, study_id, cell_type,
                                           node_degree_normalization)
        }
        if (symmetric_output) {
            cell_NV <- (cell_NV+t(cell_NV))/2
        }
    } else {
        cell_NV <-  MetaNeighborUS_from_trained(trained_model, dat, study_id, cell_type,
                                                node_degree_normalization, one_vs_best)
    }
    return(cell_NV)
}

MetaNeighborUSDefault <- function(dat, study_id, cell_type, node_degree_normalization = TRUE) {
    dat <- as.matrix(dat)
    pheno <- as.data.frame(cbind(study_id,cell_type), stringsAsFactors = FALSE)
    pheno$StudyID_CT <- paste(pheno$study_id, pheno$cell_type, sep = "|")
    celltypes   <- unique(pheno$StudyID_CT)
    cell_labels <- matrix(0, ncol=length(celltypes), nrow=dim(pheno)[1])
    rownames(cell_labels) <-colnames(dat)
    colnames(cell_labels) <- celltypes

    for(i in seq_along(celltypes)){
        type <- celltypes[i]
        matching_celltype <- match(pheno$StudyID_CT, type)
        cell_labels[!is.na(matching_celltype),i]  <- 1
    }

    cor_data    <- stats::cor(dat, method="s")
    rank_data   <- cor_data*0
    rank_data[] <- rank(cor_data, ties.method = "average", na.last = "keep")
    rank_data[is.na(rank_data)] <- 0
    rank_data   <- rank_data/max(rank_data)
    sum_in      <- (rank_data) %*% cell_labels
    
    if (node_degree_normalization) {
        sum_all     <- matrix(apply(rank_data, MARGIN = 2, FUN = sum),
                              ncol = dim(sum_in)[2],
                              nrow = dim(sum_in)[1])
        predicts    <- sum_in/sum_all
    } else {
        predicts <- sum_in        
    }

    cell_NV     <- matrix(0, ncol=length(celltypes), nrow=length(celltypes))
    colnames(cell_NV) <- colnames(cell_labels)
    rownames(cell_NV) <- colnames(cell_labels)

    for(i in seq_len(dim(cell_labels)[2])){
        predicts_temp <- predicts

        matching_celltype <- match(pheno$StudyID_CT, colnames(cell_labels)[i])
        unique_studyID    <- unique(pheno[!is.na(matching_celltype),"study_id"])
        matching_studyID  <- match(pheno$study_id, unique_studyID)
        pheno2            <- pheno[!is.na(matching_studyID),]
        predicts_temp     <- predicts_temp[!is.na(matching_studyID),]
        predicts_temp     <- apply(abs(predicts_temp),
                                   MARGIN = 2,
                                   FUN = rank,
                                   na.last= "keep",
                                   ties.method="average")


        filter  <- matrix(0, ncol=length(celltypes), nrow=dim(pheno2)[1])
        matches <- match(pheno2$StudyID_CT, colnames(cell_labels)[i])
        filter[!is.na(matches),seq_along(celltypes)] <- 1

        negatives = which(filter == 0, arr.ind = TRUE)
        positives = which(filter == 1, arr.ind = TRUE)

        predicts_temp[negatives] <- 0

        np <- colSums(filter, na.rm = TRUE)
        nn <- apply(filter, MARGIN = 2, FUN = function(x) sum(x==0,na.rm=TRUE))
        p  <- apply(predicts_temp, MARGIN = 2, FUN = sum, na.rm = TRUE)

        cell_NV[i,]= (p/np - (np+1)/2)/nn
    }
    return(cell_NV)
}

# The fast version is vectorized according to the following equations
# (Note that the point of these equations is to *never* compute the cell-cell network
#  by reordering the matrix operations):
#  - INPUTS:
#    + Q = test (Query) data (genes x cells)
#    + R = train (Ref) data (genes x cells)
#    + L = binary encoding of training cell types (Labels) (cells x cell types)
#    + S = binary encoding of train Studies (cells x studies)
#  - NOTATIONS:
#    + X* = normalize_cols(X) ~ scale(colRanks(X)) denotes normalized data
#           (Spearman correlation becomes a simple dot product on normalized data)
#    + N = Spearman(Q,R) = t(Q*).R* is the cell-cell similarity network
#    + CL = R*.L are the cell type centroids (in the normalized space)
#    + CS = R*.S are the study centroids (in the normalized space)
#    + 1.L = colSums(L) = number of cells per (train) cell type
#    + 1.S = colSums(S) = number of cells per (train) study
#  - WITHOUT node degree normalization
#    + Votes = N.L = t(Q*).R*.L = t(Q*).CL
#  - WITH node degree normalization
#    + Network becomes N+1 to avoid negative values
#    + Votes = (N+1).L = N.L + 1.L = t(Q*).CL + 1.L
#    + Node degree = (N+1).S = t(Q*).CS + 1.S
#    + Note: Node degree is computed independently for each train study.
#
MetaNeighborUSLowMem <- function(dat, study_id, cell_type,
                                 node_degree_normalization = TRUE, one_vs_best = FALSE) {
  dat <- normalize_cols(dat)
  label_matrix <- design_matrix(paste(study_id, cell_type, sep = "|"))
  cluster_centroids <- dat %*% label_matrix
  n_cells_per_cluster <- colSums(label_matrix)
    
  result <- predict_and_score(dat, study_id, cell_type,
                              cluster_centroids, n_cells_per_cluster,
                              node_degree_normalization, one_vs_best)
  result <- result[, rownames(result)]
  return(result)
}
         
# WARNING: function assumes that data have been normalized with normalize_cols
predict_and_score <- function(dat, study_id, cell_type,
                              cluster_centroids, n_cells_per_cluster,
                              node_degree_normalization = TRUE, one_vs_best = FALSE) {
  colnames(dat) <- paste(study_id, cell_type, sep = "|")
  if (node_degree_normalization) {
    centroid_study_label <- getStudyId(colnames(cluster_centroids))
    study_matrix <- design_matrix(centroid_study_label)
    study_centroids <- cluster_centroids %*% study_matrix
    n_cells_per_study <- n_cells_per_cluster %*% study_matrix
    train_study_id <- colnames(study_matrix)
  }

  result <- c()
  for (test_study in unique(study_id)) {
    is_test <- study_id == test_study
    test_dat <- dat[, is_test]
    votes <- crossprod(test_dat, cluster_centroids)
      
    if (node_degree_normalization) {
      # shift to positive values and normalize node degree
      votes <- sweep(votes, 2, n_cells_per_cluster, FUN = "+")
      node_degree <- crossprod(test_dat, study_centroids)
      node_degree <- sweep(node_degree, 2, n_cells_per_study, FUN = "+")
      for (train_study in unique(train_study_id)) {
        is_train <- centroid_study_label == train_study
        votes[, is_train] <- votes[, is_train] / node_degree[, train_study]
      }
    }
      
    aurocs <- compute_aurocs(votes, design_matrix(rownames(votes)))
    if (one_vs_best) {
      result <- rbind(result, compute_1v1_aurocs(votes, aurocs))
    } else {
      result <- rbind(result, aurocs)
    }
  }
  return(result)
}

MetaNeighborUS_from_trained <- function(trained_model, test_dat, study_id, cell_type,
                                        node_degree_normalization = TRUE, one_vs_best = FALSE) {
  dat <- normalize_cols(test_dat)
  cluster_centroids <- trained_model[-1,]
  n_cells_per_cluster <- trained_model[1,]
    
  result <- predict_and_score(dat, study_id, cell_type,
                              cluster_centroids, n_cells_per_cluster,
                              node_degree_normalization, one_vs_best)
  return(result)
}
