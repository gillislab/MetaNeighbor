#' Pretrains model for the unsupervised version of MetaNeighbor
#'
#' When comparing clusters to a large reference dataset, this function summarizes
#' the gene-by-cell matrix into a much smaller highly variable gene-by-cluster matrix
#' which can be fed as training data into MetaNeighborUS, resulting in substantial
#' time and memory savings.
#'
#' @param var_genes vector of high variance genes.
#' @param dat SummarizedExperiment object containing gene-by-sample
#' expression matrix.
#' @param i default value 1; non-zero index value of assay containing the matrix
#' data
#' @param study_id a vector that lists the Study (dataset) ID for each sample
#' @param cell_type a vector that lists the cell type of each sample
#'
#' @return The output is a gene-by-cluster matrix that contains all the information
#' necessary to run MetaNeighborUS from a pre-trained model.
#'
#' @examples
#' data(mn_data)
#' var_genes = variableGenes(dat = mn_data, exp_labels = mn_data$study_id)
#' trained_model = trainModel(var_genes = var_genes,
#'                            dat = mn_data,
#'                            study_id = mn_data$study_id,
#'                            cell_type = mn_data$cell_type))
#' celltype_NV = MetaNeighborUS(trained_model = trained_model,
#'                              dat = mn_data,
#'                              study_id = mn_data$study_id,
#'                              cell_type = mn_data$cell_type)
#' celltype_NV
#'

#' @export
trainModel <- function(var_genes, dat, i = 1, study_id, cell_type) {

    dat    <- SummarizedExperiment::assay(dat, i = i)
    samples <- colnames(dat)
    
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

    study_id <- as.character(study_id)
    cell_type <- as.character(cell_type)

    dat <- normalize_cols(dat)
    label_matrix <- design_matrix(paste(study_id, cell_type, sep = "|"))
    result <- dat %*% label_matrix
    result <- rbind(n_cells = colSums(label_matrix), result)
    return(result)
}
