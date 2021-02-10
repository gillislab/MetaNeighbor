
#' Merge multiple SingleCellExperiment objects.
#'
#' @param sce_list A *named* list, where values are SingleCellExperiment objects
#' and names are SingleCellExperiment objects.
#'
#' @return A SingleCellExperiment object containing the input datasets with the
#' following limitations: (i) only genes common to all datasets are kept,
#' (ii) only colData columns common to all datasets are kept,
#' (iii) only assays common to all datasets (i.e., having the same name)
#' are kept, (iv) all other slots (e.g., reducedDims or rowData) will
#' be ignored and left empty.
#' The SingleCellExperiment object contains a "study_id" column,
#' mapping each cell to its original dataset (names in "sce_list").
#'
#' @export
mergeSCE <- function(sce_list) {
    # test if the list has names for the datasets
    if (is.null(names(sce_list))) {
        stop("The list of SCE datasets must be named (e.g. list(dataset_1=..., dataset2=...)).")
    }
    
    # gather expression data
    common_assays <- Reduce(
        intersect,
        lapply(sce_list, function(sce) {
            names(SummarizedExperiment::assays(sce))
        })
    )
    common_genes <- Reduce(intersect, lapply(sce_list, rownames))
    new_assays <- lapply(common_assays, aggregate_assay, sce_list, common_genes)
    names(new_assays) <- common_assays
    
    # gather common coldata
    common_coldata <- Reduce(
        intersect,
        lapply(sce_list, function(sce) {
            names(SingleCellExperiment::colData(sce))
        })
    )
    new_coldata <- lapply(sce_list, function(sce) {
        SingleCellExperiment::colData(sce)[, common_coldata, drop=FALSE]
    })
    new_coldata <- do.call(rbind, new_coldata)
    gc()
    
    # create fused SCE object
    result <- SingleCellExperiment::SingleCellExperiment(assays = new_assays,
                                                         colData = new_coldata)
    result$study_id <- rep(names(sce_list), sapply(sce_list, ncol))
    return(result)
}

# Concatenate count matrices from multiple SingleCellExperiment objects.
aggregate_assay <- function(assay_name, sce_list, common_genes) {
    result <- lapply(sce_list, function(sce) {
        SummarizedExperiment::assays(sce)[[assay_name]][common_genes,]
    })
    gc()
    result <- do.call(cbind, result)
    gc()
    return(result)
}