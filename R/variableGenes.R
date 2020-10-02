#' Identify a highly variable gene set
#'
#' Identifies genes with high variance compared to their median expression
#' (top quartile) within each experimentCertain function
#'
#' @param dat SummarizedExperiment object containing gene-by-sample
#' expression matrix.
#' @param i default value 1; non-zero index value of assay containing the matrix data
#' @param exp_labels character vector that denotes the source (Study ID) of
#' each sample.
#' @param min_recurrence Number of studies across which a gene must be detected
#' as highly variable to be kept. By default, only genes that are variable
#' across all studies are kept (intersection).
#' @param downsampling_size Downsample each study to downsampling_size
#' samples without replacement. If set to 0 or value exceeds dataset size,
#' no downsampling is applied.
#'
#' @return The output is a vector of gene names that are highly variable in
#' every experiment (intersect)
#'
#' @examples
#' data(mn_data)
#' var_genes = variableGenes(dat = mn_data, exp_labels = mn_data$study_id)
#' var_genes
#'
#' @export
#'

variableGenes <- function(dat, i = 1, exp_labels,
                          min_recurrence = length(unique(exp_labels)),
                          downsampling_size = 10000) {

    dat <- SummarizedExperiment::assay(dat, i = i)
    var_genes <- vector("list")
    experiments <- unique(exp_labels)

    #check length of exp_labels equal # of samples
    if(length(exp_labels) != length(colnames(dat))){
        stop('experiment_labels length does not match number of samples')
    }
    if (min_recurrence > length(experiments)) {
        stop('min_recurrence should be smaller or equal to the number of datasets')
    }

    j <- 1
    for(exp in experiments){
        keep <- which(exp_labels == exp)
        if (downsampling_size > 0 & downsampling_size < length(keep)) {
            keep <- sample(keep, downsampling_size, replace = FALSE)
        }
        data_subset <- as.matrix(dat[, keep])
        median_data <- matrixStats::rowMedians(data_subset)
        variance_data <- MN_rowVars(data_subset)
        quant_med <- unique(stats::quantile(
            median_data,
            probs = seq(from = 0, to = 1, length = 11),
            type = 5
        ))
        genes_list <- vector("list", length = length(quant_med))

        for(i in seq_along(quant_med)){
            if(i == 1){
                filt1     <- median_data <= quant_med[i]
            } else {
                filt1     <- median_data <= quant_med[i] & median_data > quant_med[i-1]
            }
            var_temp  <- variance_data[filt1]
            quant_var <- stats::quantile(var_temp, na.rm = TRUE)
            filt2     <- var_temp > quant_var[4]
            genes_list[[i]] <- names(var_temp)[filt2]
        }
        temp <- length(genes_list)
        var_genes[[j]] <- unlist(genes_list[1:temp-1])
        j <- j+1
    }
    pooled_genes <- factor(unlist(var_genes))
    recurrence <- tabulate(pooled_genes)
    result <- levels(pooled_genes)[recurrence >= min_recurrence]
    return(result)
}

MN_rowVars <- function(M) {
    if (is(M, "dgCMatrix")) {
        M <- Matrix::t(M)
        result <- Matrix::colMeans(M**2) - Matrix::colMeans(M)**2
        result <- result * nrow(M) / (nrow(M)-1)
    } else {
        result <- matrixStats::rowVars(as.matrix(M))
        names(result) <- rownames(M)
    }
    return(result)
}
