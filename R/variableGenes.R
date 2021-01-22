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
    experiments <- unique(exp_labels)

    #check length of exp_labels equal # of samples
    if(length(exp_labels) != length(colnames(dat))){
        stop('experiment_labels length does not match number of samples')
    }
    if (min_recurrence > length(experiments)) {
        stop('min_recurrence should be smaller or equal to the number of datasets')
    }

    gene_stats <- list()
    for(exp in experiments){
        keep <- which(exp_labels == exp)
        if (downsampling_size > 0 & downsampling_size < length(keep)) {
            keep <- sample(keep, downsampling_size, replace = FALSE)
        }
        gene_stats[[exp]] <- variable_genes_single_exp(dat[, keep])
    }
    gene_stats <- dplyr::bind_rows(gene_stats, .id = "study_id")
    `%>%` <- dplyr::`%>%`
    gene_stats <- gene_stats %>%
        dplyr::group_by(gene) %>%
        dplyr::summarize(recurrence = sum(is_hvg),
                         score = mean(var_quant)) %>%
        dplyr::filter(recurrence >= min_recurrence) %>%
        dplyr::arrange(desc(recurrence), desc(score))
    result <- dplyr::pull(gene_stats, gene)
    return(result)
}

variable_genes_single_exp = function(data_subset) {
    data_subset <- as.matrix(data_subset)
    variance_data <- MN_rowVars(data_subset)
    median_data <- matrixStats::rowMedians(data_subset)
    quant_med <- unique(stats::quantile(
        median_data, probs = seq(0, 1, length = 11), type = 5
    ))
    
    # assign genes to 5 quantile bins based on median expression
    # remove bin with high expressing genes
    `%>%` <- dplyr::`%>%`
    result <- data.frame(gene = names(variance_data),
                         variance = variance_data,
                         bin_med = cut(median_data, c(-1,quant_med))) %>%
        dplyr::filter(bin_med != levels(bin_med)[length(levels(bin_med))]) %>%
        dplyr::group_by(bin_med) %>%
        dplyr::mutate(var_quant = (rank(variance)-1) / (length(variance)-1)) %>%
        dplyr::ungroup() %>%
        dplyr::select(-variance, -bin_med) %>%
        dplyr::mutate(is_hvg = var_quant > 0.75)
    return(result)
}

MN_rowVars <- function(M) {
    if (methods::is(M, "dgCMatrix")) {
        M <- Matrix::t(M)
        result <- Matrix::colMeans(M**2) - Matrix::colMeans(M)**2
        result <- result * nrow(M) / (nrow(M)-1)
    } else {
        result <- matrixStats::rowVars(as.matrix(M))
        names(result) <- rownames(M)
    }
    return(result)
}
