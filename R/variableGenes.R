#' Identify a highly variable gene set
#'
#' Identifies genes with high variance compared to their median expression
#' (top quartile) within each experimentCertain function
#'
#' @param dat SummarizedExperiment object containing gene-by-sample
#' expression matrix.
#' @param exp_labels character vector that denotes the source (Study ID) of
#' each sample.
#' @param i default value 1; non-zero index value of assay containing the matrix data

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

variableGenes <- function(dat, exp_labels, i = 1) {

    dat <- SummarizedExperiment::assay(dat, i = i)
    var_genes1 <- vector("list")

    #check length of exp_labels equal # of samples
    if(length(exp_labels) != length(colnames(dat))){
        stop('experiment_labels length does not match number of samples')
    }

    #check obj contains more than 1 unique study_id
    if(length(unique(exp_labels)) < 2){
        stop('Found only 1 unique exp_labels. Please use data from more than 1 study!')
    }


    experiments <- unique(exp_labels)
    j <- 1
    for(exp in experiments){
        data_subset   <- as.matrix(dat[ , exp_labels == exp])
        genes_list    <- vector("list")
        median_data <- matrixStats::rowMedians(data_subset)
        variance_data <- matrixStats::rowVars(data_subset)
        names(variance_data) <- rownames(data_subset)
        quant_med     <- unique(stats::quantile(median_data,
                                                probs = seq(from = 0,
                                                            to = 1,
                                                            length = 11),
                                                type = 5))
        genes_list    <- vector("list", length = length(quant_med))

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
        var_genes1[[j]] <- unlist(genes_list[1:temp-1])
        j <- j+1
    }
    var_genes <- Reduce(intersect, var_genes1)
    return(var_genes)
}
