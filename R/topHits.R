#' Find reciprocal top hits
#'
#' Identifies reciprocal top hits and high scoring cell type pairs. This
#' function calls the newer topHitsByStudy function. 
#'
#' @param cell_NV matrix of celltype-to-celltype AUROC scores
#' (output from \code{\link{MetaNeighborUS}})
#' @param dat a SummarizedExperiment object containing gene-by-sample
#' expression matrix.
#' @param i default value 1; non-zero index value of assay containing the matrix
#' data
#' @param study_id a vector that lists the Study (dataset) ID for each sample
#' @param cell_type a vector that lists the cell type of each sample
#' @param threshold default value 0.95. Must be between [0,1]
#'
#' @return Function returns a dataframe with cell types that are either reciprocal best 
#' matches, and/or those with AUROC values greater than or equal to threshold 
#' value

#' @examples
#' data(mn_data)
#' var_genes = variableGenes(dat = mn_data, exp_labels = mn_data$study_id)
#' celltype_NV = MetaNeighborUS(var_genes = var_genes, 
#'                              dat = mn_data, 
#'                              study_id = mn_data$study_id,
#'                              cell_type = mn_data$cell_type)
#' top_hits = topHits(cell_NV = celltype_NV,
#'                    dat = mn_data,
#'                    study_id = mn_data$study_id,
#'                    cell_type = mn_data$cell_type,
#'                    threshold = 0.9)
#' top_hits
#'
#' @seealso \code{\link{topHitsByStudy}}
#' @export
#'

topHits <- function(cell_NV, dat, i = 1, study_id, cell_type, threshold=0.95){
  by_study_results <- topHitsByStudy(cell_NV, threshold=threshold, n_digits = 2, collapse_duplicates = TRUE)
  by_study_results <- dplyr::rename(by_study_results, Mean_AUROC = AUROC)
  return (by_study_results)
}

#' Find reciprocal top hits, stratifying results by study.
#'
#' This function looks for reciprocal top hits in each target study
#' separately, allowing for as many reciprocal top hits as target studies.
#' This is the recommended function for extracting top hits.
#'
#' @param auroc matrix of celltype-to-celltype AUROC scores
#' (output from \code{\link{MetaNeighborUS}})
#' @param threshold AUROC threshold, must be between [0,1]. Default is 0.9.
#'  Only top hits above this threshold are included in the result table.
#' @param n_digits Number of digits for AUROC values in the result table. Set
#'  to "Inf" to skip rounding.
#' @param collapse_duplicates Collapse identical pairs of cell types (by
#'  default), effectively averaging AUROCs when reference and target roles are
#'  reversed. Setting this option to FALSE makes it easier to filter results
#'  by study or cell type.
#'  If collapse_duplicates is set to FALSE, "Celltype_1" is the
#'  reference cell type and "Celltype_2" is the target cell type (relevant
#'  if MetaNeighborUS was run with symmetric_output = FALSE). 
#'
#' @return Function returns a dataframe with cell types that are either reciprocal best 
#' matches, and/or those with AUROC values greater than or equal to threshold 
#' value
#'
#' @examples
#' data(mn_data)
#' var_genes = variableGenes(dat = mn_data, exp_labels = mn_data$study_id)
#' aurocs = MetaNeighborUS(var_genes = var_genes, 
#'                         dat = mn_data, 
#'                         study_id = mn_data$study_id,
#'                         cell_type = mn_data$cell_type)
#' top_hits = topHitsByStudy(aurocs)
#' top_hits
#'
#' @seealso \code{\link{topHits}}
#' @export
#'
topHitsByStudy = function(auroc, threshold = 0.9, n_digits = 2, collapse_duplicates = TRUE) {
    `%>%` <- dplyr::`%>%`
    #could be sped up by finding same study AUROC's when in matrix form and masking them to NA (as in the old topHits function)
    result <- tibble::as_tibble(auroc, rownames = "ref_cell_type") %>%
        tidyr::pivot_longer(cols = -ref_cell_type,
                            names_to = "target_cell_type",
                            values_to = "auroc") %>%
        dplyr::mutate(ref_study = getStudyId(ref_cell_type),
                      target_study = getStudyId(target_cell_type)) %>%
        dplyr::filter(ref_study != target_study) %>%
        dplyr::group_by(ref_cell_type, target_study) %>%
        dplyr::filter(auroc == max(auroc, na.rm = TRUE)) %>%
        dplyr::ungroup() %>%
        dplyr::select(-ref_study, -target_study) %>%
        dplyr::mutate(is_reciprocal = is_reciprocal_top_hit(.)) %>%
        dplyr::filter(auroc >= threshold) 
    
    if (collapse_duplicates) {
        result <- result %>%
            dplyr::group_by(ref_cell_type, target_cell_type) %>%
            dplyr::mutate(pair_id = paste(sort(c(ref_cell_type, target_cell_type)),
                                          collapse = "")) %>%
            dplyr::group_by(pair_id) %>%
            dplyr::summarize(ref_cell_type = dplyr::first(ref_cell_type),
                             target_cell_type = dplyr::first(target_cell_type),
                             auroc = mean(auroc),
                             is_reciprocal = dplyr::first(is_reciprocal)) %>%
            dplyr::ungroup() %>%
            dplyr::select(-pair_id)
    }
    
    # final formatting
    result <- result %>%
        dplyr::arrange(desc(auroc)) %>%
        dplyr::mutate(auroc = round(auroc, n_digits)) %>%    
        dplyr::mutate(Match_type = ifelse(is_reciprocal,
                                          "Reciprocal_top_hit",
                                          paste0("Above_", threshold))) %>%
        dplyr::select("Study_ID|Celltype_1" = ref_cell_type,
                      "Study_ID|Celltype_2" = target_cell_type,
                      "AUROC" = auroc,
                      Match_type)
    return(result)
}

is_reciprocal_top_hit = function(best_hits) {
    `%>%` <- dplyr::`%>%`
    best_hits <- best_hits %>%
        dplyr::select(-auroc)
    reverse_hits <- best_hits %>%
        dplyr::select(target_cell_type = ref_cell_type,
                      reciprocal_cell_type = target_cell_type)
    reciprocal_best_hits <- dplyr::inner_join(best_hits, reverse_hits,
                                              by = "target_cell_type") %>%
        dplyr::filter(ref_cell_type == reciprocal_cell_type) %>%
        dplyr::select(-reciprocal_cell_type) %>%
        tibble::add_column(is_reciprocal = TRUE)
    result <- dplyr::left_join(best_hits, reciprocal_best_hits,
                              by = c("ref_cell_type", "target_cell_type")) %>%
        tidyr::replace_na(replace = list(is_reciprocal = FALSE)) %>%
        dplyr::pull(is_reciprocal)
    return(result)
}

