#' Find reciprocal top hits
#'
#' Identifies reciprocal top hits and high scoring cell type pairs
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
#' @export
#'

topHits <- function(cell_NV, dat, i = 1, study_id, cell_type, threshold=0.95){
    samples <- colnames(dat)
    study_id <- as.character(study_id)
    cell_type <- as.character(cell_type)
    
    #check obj contains study_id
    if(length(study_id)!=length(samples)){
        stop('study_id length does not match number of samples')
    }
    
    #check obj contains cell_type
    if(length(cell_type)!=length(samples)){
        stop('cell_type length does not match number of samples')
    } 
    
    pheno <- as.data.frame(cbind(study_id,cell_type), stringsAsFactors = FALSE)
    pheno$StudyID_CT <- makeClusterName(pheno$study_id, pheno$cell_type)
    
    type_by_study <- table(pheno[,"StudyID_CT"])
    m <- match(rownames(cell_NV), rownames(type_by_study))
    f_a <- !is.na(m)
    f_b <- m[f_a]
    cell_NV <- cell_NV[f_a,f_a]
    type_by_study <- type_by_study[f_b]

    # remove within-dataset scores
    for(i in unique(pheno$study_id)){
        filt <- grepl(i, row.names(type_by_study != 0))
        cell_NV[filt,filt] <- 0
    }

    # remove self-scores
    diag(cell_NV) <-0
    temp <- vector(length = length(rownames(cell_NV)))
    geneInd <- vector(length = length(rownames(cell_NV)))
    
    # identify top hits
    for(i in seq_len(dim(cell_NV)[1])){
        val <- which.max(cell_NV[i,])
        temp[i] <- val
        geneInd[i] <- names(val)
    }

    temp <- cbind(rownames(cell_NV), temp)
    for(i in seq_len(dim(cell_NV)[1])){
        temp[i,2]=cell_NV[i,as.numeric(temp[i,2])]
    }
    
    rownames(temp) <- geneInd
    
    recip <- temp[duplicated(temp[,2]),]
    filt  <- as.numeric(temp[,2]) >= threshold
    recip <- rbind(recip,temp[filt,])
    recip <- cbind(recip, c(rep("Reciprocal_top_hit",
                                each=dim(recip)[1]-sum(filt)),
                            rep(paste("Above",threshold,sep="_"),
                                each=sum(filt))))
    recip <- recip[!duplicated(recip[,2]),]

    #tidy results
    recip2  <- cbind(rownames(recip),recip[,1:3])
    colnames(recip2) <- c("Study_ID|Celltype_1","Study_ID|Celltype_2","Mean_AUROC","Match_type")
    rownames(recip2) <- NULL

    recip     <- recip2[order(recip2[,3],decreasing=TRUE),]
    recip2    <- as.data.frame(recip)
    recip2[,3]<- round(as.numeric(as.character(recip2[,3])),2)
    return(recip2)
}
