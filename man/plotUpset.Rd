% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{plotUpset}
\alias{plotUpset}
\title{Plot Upset plot showing how replicability depends on input dataset.}
\usage{
plotUpset(metaclusters, min_recurrence = 2, outlier_name = "outliers")
}
\arguments{
\item{metaclusters}{Metaclusters extracted from MetaNeighborUS analysis.}

\item{min_recurrence}{Only show replicability structure for metaclusters
that are replicable across at least min_recurrence datasets.}

\item{outlier_name}{In metaclusters, name assigned to outliers (clusters that
did not match with any other cluster)}
}
\description{
Plot Upset plot showing how replicability depends on input dataset.
}
\examples{
data(mn_data)
var_genes = variableGenes(dat = mn_data, exp_labels = mn_data$study_id)
celltype_NV = MetaNeighborUS(var_genes = var_genes,
                             dat = mn_data,
                             study_id = mn_data$study_id,
                             cell_type = mn_data$cell_type,
                             fast_version = TRUE, one_vs_best = TRUE)
mclusters = extractMetaClusters(celltype_NV)
plotUpset(mclusters)

}
