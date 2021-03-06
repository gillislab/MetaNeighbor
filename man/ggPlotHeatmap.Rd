% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{ggPlotHeatmap}
\alias{ggPlotHeatmap}
\title{Plots symmetric AUROC heatmap, clustering cell types by similarity.}
\usage{
ggPlotHeatmap(aurocs, label_size = 10)
}
\arguments{
\item{aurocs}{A square AUROC matrix as returned by MetaNeighborUS.}

\item{label_size}{Font size of cell type labels along the heatmap
(default is 10).}
}
\value{
A ggplot object.
}
\description{
This function is a ggplot alternative to plotHeatmap (without the cell type
dendrogram).
}
\seealso{
\code{\link{plotHeatmap}}
}
