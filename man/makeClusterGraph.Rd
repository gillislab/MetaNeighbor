% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/graph_visualization.R
\name{makeClusterGraph}
\alias{makeClusterGraph}
\title{Convert AUROC matrix into a graph.}
\usage{
makeClusterGraph(best_hits, low_threshold = 0, high_threshold = 1)
}
\arguments{
\item{best_hits}{Matrix of AUROCs produced by MetaNeighborUS.}

\item{low_threshold}{AUROC threshold value. An edge is drawn between two
clusters only if their similarity exceeds low_threshold.}

\item{high_threshold}{AUROC threshold value. An edge is drawn between two
clusters only if their similarity is lower than high_threshold
(enables focusing on close calls).}
}
\value{
A graph in igraph format, where nodes are clusters and edges are
AUROC similarities.
}
\description{
This representation is a useful alternative for heatmaps for large datasets
and sparse AUROC matrices (MetaNeighborUS with one_vs_best = TRUE)
}
