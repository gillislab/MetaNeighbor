% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/split_data.R
\name{splitTrainClusters}
\alias{splitTrainClusters}
\title{Split train clusters according to AUROC similarity to test clusters.}
\usage{
splitTrainClusters(mn_scores, k)
}
\arguments{
\item{mn_scores}{An AUROC matrix as generated by MetaNeighborUS, usually with
the "trained_model" option.}

\item{k}{The number of desired cluster sets.}
}
\value{
A list of cluster sets, each cluster set is a character vector
containg cluster labels.
}
\description{
This function computes hierarchical clustering to group similar train
clusters, using similarity to test clusters as features, then uses a standard
tree cutting algorithm to obtain groups of similar clusters. Note that the
cluster hierarchy corresponds exactly to the column dendrogram shown when
using the plotHeatmapPretrained function.
}
\seealso{
\code{\link{plotHeatmapPretrained}}
}
