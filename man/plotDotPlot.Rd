% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{plotDotPlot}
\alias{plotDotPlot}
\title{Plot dot plot showing expression of a gene set across cell types.}
\usage{
plotDotPlot(
  dat,
  experiment_labels,
  celltype_labels,
  gene_set,
  i = 1,
  normalize_library_size = TRUE,
  alpha_row = 10,
  average_expressing_only = FALSE
)
}
\arguments{
\item{dat}{A SummarizedExperiment object containing gene-by-sample
expression matrix.}

\item{experiment_labels}{A vector that indicates the source/dataset
of each sample.}

\item{celltype_labels}{A character vector that indicates the cell type of
each sample.}

\item{gene_set}{Gene set of interest provided as a vector of genes.}

\item{i}{Default value 1; non-zero index value of assay containing the matrix
data.}

\item{normalize_library_size}{Whether to apply library size normalization
before computing average expression (set this value to FALSE if data are
already normalized).}

\item{alpha_row}{Parameter controling row ordering: a higher value of
alpha_row gives more weight to extreme AUROC values (close to 1).}

\item{average_expressing_only}{Whether average expression should be computed
based only on expressing cells (Seurat default) or taking into account zeros.}
}
\value{
a ggplot object.
}
\description{
The size of each dot reflects the number of cell that express a gene,
the color reflects the average expression.
Expression of genes is first average and scaled in each dataset
independently. The final value is obtained by averaging across datasets.
}
