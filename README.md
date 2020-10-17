MetaNeighbor: a method to rapidly assess cell type identity using both functional and random gene sets
================
MetaNeighbor allows users to quantify cell type replicability across datasets using neighbor voting.

Please refer to its [online](./Documentation.md) or [pdf](https://bioconductor.org/packages/release/bioc/vignettes/MetaNeighbor/inst/doc/MetaNeighbor.pdf) documentation and consider citing [Crow et al (2018) Nature Communications](https://www.nature.com/articles/s41467-018-03282-0) if you find MetaNeighbor useful in your research.

## Quick installation procedure

MetaNeighbor has been tested on Windows 10, MacOS Catalina 10.15 and Linux RHEL7 and is expected to run on reasonably up-to-date R (tested on versions 3.6 and 4.0). The main dependencies are the tidyverse, igraph and SingleCellExperiment libraries (full list can be found in [DESCRIPTION](./DESCRIPTION)), all missing dependencies will be automatically installed by running the following commands.

To install the stable version of MetaNeighbor, we recommend using [Bioconductor](https://www.bioconductor.org/install/):

```{r}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()
BiocManager::install("MetaNeighbor")
```

To install the development version of MetaNeighbor, we recommend installing the Github version:

```{r}
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("gillislab/MetaNeighbor")
```

Installation usually completes in 1 or 2 minutes, but can take up to 20 minutes if you are starting with an empty R distribution.

## MetaNeighbor demos
To run a demo of MetaNeighbor, we recommend running the [vignette](./vignettes/MetaNeighbor.Rmd) used for the [documentation](./Documentation.md) or try one of our [protocols](https://github.com/gillislab/MetaNeighbor-Protocol).