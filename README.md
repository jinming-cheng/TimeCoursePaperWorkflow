# Unraveling the timeline of gene expression: a pseudo-temporal trajectory analysis of single-cell RNA sequencing data

**Authors**: Jinming Cheng, Gordon K. Smyth and Yunshun Chen

The workflow uses open-source R software packages and covers all steps of the analysis pipeline, including quality control, doublet prediction, normalization, integration, dimension reduction, cell clustering, trajectory inference, and pseudo-bulk time course analysis.
Sample integration and cell clustering follows the *Seurat* pipeline while the trajectory inference is conducted using the *monocle3* package.
The pseudo-bulk time course analysis uses the quasi-likelihood framework of *edgeR*.


# Package installation

The packages used for this workflow can be installed as follows.

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("edgeR")
BiocManager::install("org.Mm.eg.db")
BiocManager::install("scDblFinder")
BiocManager::install("ComplexHeatmap")

install.packages("Seurat")
install.packages("remotes")
remotes::install_github("satijalab/seurat-wrappers")
install.packages("ggplot2")
install.packages("patchwork")

BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'Matrix.utils'))

install.packages("devtools")
devtools::install_github('cole-trapnell-lab/leidenbase')
devtools::install_github('cole-trapnell-lab/monocle3') 
```
