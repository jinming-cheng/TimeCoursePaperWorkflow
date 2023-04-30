# Unraveling the timeline of gene expression: a pseudo-temporal trajectory analysis of single-cell RNA sequencing data

The workflow uses open-source R software packages and covers all steps of the analysis pipeline, including quality control, doublet prediction, normalization, integration, dimension reduction, cell clustering, trajectory inference, and pseudo-bulk time course analysis.
Sample integration and cell clustering follows the *Seurat* pipeline while the trajectory inference is conducted using the *monocle3* package.
The pseudo-bulk time course analysis uses the quasi-likelihood framework of *edgeR*.
