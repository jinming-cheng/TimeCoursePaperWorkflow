## ----style, include=FALSE---------------------------------------------------------------
knitr::opts_chunk$set(error=FALSE,
                      eval=TRUE, # use FALSE for test
                      prompt=TRUE, 
                      message=FALSE,
                      warning=FALSE,
                      comment=NA)

knitr::opts_chunk$set(fig.width=7, fig.height=7, out.width="100%", fig.align="center", fig.path="Figures/") 
knitr::opts_chunk$set(dpi=300, dev="png", 
                      #dev.args=list(pointsize=15), 
                      cache=FALSE, cache.lazy=FALSE)
options(width=83,digits=3)

options(kableExtra.latex.load_packages = FALSE)

options(timeout=300) # to slove "Timeout of 60 seconds was reached" in download.file 

# replace 'xcolor' with 'color' in the tex file
knitr::knit_hooks$set(document = function(x) {sub('\\usepackage[]{xcolor}', '\\usepackage{color}', x, fixed = TRUE)})


## ----chunkEval, include=FALSE-----------------------------------------------------------
eval_download = TRUE
eval_single_sample_Seurat = TRUE
eval_1st_integration = TRUE
eval_prepare_pseudobulk = TRUE
eval_cds = TRUE
eval_save_rdata = TRUE




## ----data dir---------------------------------------------------------------------------
data_dir <- "data"
if(!dir.exists(data_dir)){dir.create(data_dir, recursive=TRUE)}


## ----rds file dir,include=FALSE---------------------------------------------------------
rds_res_dir <- "results/rds_files"
if(!dir.exists(rds_res_dir)){dir.create(rds_res_dir, recursive=TRUE)}


## ----barcode and matrix files, include=FALSE--------------------------------------------
accessions <-c("GSM4994960","GSM4994962","GSM4994963","GSM2759554","GSM4994967")
stages <- c("E18-ME", "Pre-D5-BL6", "Pre-BL6", "5wk-1", "Adult-BL6")
file_suffixes <- c("barcodes.tsv.gz", "matrix.mtx.gz")




## ----featureFiles,  include=FALSE-------------------------------------------------------
GSE <- c("GSE164017", "GSM2759554")
feature_filenames <- c("GSE164017_features.tsv.gz",
                       "GSM2759554_5wk-1-genes.tsv.gz")




## ----targets----------------------------------------------------------------------------
samples <- c("E18.5-epi", "P5", "Pre-puberty", "Puberty", "Adult")
targets <- data.frame(
    samples=samples,
    stages=stages, 
    accessions=accessions,
    matrix.file = paste0("data/",accessions[1:5],"_",stages[1:5],"-","matrix.mtx.gz"), 
    barcode.file = paste0("data/",accessions[1:5],"_",stages[1:5],"-","barcodes.tsv.gz"), 
    feature.file = paste0("data/",feature_filenames[c(1,1,1,2,1)]))
targets


## ----read and merge data----------------------------------------------------------------
library(edgeR)
dge_all <- list()
for ( i in 1:5 ) {
  y <- read10X(mtx = targets$matrix.file[i], 
      barcodes = targets$barcode.file[i], genes = targets$feature.file[i])
  y$samples$group <- targets$samples[i]
  colnames(y) <- paste0(targets$accessions[i],"-",y$samples$Barcode)
  y$genes$Ensembl_geneid <- rownames(y)
  y$genes <- y$genes[,c("Ensembl_geneid","Symbol")]
  y <- y[!duplicated(y$genes$Symbol),]
  rownames(y) <- y$genes$Symbol
  dge_all[[i]] <- y
}
rm(y)
common.genes <- Reduce(intersect, lapply(dge_all, rownames))
for(i in 1:5) dge_all[[i]] <- dge_all[[i]][common.genes, ]
dge_merged <- do.call("cbind", dge_all)


## ----relevel groups---------------------------------------------------------------------
dge_merged$samples$group <- factor(dge_merged$samples$group, levels=samples)


## ----dim of merged data-----------------------------------------------------------------
dim(dge_merged)
table(dge_merged$samples$group)


## ----quality control--------------------------------------------------------------------
dge_merged$samples$num_exp_gene <- colSums(dge_merged$counts>0)
mito_genes <- rownames(dge_merged)[grep("^mt-",rownames(dge_merged))]
dge_merged$samples$mito_percentage <- 
  colSums(dge_merged$counts[mito_genes,])/
  colSums(dge_merged$counts)*100


## ----set theme and colors for ggplot2---------------------------------------------------
library(ggplot2)
my_theme_ggplot <- theme_classic() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=15,face="bold"),
        plot.title=element_text(size=15,face="bold",hjust=0.5),
        plot.margin=margin(0.5, 0.5, 0.5, 0.5, "cm"))
my_theme_facet <- 
  theme(strip.background=element_rect(colour="white",fill="white"),
        strip.text=element_text(size=15, face="bold",color="black")) 
my_colors_15 <- c("cornflowerblue", "darkorchid1", "firebrick1", "gold",
                  "greenyellow", "mediumspringgreen", "mediumturquoise",
                  "orange1", "pink", "deeppink3", "violet", "magenta",
                  "goldenrod4", "cyan", "gray90")


## ----Figure1, fig.cap="Scatter plots of quality control metrics across all the samples. The plots on the left show library size vs number of genes detected, whereas those on the right show library size vs mitochondria read percentage.",fig.height = 2*5*0.8,fig.width=4*2*0.8,out.width="70%"----
p1 <- ggplot(data = dge_merged$samples,
             aes(x=num_exp_gene, y=lib.size, color = group ) ) + 
      geom_point(size=0.5, show.legend=FALSE) + 
      facet_wrap(group~., ncol=1) + 
      scale_color_manual(values=my_colors_15 ) + 
      labs(x="Number of genes", y="Library size") +
      my_theme_ggplot + my_theme_facet
p2 <- ggplot(data = dge_merged$samples,
             aes(x = mito_percentage, y=lib.size, color = group ) ) + 
      geom_point(size = 0.5, show.legend = FALSE) + 
      facet_wrap(group~., ncol=1) + 
      scale_color_manual(values=my_colors_15) + 
      labs(x="Mito-percentage", y="Library size") +
      my_theme_ggplot + my_theme_facet
patchwork::wrap_plots(p1, p2, ncol=2)












## ----Figure2, fig.cap="UMAP visualization of each individual samples. The UMAP plots, in sequence from the top row to the bottom row, correspond to E18.5-epi, P5, Pre-puberty, Puberty, and Adult, respectively. In each row, cells are coloured by cluster on the left, by Epcam expression level in the middle, and by doublet prediction on the right.", fig.height = 16, fig.width=12, out.width="85%"----
p1 <- lapply(data_seurat,function(x){DimPlot(x, pt.size=0.1, cols=my_colors_15) + 
                ggtitle(x$group[1]) + theme(plot.title=element_text(hjust=0.5))})
p2 <- lapply(data_seurat, FeaturePlot, feature="Epcam", pt.size=0.1)
p3 <- lapply(data_seurat, DimPlot, group.by="db_type", pt.size=0.1, 
             cols=c("gray90", "firebrick1"))
patchwork::wrap_plots(c(p1,p2,p3), nrow=5, byrow=FALSE)


## ----epi clusters-----------------------------------------------------------------------
epi_clusters <- list(
  "E18.5-epi" = 0,
  "P5" = c(1,3),
  "Pre-puberty" = c(0:2, 5),
  "Puberty" = 0:6,
  "Adult" = 0:3
)


## ----epi cell filtering-----------------------------------------------------------------
epi_cells <- list()
for (i in samples) {
  epi_cells[[i]] <- rownames(
    subset(data_seurat[[i]]@meta.data,
      (db_type == "singlet") & (seurat_clusters %in% epi_clusters[[i]])))
}
do.call(c, lapply(epi_cells, length))


## ----merge epi cells--------------------------------------------------------------------
epi_cells <- do.call(c, epi_cells)
dge_merged_epi <- dge_merged[, epi_cells]
seurat_merged <- CreateSeuratObject(counts = dge_merged_epi$counts,
                      meta.data = dge_merged_epi$samples,
                      min.cells = 3, min.features = 0, project = "mammary_epi")


## ----split seurat epi-------------------------------------------------------------------
seurat_epi <- SplitObject(seurat_merged, split.by = "group")
seurat_epi <- lapply(seurat_epi, NormalizeData)
seurat_epi <- lapply(seurat_epi, FindVariableFeatures, nfeatures = 2000)










## ----Figure3, fig.cap="UMAP visualization of the integrated data. Cells are coloured by cluster on the left and by original sample on the right.", fig.height = 5, fig.width=12, out.width="100%"----
seurat_int$group <- factor(seurat_int$group, levels = samples)
p1 <- DimPlot(seurat_int, pt.size = 0.1, cols = my_colors_15)
p2 <- DimPlot(seurat_int, pt.size = 0.1, group.by = "group", 
          shuffle = TRUE, cols = my_colors_15) + labs(title="")
p1 | p2


## ----Figure4,  fig.cap="Feature plots of the integrated data. Genes from the top row to the bottom rows are the markers of basal, LP, ML, cycling cells, and stroma, respectively.", fig.width=4*2, fig.height = 3*5, out.width="60%"----
markers <- c("Krt14", "Acta2", "Csn3","Elf5", "Prlr","Areg", 
             "Hmgb2", "Mki67", "Igfbp7","Fabp4")
DefaultAssay(seurat_int) <- "RNA"
FeaturePlot(seurat_int, order = TRUE, pt.size = 0.1, features = markers, ncol = 2)


## ----table cell number------------------------------------------------------------------
tab_number <- table(seurat_int$group, seurat_int$seurat_clusters)
tab_number


## ----table cell proportion--------------------------------------------------------------
tab_ratio <- round(100*tab_number/rowSums(tab_number), 2)
tab_ratio <- as.data.frame.matrix(tab_ratio)
tab_ratio


## ----Figure5, fig.cap="Bar plot of cell proportion of each cluster in each sample.",fig.width = 6, fig.height = 4, out.width="55%"----
par(mar=c(5, 7, 1, 7), xpd=TRUE)
barplot(t(tab_ratio), col=my_colors_15, xlab="Cell proportion (%)", 
    horiz = TRUE, las=1)
legend("right", inset = c(-0.3,0), legend = 0:5, pch = 15, 
    col=my_colors_15, title="Cluster")










## ----Figure6, fig.cap="UMAP visualization of trajectory and pseudotime computed by monocle3. Cells are coloured by cluster on the left and by pseudotime on the right.", fig.width=10, fig.height = 4, out.width="100%"----
p1 <- plot_cells(cds_obj, color_cells_by="seurat_clusters",
                 group_label_size=4, graph_label_size=3,
                 label_cell_groups=FALSE, label_principal_points=TRUE,
                 label_groups_by_cluster=FALSE) + 
      scale_color_manual(values = my_colors_15)
cds_obj <- order_cells(cds_obj, root_pr_nodes="Y_65")
p2 <- plot_cells(cds_obj, color_cells_by="pseudotime",
                 label_groups_by_cluster=FALSE, label_leaves=FALSE,
                 label_branch_points=FALSE)
p1 | p2


## ----obtain pseudotime------------------------------------------------------------------
seurat_int$pseudotime <- pseudotime(cds_obj)












## ----sample info pseudobulk-------------------------------------------------------------
y$samples[, c("lib.size", "pseudotime", "cell_number")]


## ----cell filtering pseudobulk----------------------------------------------------------
keep_samples <- y$samples$cell_number > 30
y <- y[, keep_samples]


## ----gene filtering pseudobulk----------------------------------------------------------
keep_genes <- filterByExpr(y)
y <- y[keep_genes, , keep.lib.sizes=FALSE]


## ----dim after filtering pseudobulk-----------------------------------------------------
dim(y)


## ----TMM normalization------------------------------------------------------------------
y <- calcNormFactors(y)


## ----Figure7, fig.cap="Multi-dimensional scaling (MDS) plot of the pseudo-bulk samples labelled by pseudotime. Samples are coloured by origianl cell cluster on the left and by developmental stage on the right.", fig.width = 11, fig.height = 5, out.width="100%"----
par(mar = c(5.1, 5.1, 2.1, 2.1), mfrow=c(1,2))
cluster <- y$samples$seurat_clusters
group <- y$samples$group
plotMDS(y, labels = round(y$samples$pseudotime, 2),
    xlim=c(-6,4), ylim=c(-3,3), col=my_colors_15[cluster])
legend("topleft", legend=levels(cluster), col=my_colors_15, pch=16)
plotMDS(y, labels = round(y$samples$pseudotime, 2),
    xlim=c(-6,4), ylim=c(-3,3), col=my_colors_15[group])
legend("topleft", legend=levels(group), col=my_colors_15, pch=16)


## ----construct spline design matrix-----------------------------------------------------
t1 <- y$samples$pseudotime
X <- splines::ns(as.numeric(t1),df = 3)
A <- cbind(1,t1,X)
QR <- qr(A)
r <- QR$rank
R_rank <- QR$qr[1:r,1:r]
Z <- t(backsolve(R_rank,t(A),transpose=TRUE))
Z <- Z[,-1]
design <- model.matrix(~ Z)
design


## ----Figure8, fig.cap="A scatter plot of the biological coefficient of variation (BCV) against the average abundance of each gene. The square-root estimates of the common, trended and gene-wise NB dispersions are shown.", fig.width = 5, fig.height = 4.5, out.width="55%"----
y <- estimateDisp(y, design)
sqrt(y$common.dispersion)
plotBCV(y)


## ----Figure9, fig.cap="A scatter plot of the quarter-root QL dispersion against the average abundance of each gene. Estimates are shown for the raw, trended and squeezed dispersions.", fig.width=5, fig.height = 4.5, out.width="55%"----
fit <- glmQLFit(y, design, robust=TRUE)
plotQLDisp(fit)


## ----F test-----------------------------------------------------------------------------
res <- glmQLFTest(fit, coef=2:4)


## ----DE gene number---------------------------------------------------------------------
summary(decideTests(res))


## ----all DE genes-----------------------------------------------------------------------
topTags(res, n=10L) 


## ----top DE genes-----------------------------------------------------------------------
tab <- topTags(res, n=Inf)$table
tab$trend <- ifelse(tab$logFC.Z1 > 0, "Up", "Down")
tab.up <- tab[tab$trend == "Up", ]
tab.down <- tab[tab$trend == "Down", ]
head(tab.up)
head(tab.down)


## ----Figure10, fig.cap="Scatter plots of expression of top genes along pseudotime. The black dots indicate the observed values, while the red line represents the fitted values calculated along pseudotime.", fig.width = 8.4, fig.height = 5.6, out.width="100%"----
logCPM.obs <- edgeR::cpm(y, log=TRUE, prior.count=fit$prior.count)
logCPM.fit <- edgeR::cpm(fit, log=TRUE)
topGenes <- c(rownames(tab.up)[1:3], rownames(tab.down)[1:3])
par(mfrow=c(2,3))
for(i in 1:6) {
  Symbol <- topGenes[i]
  logCPM.obs.i <- logCPM.obs[Symbol, ]
  logCPM.fit.i <- logCPM.fit[Symbol, ]
  plot(y$samples$pseudotime, logCPM.obs.i, xlab="pseudotime",
       ylab="log-CPM", main=Symbol, pch=16, frame=FALSE)
  lines(y$samples$pseudotime, logCPM.fit.i, col="red", lwd=2)
}


## ----Figure11, fig.cap="Heatmap of top 20 up and top 20 down genes. Rows are genes and columns are pseudo-bulk samples.", fig.width=6, fig.height = 8, out.width="60%"----
topGenes <- c(rownames(tab.up)[1:20], rownames(tab.down)[1:20])
z <- logCPM.obs[topGenes, ]
z <- t(scale(t(z)))
ComplexHeatmap::Heatmap(z, name = "Z score",
    cluster_rows = FALSE,cluster_columns = FALSE)


## ----Test 2nd coefficient---------------------------------------------------------------
res_2 <- glmQLFTest(fit, coef=2)
summary(decideTests(res_2))


## ----goana------------------------------------------------------------------------------
go <- goana(res_2, geneid="ENTREZID", species="Mm")
topGO(go, truncate.term = 30, n=15)


## ----Figure12, fig.cap="Barplot of $-\\log_{10}$ p-values of the top 10 down-regulated GO terms under each GO category.", fig.width=10, fig.height = 7.5, out.width="70%"----
top_go <- rbind.data.frame(topGO(go, ont =c("BP"), sort="Down",n=10),
                           topGO(go, ont =c("CC"), sort="Down",n=10),
                           topGO(go, ont =c("MF"), sort="Down",n=10))
d <- transform(top_go, P_DE = P.Down, neg_log10_P = -log10(P.Down))
d$Term <- factor(d$Term,levels = d$Term)
ggplot(data = d, aes(x = neg_log10_P, y = Term, fill = Ont) ) +
  geom_bar(stat = "identity", show.legend = TRUE) + 
  labs(x="-log10 (P value)", y="", title = "Down") + 
  facet_grid(Ont~.,scales = "free",space = "free") + 
  my_theme_ggplot + my_theme_facet + 
  scale_fill_manual(values = my_colors_15[-2]) +
  theme(strip.text = ggplot2::element_blank())


## ----kegga------------------------------------------------------------------------------
kegg <- kegga(res_2, geneid="ENTREZID", species="Mm")
topKEGG(kegg, truncate.path=40, n=15)


## ----Figure13, fig.cap="Barplot of $-\\log_{10}$ p-values of the top 15 down-regulated KEGG pathways.", fig.width=8, fig.height = 4.5, out.width="55%"----
top_path <- topKEGG(kegg,sort="Down",n=15)
data_for_barplot <- transform(top_path, P_DE=P.Down, neg_log10_P=-log10(P.Down))
data_for_barplot$Pathway <- factor(data_for_barplot$Pathway,
                                   levels=data_for_barplot$Pathway)
ggplot(data=data_for_barplot,aes(x=neg_log10_P, y=Pathway) ) +
    geom_bar(stat="identity", show.legend=FALSE, fill=my_colors_15[1]) + 
    labs(x="-log10 (P value)", y="", title="Down" ) + 
    my_theme_ggplot 


## ----keggLinks--------------------------------------------------------------------------
kegg_links <- getGeneKEGGLinks("mmu") 
p_names <- getKEGGPathwayNames("mmu")
p1 <- p_names[grep("PI3K", p_names$Description), ] 
p1_GeneIDs <- subset(kegg_links, PathwayID == p1$PathwayID)$GeneID
tab_p1 <- tab[tab$ENTREZID %in% p1_GeneIDs, ]
d <- logCPM.obs[tab_p1$Symbol,]
d <- apply(d, 2, mean)
d <- data.frame(avg_logCPM = d, avg_pseudotime = y$samples$pseudotime)
head(d)


## ----Figure14, fig.cap="A smooth curve of PI3K-Akt signaling pathway expression level against pseudotime.", fig.width = 6.5, fig.height = 5, out.width="50%"----
ggplot(data = d,aes(x = avg_pseudotime, y = avg_logCPM) ) +
    geom_smooth(color=my_colors_15[1],se = FALSE) + 
    labs(x="Pseudotime", y="Average log-CPM", 
         title = "PI3K-Akt signaling pathway" ) + 
    my_theme_ggplot 


## ----sessionInfo, eval=TRUE-------------------------------------------------------------
sessionInfo()

