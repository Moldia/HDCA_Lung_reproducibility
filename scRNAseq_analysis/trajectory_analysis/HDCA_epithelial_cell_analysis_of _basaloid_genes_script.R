# load packages
library(plyr)
library(dplyr)
library(Seurat)
library(Biobase)
library(cluster)
library(fastcluster)
library(BBmisc)
library(MAST)
library(clustree)
library(readxl)
#library(loomR)
library(harmony)
library(pheatmap)
library(RColorBrewer)
library(sctransform)
library(SeuratWrappers)
#library(velocyto.R)
library(monocle3)
library(slingshot)
library(BUSpaRse)
library(tidyverse)
library(tidymodels)
library(scales)
library(viridis)
library(Matrix)
library(plot3D)
library(plot3Drgl)
library(plotly)
library(orca)
library(SeuratDisk)
library(future)
library(niceRplots)
library(igraph)
library(scales)
library(rafalib)

DefaultAssay(lung) <- "RNA"

# the list of basaloid genes was derived from the https://www.science.org/doi/10.1126/sciadv.aba1983
# which correspond to the results of Wilcoxon rank-sum test of each epithelial cell-type against the other epithelial
# varieties, using the average gene expression per subject, per cell type. 

aba1983_data_s2 <- read.delim("C:/Users/alex/Desktop/HDCA_lung_project/163k_epi_final_strict/basaloid_markers/aba1983_data_s2.txt", row.names=1)
basaloid <- aba1983_data_s2
basaloid$Dpct <- basaloid$pct.inGroup - basaloid$pct.outGroup
# 6495 genes
basaloid <- basaloid[basaloid$gene %in% all.genes,]
# in our dataset 6025 genes 
basaloid <- basaloid[basaloid$wilcox.fdr <0.0001,] #4781 genes
basaloid <- basaloid[basaloid$Dpct > 0.1,] #2236
basaloid <- basaloid[basaloid$log1p_FC > 0.5,] #1135

cell_types <- levels(as.factor(basaloid$cellType))
cell_types

DefaultAssay(lung) <- "RNA"
a1 <- basaloid[basaloid$cellType=="Aberrant Basaloid",]  #117 genes
levels(lung) <- rev(c("14", "11", "12", "7", "6", "0", "4", "1", "10", "2", "9", "3", "8", "13", "5"))
DotPlot(lung, features = a1$gene, scale = TRUE) + RotatedAxis()

pdf("epi_cells_basaloid_genes_unfiltered_dotplot_20220310.pdf", width=27,height=5,paper="special")
a <- DotPlot(lung, features = a1$gene, scale = TRUE) + RotatedAxis()
plot(a)
rm(a)
dev.off()

# IMPORTANT for that step load the dataset and use the "m" with mean values of all clusters. Then select the genes with mean log expression >0.25, in any of the clusters
m$max <- rowMaxs(as.matrix(m[,c(2:16)]))
high_genes <- m[m$max>0.25, ]
a2 <- a1[a1$gene %in% high_genes$gene,] #96 genes
DotPlot(lung, features = a2$gene, scale = TRUE) + RotatedAxis()

pdf("epi_cells_basaloid_genes_filtered_dotplot_20220310.pdf", width=22,height=5,paper="special")
a <- DotPlot(lung, features = a2$gene, scale = TRUE) + RotatedAxis()
plot(a)
rm(a)
dev.off()


pdf("epi_cells_basaloid_genes_20selected_dotplot_20220310.pdf", width=8,height=5,paper="special")
a <- DotPlot(lung, features = c("KRT17", "F3", "LBH", "FRMD5",  "PTPRE", "TACSTD2", #common
                                "ITGB8", "SCD5", "MACC1", "QSOX1", "PTGS2", "NCEH1", "SPSB1",
                                "ITGAV", "ITGA2", "ITGB6", "LAMC2", "FHL2", "CREB5"
                                ), scale = TRUE) + RotatedAxis()
plot(a)
rm(a)
dev.off()

fig_genes <- c("KRT17", "F3", "LBH", "FRMD5",  "PTPRE", "TACSTD2", #common
               "ITGB8", "SCD5", "MACC1", "QSOX1", "PTGS2", "NCEH1", "SPSB1",
               "ITGAV", "ITGA2", "ITGB6", "LAMC2", "FHL2", "CREB5"
)
write.csv(fig_genes, file="Ext_Data_Fig9_20basaloid_genes.csv")


#create and export a list with the 96 genes that are used for basaloid score calculation
basaloid_hdca_genes <- join(a2, high_genes, by="gene", type="inner")
write.csv(basaloid_hdca_genes, file="epi_163k_basaloid_hdca_96_gene_list.csv")

a3 <- list(a2$gene)
lung <- AddModuleScore(
    object = lung,
    features = a3,
    name =  "Basaloid_features")

pdf("epi_cells_umap2d_basaloid_score_20220310.pdf",width = 8.87, height = 7,paper="special")
a <- FeaturePlot(lung, reduction = "umap2d", features = "Basaloid_features1", pt.size = 1, order = F)+ scale_colour_gradientn(colors = c("orange", "gray", "blue"))
plot(a)
rm(a)
dev.off()

pdf("epi_cells_umap2d_basaloid_score_20220611_blue.pdf",width = 8.87, height = 7,paper="special")
a <- FeaturePlot(lung, reduction = "umap2d", features = "Basaloid_features1", pt.size = 1, order = F)
plot(a)
rm(a)
dev.off()


levels(lung) <- c("14", "11", "12", "7", "6", "0", "4", "1", "10", "2", "9", "3", "8", "13", "5")
n_clusters <- nlevels(lung$seurat_clusters)
colors <- DiscretePalette(n_clusters, palette = "alphabet2")
colors <- colors[c(15, 12, 13, 8, 7, 1, 5, 2, 11, 3, 10, 4, 9, 14, 6)]

pdf("epi_cellst_basaloid_score_vlnplot_20220310.pdf",width = 12, height = 4,paper="special")
a <- VlnPlot(lung, features = "Basaloid_features1", pt.size = 0, cols = colors) + NoLegend()
plot(a)
rm(a)
dev.off()

pdf("epi_cellst_Krt8_basaloid_score_vlnplot_dots_20220310.pdf",width = 12, height = 4,paper="special")
a <- VlnPlot(lung, features = "Basaloid_features1", pt.size = 0.1, cols = colors) + NoLegend()
plot(a)
rm(a)
dev.off()

# select the selective markers of cluster3 and 4 and ask, which of the 117 genes are common
activated_epi <- lung_cluster.markers[lung_cluster.markers$cluster==c(3, 4),]
activated_epi <- activated_epi[activated_epi$p_val_adj<0.001,]
a3 <- unique(activated_epi$gene) 

activated_common <- a1[a1$gene %in% a3,]
DotPlot(lung, features = activated_common$gene, scale = TRUE) + RotatedAxis()
DotPlot(lung, features = activated_common$gene, scale = F) + RotatedAxis()

