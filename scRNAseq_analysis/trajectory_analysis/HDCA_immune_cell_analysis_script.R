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

# increase the maximum of allowed RAM
options(future.globals.maxSize = 30000 * 1024^2)

load("immune_163k_org_data_with.velo.RData")
#=============================================

# create an object list based on the original datasets
lung <- immune
rm(immune)

DefaultAssay(lung) <-"RNA"

lung <- SetIdent(lung, value = lung@meta.data[["donor"]]) 
all.genes <- rownames(lung)

# human orthologues of stress induced genes by tissue enzymatic dissociation at 37oC (https://doi.org/10.1186/s13059-020-02048-6)
stress_den <- c("KLF6", "GEM", "TNFAIP3", "EGR1", "AC005192.1", "IFRD1", "RNU12", "GADD45G", "SNORD13", "RNU6-50P", "ARC", "MIR28", "FOSB",
                "DUSP1", "KLF4", "MAFF", "GDF15", "ZFP36", "FOS", "RNU6-214P", "NFKBID", "NFKBIZ", "JUN", "MIR23A", "IER2", "RASD1", "JUNB",
                "HSPA1B", "HSPA1A", "HSPA1B", "HSPA1A", "PLK2", "TNF", "DUSP5", "GADD45B", "RRAD", "NR4A2", "PPP1R15A", "NR4A1", "MIR10B", "BTG2",
                "ELF3", "SNORD92", "RGS1", "CCL18", "CCL3", "CCL3L3", "CCL4L2", "CCL4", "HCAR2", "HCAR3", "DUSP2", "ATF3", "CSRNP1", "CLDN4", "SOCS3")

stress_den <- stress_den[stress_den %in% all.genes]

stress_features <- list(stress_den)
lung <- AddModuleScore(
  object = lung,
  features = stress_features,
  name = 'stress_Features'
)

lung.list <- SplitObject(lung, split.by = "donor")

rm(lung)

# run the SCTransform regressing out stress induced genes
for (i in 1:length(lung.list)) {
  lung.list[[i]] <- SCTransform(lung.list[[i]], vars.to.regress = c("stress_Features1"), variable.features.n = 5000, verbose = TRUE)
}
#=============================================
# Select features for downstream integration, and run PrepSCTIntegration, which ensures that all necessary Pearson residuals have been calculated.
lung.features <- SelectIntegrationFeatures(object.list = lung.list, nfeatures = 5000)
lung.list <- PrepSCTIntegration(object.list = lung.list, anchor.features = lung.features,
                                verbose = TRUE)

# Next, identify anchors and integrate the datasets. The small number of immune cells from a specific donor made necessary the change of default parameters
lung.anchors <- FindIntegrationAnchors(object.list = lung.list, normalization.method = "SCT", dims = 1:26, k.anchor = 26, k.filter = 26,
                                       k.score= 26, anchor.features = lung.features, verbose = TRUE)


rm(lung.list)
lung <- IntegrateData(anchorset = lung.anchors, normalization.method = "SCT", 
                      verbose = TRUE, k.weight = 26)
rm(lung.anchors, lung.features)
save.image("immune_batch_integrated.RData")

load("immune_batch_integrated.RData")

#Perform linear dimensional reduction
# Next we perform PCA on the scaled data. By default, only the previously determined variable features are used as input, but can be defined using 
# features argument if you wish to choose a different subset.
lung <- RunPCA(lung, npcs = 100)

# set the number of principle components to use for dimension reduction
pcs=50
dim_n = 1:pcs

lung <- FindNeighbors(lung, dims = dim_n)

lung <- FindClusters(lung, resolution = 0.8)


#Run non-linear dimensional reduction (UMAP)
# 1.
lung <- RunUMAP(lung, dims = dim_n)

# Plot on Umap 
DimPlot(lung, reduction = "umap", label = TRUE, pt.size = 1)

#======================================================================
# change the active assay from "integrated" to "RNA"
DefaultAssay(lung) <- "RNA"

# to be able to do it with the correct values, we have to normalize and scale the RNA counts
lung <- NormalizeData(lung, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(lung)
lung <- ScaleData(lung, features = all.genes)

FeaturePlot(lung, reduction = 'umap', features = "COL1A2", pt.size = 2, order = TRUE)
FeaturePlot(lung, reduction = 'umap', features = "PTPRC", pt.size = 2, order = TRUE)
FeaturePlot(lung, reduction = 'umap', features = "HBA1", pt.size = 2, order = TRUE)

VlnPlot(lung, features = "PTPRC", pt.size = 0)
#====================================================================

# remove the artifact cluster with the doublets
lung <- subset(lung, subset = HBA1 < 4)
lung <- subset(lung, idents= c(0:9, 11:15) )
DefaultAssay(lung) <-"integrated"

lung #4113 cells
save(lung, file = "immune_batch_integrated_filtered_org.RData")

#====================================================================
#====================================================================
load("immune_batch_integrated_filtered_org.RData")

# add numeric age to the metadata
ages <- data.frame(lung@meta.data[["age"]])
row.names(ages) <- colnames(lung)
names(ages)[] <- "age_n"
ages$age_n<- as.numeric(str_replace(ages$age_n, "w", ""))
lung <- AddMetaData(lung, ages$age_n, col.name = "age_n")

# after filtering of the cells we run again pca 
DefaultAssay(lung) <- "integrated"
lung <- RunPCA(lung, npcs = 50)

# set the number of principle components to use for dimension reduction
pcs=50
dim_n = 1:pcs
lung <- FindNeighbors(lung, dims = dim_n, reduction = "pca", assay = "integrated",k.param = 25, n.trees = 500)

# Run UMAP on a graph
lung <- FindClusters(lung, resolution = 1.64)
# run tSNE
lung <- RunTSNE(lung, reduction = "pca", dims = dim_n)
# Plot on tSNE
DimPlot(lung, reduction = "tsne", label = TRUE, pt.size = 1) + NoLegend()

# run UMAP
lung <- RunUMAP(lung, 
                assay = "integrated",
                n.neighbors = 25L, 
                min.dist = 0.5, 
                dims = dim_n, 
                spread = 1,
                metric = "cosine",
                repulsion.strength = 0.5,
                negative.sample.rate = 10,
                n.epochs = NULL,
                seed.use = 42L,
                learning.rate = 15,
                n.components = 2L,
                reduction.name = "umap2d",
                reduction.key='umap2d_')


DimPlot(lung, reduction = "umap2d", label = TRUE, pt.size = 1)

#================================================================================
# set the resolution and run 3d UMAP

lung <- FindClusters(lung, resolution = 1.6)
DimPlot(lung, reduction = "umap2d", label = TRUE, pt.size = 1)

lung <- RunUMAP(lung, 
                assay = "integrated",
                n.neighbors = 25L, 
                min.dist = 0.5, 
                dims = dim_n, 
                spread = 0.5,
                metric = "cosine",
                repulsion.strength = 15,
                negative.sample.rate = 0.45,
                seed.use = 42,
                learning.rate = 15,
                n.epochs = NULL,
                n.components = 3L)

# plot of 3D UMAP-plots
umap_1 <- lung[["umap"]]@cell.embeddings[,1]
umap_2 <- lung[["umap"]]@cell.embeddings[,2]
umap_3 <- lung[["umap"]]@cell.embeddings[,3]

# Visualize what headings are called so that you can extract them to form a dataframe
Embeddings(object = lung, reduction = "umap")

# obtain cluster colors, similar to 2D default Seurat plot
n_clusters <- nlevels(lung$seurat_clusters)
#colors <- hue_pal()(n_clusters)
colors <- DiscretePalette(n_clusters, palette = "alphabet2")

# Prepare a dataframe for cell plotting
plot.data <- FetchData(object = lung, vars = c("UMAP_1", "UMAP_2", "UMAP_3", "seurat_clusters"))

# Make a column of row name identities (these will be your cell/barcode names)
plot.data$label <- paste(rownames(plot.data))

dev.off()

# Plot your data
plot_ly(data = plot.data, 
        x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, 
        color = ~seurat_clusters, 
        colors = colors,
        type = "scatter3d", 
        mode = "markers",
        marker = list(size = 4, width=4), # controls size of points
        text=~label, #This is that extra column we made earlier for which we will use for cell ID
        hoverinfo="text") #When you visualize your plotly object, hovering your mouse pointer over a point shows cell names


# plot the 3D UMAP using plot3D package, theta and phi change the orientation of the graph 
scatter3D(x = plot.data[,1], y = plot.data[,2], z = plot.data[,3], 
          xlab= "UMAP-1", ylab= "UMAP-2", zlab= "UMAP-3",
          colvar= as.integer(plot.data[,4]),
          col = colors,
          bty = "b2", 
          pch = 20, cex = 1, ticktype = "detailed",
          theta = 185, phi = 5,
          clab= "Clusters"
)



pdf("immune_163k_umap3d_clusters.pdf", width = 11.69, height = 8.27, paper = "special")
pdf <-scatter3D(x = plot.data[,1], y = plot.data[,2], z = plot.data[,3], 
                xlab= "UMAP-1", ylab= "UMAP-2", zlab= "UMAP-3",
                colvar= as.integer(plot.data[,4]),
                col = colors,
                bty = "b2", 
                pch = 20, cex = 1, ticktype = "detailed",
                theta = 185, phi = 5,
                clab= "Clusters"
)
plot(pdf)
dev.off()
rm(pdf)
#-----------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------
# reorder the ages
lung@meta.data[["age"]] <- factor(lung@meta.data[["age"]], levels= c("5w", "5.5w", "6w", "7w", "8w", "8.5w", "10w", "11.5w", "12w", "13w", "14w"))

# plot the cells according to metadata parameters
# set the metadata to plot
metadata_to_plot <- "age"

# obtain colors, similar to 2D default Seurat plot
n_vars <- nlevels(as.factor(lung@meta.data[[metadata_to_plot]]))
#colors <- hue_pal()(n_vars)

colors1 <- DiscretePalette(n_vars, palette = "alphabet2")


# Prepare a dataframe for cell plotting
plot.data_meta <- FetchData(object = lung, vars = c("UMAP_1", "UMAP_2", "UMAP_3", metadata_to_plot))

# Make a column of row name identities (these will be your cell/barcode names)
plot.data_meta$label <- paste(rownames(plot.data_meta))

# Plot your data, in this example my Seurat object had 21 clusters (0-20)
plot_ly(data = plot.data_meta, 
        x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, 
        color = plot.data_meta[,4], 
        colors = colors1,
        type = "scatter3d", 
        mode = "markers",
        marker = list(size = 4, width=4), # controls size of points
        text=~label, #This is that extra column we made earlier for which we will use for cell ID
        hoverinfo="text") #When you visualize your plotly object, hovering your mouse pointer over a point shows cell names


# plot the 3D UMAP using plot3D package, theta and phi change the orinetation of the graph 
scatter3D(x = plot.data_meta[,1], y = plot.data_meta[,2], z = plot.data_meta[,3], 
          xlab= "UMAP-1", ylab= "UMAP-2", zlab= "UMAP-3",
          colvar= as.integer(as.factor(plot.data_meta[,4])),
          col = colors1,
          bty = "b2", 
          pch = 20, cex = 1, ticktype = "detailed",
          theta = 185, phi = 5,
          clab= metadata_to_plot
)

pdf(paste0("immune_163k_umap3d_",metadata_to_plot,".pdf"), width = 11.69, height = 8.27, paper = "special")
pdf <-scatter3D(x = plot.data_meta[,1], y = plot.data_meta[,2], z = plot.data_meta[,3], 
                xlab= "UMAP-1", ylab= "UMAP-2", zlab= "UMAP-3",
                colvar= as.integer(as.factor(plot.data_meta[,4])),
                col = colors1,
                bty = "b2", 
                pch = 20, cex = 1, ticktype = "detailed",
                theta = 185, phi = 5,
                clab= metadata_to_plot
)
plot(pdf)
dev.off()
rm(pdf)

#=======================================================================================================
DefaultAssay(lung) <- "integrated"

# export Umap plots with the identified clusters, the age and sample name
pdf("immune_cells_umap2d_clusters.pdf",width=11.69,height=8.27,paper="special")
immune_cells_umap2d_clusters <- DimPlot(lung, reduction = "umap2d", cols =  "alphabet2", label = TRUE, pt.size = 3)
plot(immune_cells_umap2d_clusters)
rm(immune_cells_umap2d_clusters)
dev.off()

pdf("immune_cells_umap2d_age.pdf",width=11.69,height=8.27,paper="special")
immune_cells_umap2d_age <- DimPlot(lung, reduction = "umap2d", cols =  "alphabet2", group.by = "age", pt.size = 3)
plot(immune_cells_umap2d_age)
rm(immune_cells_umap2d_age)
dev.off()

pdf("immune_cells_umap2d_SampleName.pdf",width=11.69,height=8.27,paper="special")
immune_cells_umap2d_SampleName <- DimPlot(lung, reduction = "umap2d", cols =  "alphabet2", group.by = "SampleName", pt.size = 3)
plot(immune_cells_umap2d_SampleName)
rm(immune_cells_umap2d_SampleName)
dev.off()

pdf("immune_cells_umap2d_donor.pdf",width=11.69,height=8.27,paper="special")
immune_cells_umap2d_donor <- DimPlot(lung, reduction = "umap2d", cols =  "alphabet2", group.by = "donor", pt.size = 3)
plot(immune_cells_umap2d_donor)
rm(immune_cells_umap2d_donor)
dev.off()

pdf("immune_cells_umap2d_chemistry.pdf",width=11.69,height=8.27,paper="special")
immune_cells_umap2d_chemistry <- DimPlot(lung, reduction = "umap2d", group.by = "chemistry", pt.size = 3)
plot(immune_cells_umap2d_chemistry)
rm(immune_cells_umap2d_chemistry)
dev.off()

pdf("immune_cells_umap2d_features.pdf",width=11.69,height=8.27,paper="special")
immune_cells_umap2d_features <- FeaturePlot(lung, reduction = "umap2d",  features = "nFeature_RNA", pt.size = 3)
plot(immune_cells_umap2d_features)
rm(immune_cells_umap2d_features)
dev.off()

pdf("immune_cells_umap2d_percent_mt.pdf",width=11.69,height=8.27,paper="special")
immune_cells_umap2d_percent_mt <- FeaturePlot(lung, reduction = "umap2d",  features = "MT_ratio", pt.size = 3)
plot(immune_cells_umap2d_percent_mt)
rm(immune_cells_umap2d_percent_mt)
dev.off()

pdf("immune_cells_tsne_clusters.pdf",width=11.69,height=8.27,paper="special")
s10X140_cells_tsne_clusters <- DimPlot(lung, reduction = "tsne", cols =  "alphabet2", label = TRUE, pt.size = 3)
plot(s10X140_cells_tsne_clusters)
rm(s10X140_cells_tsne_clusters)
dev.off()

pdf("immune_cells_tsne_age.pdf",width=11.69,height=8.27,paper="special")
immune_cells_tsne_age <- DimPlot(lung, reduction = "tsne", cols =  "alphabet2", group.by = "age", pt.size = 3)
plot(immune_cells_tsne_age)
rm(immune_cells_tsne_age)
dev.off()

pdf("immune_cells_tsne_SampleName.pdf",width=11.69,height=8.27,paper="special")
immune_cells_tsne_SampleName <- DimPlot(lung, reduction = "tsne", cols =  "alphabet2", group.by = "SampleName", pt.size = 3)
plot(immune_cells_tsne_SampleName)
rm(immune_cells_tsne_SampleName)
dev.off()

pdf("immune_cells_tsne_donor.pdf",width=11.69,height=8.27,paper="special")
immune_cells_tsne_donor <- DimPlot(lung, reduction = "tsne", cols =  "alphabet2", group.by = "donor", pt.size = 3)
plot(immune_cells_tsne_donor)
rm(immune_cells_tsne_donor)
dev.off()

pdf("immune_cells_tsne_chemistry.pdf",width=11.69,height=8.27,paper="special")
immune_cells_tsne_chemistry <- DimPlot(lung, reduction = "tsne", group.by = "chemistry", pt.size = 3)
plot(immune_cells_tsne_chemistry)
rm(immune_cells_tsne_chemistry)
dev.off()

pdf("immune_cells_tsne_percent_mt.pdf",width=11.69,height=8.27,paper="special")
immune_cells_tsne_percent_mt <- FeaturePlot(lung, reduction = "tsne", features = "MT_ratio", pt.size = 3)
plot(immune_cells_tsne_percent_mt)
rm(immune_cells_tsne_percent_mt)
dev.off()

pdf("immune_cells_tsne_features.pdf",width=11.69,height=8.27,paper="special")
immune_cells_tsne_features <- FeaturePlot(lung, reduction = "tsne", features = "nFeature_RNA", pt.size = 3)
plot(immune_cells_tsne_features)
rm(immune_cells_tsne_features)
dev.off()

#=============================================================================================
# create *.tiff files with individual clusters being colored with red

for (i in 0:(nlevels(lung@active.ident)-1)) {
  cluster <- WhichCells(lung, ident = i)
  name = sprintf("immune_umap2d_cluster%s.tiff", i)
  tiff(name, width = 707, height = 500, compression = "lzw")
  image <- DimPlot(lung, reduction = "umap2d",  cells.highlight = cluster, sizes.highlight= 3, pt.size = 3) + NoLegend()
  plot(image)
  rm(image)
  dev.off()
}

for (i in 0:(nlevels(lung@active.ident)-1)) {
  cluster <- WhichCells(lung, ident = i)
  name = sprintf("immune_tsne_cluster%s.tiff", i)
  tiff(name, width = 707, height = 500, compression = "lzw")
  image <- DimPlot(lung, reduction = "tsne", cells.highlight = cluster, sizes.highlight= 3, pt.size = 3) + NoLegend()
  plot(image)
  rm(image)
  dev.off()
}

for (i in 0:(nlevels(lung@active.ident)-1)){
  # Prepare a dataframe for cell plotting
  plot.data_indiv <- FetchData(object = lung, vars = c("UMAP_1", "UMAP_2", "UMAP_3", "seurat_clusters"))
  
  plot.data_indiv$Cluster <- ifelse(plot.data_indiv$seurat_clusters==i, 1, 0)
  
  # Make a column of row name identities (these will be your cell/barcode names)
  plot.data_indiv$label <- paste(rownames(plot.data_indiv))
  
  n_clusters <- nlevels(as.factor(plot.data_indiv$Cluster))
  colors <- c("gray", "red")
  
  pdf(paste0("umap3d_cluster",i,".pdf"), width = 11.69, height = 8.27, paper = "special")
  pdf <-scatter3D(x = plot.data_indiv[,1], y = plot.data_indiv[,2], z = plot.data_indiv[,3], 
                  xlab= "UMAP-1", ylab= "UMAP-2", zlab= "UMAP-3",
                  colvar= as.integer(plot.data_indiv[,5]),
                  col = colors,
                  bty = "b2", 
                  pch = 20, cex = 1, ticktype = "detailed",
                  theta = 140, phi = 20,
                  clab= paste0("Cluster", i)
  )
  plot(pdf)
  dev.off()
  rm(pdf)
}
#==============================================================================================
#==============================================================================================
# differential expression analysis to identify specific or selective genes for each cluster

# change the active assay from "integrated" to "RNA"
DefaultAssay(lung) <- "RNA"

# to be able to do it with the correct values, we have to normalize and scale the RNA counts
lung <- NormalizeData(lung, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(lung)
lung <- ScaleData(lung, features = all.genes)

# Run differential expression analysis with MAST
min_cluster <- min(table(lung@active.ident))
lung_cluster.markers <- FindAllMarkers(lung, test.use = "MAST", min.diff.pct = 0.1, max.cells.per.ident = min_cluster, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.1)
# calculate the percentage difference of the positive cells in the indicated cluster, in comparison to the others
lung_cluster.markers$Dpct <- (lung_cluster.markers$pct.1 - lung_cluster.markers$pct.2)

c <-subset(lung, idents =0 )
counts <- as.matrix(c[["RNA"]]@data)
m <- data.frame(rowMeans(counts))
gene <- row.names(c)
m <- cbind(gene, m)
names(m)[2] <- paste0("cluster", "0") 

for (i in 1:(nlevels(lung@active.ident)-1)) {
  c1 <-subset(lung, idents =i )
  counts1 <- as.matrix(c1[["RNA"]]@data)
  m1 <- data.frame(rowMeans(counts1))
  names(m1)[1] <- paste0("cluster", (i))
  m <- cbind(m, m1)
}
rm(c, counts, counts1, c1, m1)


lung_cluster.markers <-join(lung_cluster.markers, m, by="gene")
write.csv(lung_cluster.markers, file="immune_MAST_lung_markers.csv")


# filter out insignificant genes
lung_cluster.markers <- lung_cluster.markers[lung_cluster.markers$p_val_adj<0.001, ]

# obtain the top 25 genes of each cluster, based on the average logFC
top25 <- lung_cluster.markers %>% group_by(cluster) %>% top_n(n = 25, wt = avg_log2FC)
write.csv(top25, file="immune_MAST_lung_cluster_top25_markers.csv")


top5 <- lung_cluster.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
write.csv(top5, file="immune_MAST_lung_cluster_top5_markers.csv")

# create and export a tiff file with the heatmap of the top 5 genes, based on average logFC
tiff("top5_cluster_heatmap.tiff",width = 2000, height = 1500, compression = "lzw")
top5_cluster_heatmap <-DoHeatmap(lung, features = top5$gene) + NoLegend()
plot(top5_cluster_heatmap)
rm(top5_cluster_heatmap)
dev.off()

pdf("top5_cluster_heatmap.pdf",width = 11.69, height = 8.27, paper = "special")
top5_cluster_heatmap <-DoHeatmap(lung, features = top5$gene) + NoLegend()
plot(top5_cluster_heatmap)
rm(top5_cluster_heatmap)
dev.off()

#-------------------------------------------------------------------------------------------
# The transcription factor and co-factor lists were downloaded from the http://bioinfo.life.hust.edu.cn/AnimalTFDB/#!/download
# The kinase list "protein_class_kinases.tsv", cd marker list "protein_class_CD.tsv" and "protein_class_secreted.tsv" were downloaded from https://www.proteinatlas.org/


# import the list of transcription factors 
TFs <- read.csv("TFs.csv")

# create a new file by merging the TF-list and the enriched genes in each cluster, to identify the TFs, only
lung_cluster_tfs <- join(lung_cluster.markers, TFs, by="gene", type="inner")
top10_cluster_tfs <- lung_cluster_tfs %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
# export the list of the top 10 TFs based on average expression difference
write.csv(top10_cluster_tfs, file="immune_MAST_top10_cluster_tfs.csv")

# create a heatmap with the top 10 enriched TFs
tiff("top10_cluster_tfs_heatmap.tiff",width = 2500, height = 2500, compression = "lzw")
top10_cluster_tfs_heatmap <-DoHeatmap(lung, features = top10_cluster_tfs$gene) + NoLegend()
plot(top10_cluster_tfs_heatmap)
rm(top10_cluster_tfs_heatmap)
dev.off()


# import the list of transcription cofactors 
coTFs <- read.csv("TF_cofactors.csv")

# create a new file by merging the coTF-list and the enriched genes in each cluster, to identify the coTFs, only
lung_cluster_cotfs <- join(lung_cluster.markers, coTFs, by="gene", type="inner")
top10_cluster_cotfs <- lung_cluster_cotfs %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
# export the list of the top 10 coTFs based on average expression difference
write.csv(top10_cluster_cotfs, file="immune_MAST_top10_cluster_cotfs.csv")

# create a heatmap with the top 10 enriched coTFs
tiff("top10_cluster_cotfs_heatmap.tiff",width = 2500, height = 2500, compression = "lzw")
top10_cluster_cotfs_heatmap <-DoHeatmap(lung, features = top10_cluster_cotfs$gene) + NoLegend()
plot(top10_cluster_cotfs_heatmap)
rm(top10_cluster_cotfs_heatmap)
dev.off()


# import the list of kinases
kinases <- data.frame(read.delim("protein_class_kinases.tsv"))
kinases <- kinases[1:10]

# create a new file by merging the kinase-list and the enriched genes in each cluster, to identify the kinases, only
lung_cluster_kinases <- join(lung_cluster.markers, kinases, by="gene", type="inner")
top10_cluster_kinases <- lung_cluster_kinases %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
# export the list of the top 10 kinases based on average expression difference
write.csv(top10_cluster_kinases, file="immune_MAST_top10_cluster_kinases.csv")

# create a heatmap with the top 10 enriched kinases
tiff("top10_cluster_kinases_heatmap.tiff",width = 2500, height = 2500, compression = "lzw")
top10_cluster_kinases_heatmap <-DoHeatmap(lung, features = top10_cluster_kinases$gene) + NoLegend()
plot(top10_cluster_kinases_heatmap)
rm(top10_cluster_kinases_heatmap)
dev.off()

# import the list of CD markers
cds <- data.frame(read.delim("protein_class_CD.tsv"))
cds <- cds[1:10]

# create a new file by merging the CD-list and the enriched genes in each cluster, to identify the surface markers, only
lung_cluster_cds <- join(lung_cluster.markers, cds, by="gene", type="inner")
top10_cluster_cds <- lung_cluster_cds %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
# export the list of the top 10 CDs based on average expression difference
write.csv(top10_cluster_cds, file="immune_MAST_top10_cluster_cds.csv")

# create a heatmap with the top 10 enriched CDs
tiff("top10_cluster_cds_heatmap.tiff",width = 2500, height = 2500, compression = "lzw")
top10_cluster_cds_heatmap <-DoHeatmap(lung, features = top10_cluster_cds$gene) + NoLegend()
plot(top10_cluster_cds_heatmap)
rm(top10_cluster_cds_heatmap)
dev.off()

#----------------------------------------------------------------------------------------------------------------------
# Run differential expression analysis with ROC
DefaultAssay(lung) <-"RNA"
lung_roc_cluster.markers <- FindAllMarkers(lung, test.use = "roc", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.1)
# calculate the percentage difference of the positive cells in the indicated cluster, in comparison to the others
lung_roc_cluster.markers$Dpct <- (lung_roc_cluster.markers$pct.1 - lung_roc_cluster.markers$pct.2)

lung_roc_cluster.markers <-join(lung_roc_cluster.markers, m, by="gene")
write.csv(lung_roc_cluster.markers, file="immune_163k_roc_markers.csv")

lung_roc_cluster.markers <- lung_roc_cluster.markers[lung_roc_cluster.markers$Dpct >0.1, ]

# obtain the top 25 genes of each cluster, based on the power
top25_roc <- lung_roc_cluster.markers %>% group_by(cluster) %>% top_n(n = 25, wt = power)
write.csv(top25_roc, file="immune_163k_roc_cluster_top25_markers.csv")


top5_roc <- lung_roc_cluster.markers %>% group_by(cluster) %>% top_n(n = 5, wt = power)
write.csv(top5_roc, file="immune_163k_roc_cluster_top5_markers.csv")

# create and export a tiff file with the heatmap of the top 5 genes, based on power
tiff("top5_roc_cluster_heatmap.tiff",width = 2000, height = 1500, compression = "lzw")
top5_roc_cluster_heatmap <-DoHeatmap(lung, features = top5_roc$gene) + NoLegend()
plot(top5_roc_cluster_heatmap)
rm(top5_roc_cluster_heatmap)
dev.off()

pdf("top5_roc_cluster_heatmap.pdf",width = 11.69, height = 8.27, paper = "special")
top5_roc_cluster_heatmap <-DoHeatmap(lung, features = top5_roc$gene) + NoLegend()
plot(top5_roc_cluster_heatmap)
rm(top5_roc_cluster_heatmap)
dev.off()

top10_roc <- lung_roc_cluster.markers %>% group_by(cluster) %>% top_n(n = 10, wt = power)
write.csv(top10_roc, file="immune_163k_roc_cluster_top10_markers.csv")

# create and export a tiff file with the heatmap of the top 5 genes, based on power
tiff("top10_roc_cluster_heatmap.tiff",width = 2000, height = 1500, compression = "lzw")
top10_roc_cluster_heatmap <-DoHeatmap(lung, features = top10_roc$gene) + NoLegend()
plot(top10_roc_cluster_heatmap)
rm(top10_roc_cluster_heatmap)
dev.off()
#-------------------------------------------------------------------------------------------
# The transcription factor and co-factor lists were downloaded from the http://bioinfo.life.hust.edu.cn/AnimalTFDB/#!/download
# The kinase list "protein_class_kinases.tsv" was downloaded from https://www.proteinatlas.org/search/protein_class:kinases
# The cd marker list "protein_class_CD.tsv" was downloaded from https://www.proteinatlas.org/search/protein_class%3ACD+markers

# import the list of transcription factors 
TFs <- read.csv("TFs.csv")

# create a new file by merging the TF-list and the enriched genes in each cluster, to identify the TFs, only
lung_roc_cluster_tfs <- join(lung_roc_cluster.markers, TFs, by="gene", type="inner")
top10_roc_cluster_tfs <- lung_roc_cluster_tfs %>% group_by(cluster) %>% top_n(n = 10, wt = power)
# export the list of the top 10 TFs based on average expression difference
write.csv(top10_roc_cluster_tfs, file="immune_163k_roc_top10_cluster_tfs.csv")

# create a heatmap with the top 10 enriched TFs
tiff("top10_roc_cluster_tfs_heatmap.tiff",width = 2500, height = 2500, compression = "lzw")
top10_roc_cluster_tfs_heatmap <-DoHeatmap(lung, features = top10_roc_cluster_tfs$gene) + NoLegend()
plot(top10_roc_cluster_tfs_heatmap)
rm(top10_roc_cluster_tfs_heatmap)
dev.off()


# import the list of transcription cofactors 
coTFs <- read.csv("TF_cofactors.csv")

# create a new file by merging the coTF-list and the enriched genes in each cluster, to identify the coTFs, only
lung_roc_cluster_cotfs <- join(lung_roc_cluster.markers, coTFs, by="gene", type="inner")
top10_roc_cluster_cotfs <- lung_roc_cluster_cotfs %>% group_by(cluster) %>% top_n(n = 10, wt = power)
# export the list of the top 10 coTFs based on average expression difference
write.csv(top10_roc_cluster_cotfs, file="immune_163k_roc_top10_cluster_cotfs.csv")

# create a heatmap with the top 10 enriched coTFs
tiff("top10_roc_cluster_cotfs_heatmap.tiff",width = 2500, height = 2500, compression = "lzw")
top10_roc_cluster_cotfs_heatmap <-DoHeatmap(lung, features = top10_roc_cluster_cotfs$gene) + NoLegend()
plot(top10_roc_cluster_cotfs_heatmap)
rm(top10_roc_cluster_cotfs_heatmap)
dev.off()


# import the list of kinases
kinases <- data.frame(read.delim("protein_class_kinases.tsv"))
kinases <- kinases[1:10]

# create a new file by merging the kinase-list and the enriched genes in each cluster, to identify the kinases, only
lung_roc_cluster_kinases <- join(lung_roc_cluster.markers, kinases, by="gene", type="inner")
top10_roc_cluster_kinases <- lung_roc_cluster_kinases %>% group_by(cluster) %>% top_n(n = 10, wt = power)
# export the list of the top 10 kinases based on average expression difference
write.csv(top10_roc_cluster_kinases, file="immune_163k_roc_top10_cluster_kinases.csv")

# create a heatmap with the top 10 enriched kinases
tiff("top10_roc_cluster_kinases_heatmap.tiff",width = 2500, height = 2500, compression = "lzw")
top10_roc_cluster_kinases_heatmap <-DoHeatmap(lung, features = top10_roc_cluster_kinases$gene) + NoLegend()
plot(top10_roc_cluster_kinases_heatmap)
rm(top10_roc_cluster_kinases_heatmap)
dev.off()

# import the list of CD markers
cds <- data.frame(read.delim("protein_class_CD.tsv"))
cds <- cds[1:10]

# create a new file by merging the CD-list and the enriched genes in each cluster, to identify the surface markers, only
lung_roc_cluster_cds <- join(lung_roc_cluster.markers, cds, by="gene", type="inner")
top10_roc_cluster_cds <- lung_roc_cluster_cds %>% group_by(cluster) %>% top_n(n = 10, wt = power)
# export the list of the top 10 CDs based on average expression difference
write.csv(top10_roc_cluster_cds, file="top10_roc_cluster_cds.csv")

# create a heatmap with the top 10 enriched CDs
tiff("top10_roc_cluster_cds_heatmap.tiff",width = 2500, height = 2500, compression = "lzw")
top10_roc_cluster_cds_heatmap <-DoHeatmap(lung, features = top10_roc_cluster_cds$gene) + NoLegend()
plot(top10_roc_cluster_cds_heatmap)
rm(top10_roc_cluster_cds_heatmap)
dev.off()

#============================================
# create and export a table with the cell names and their corresponding cluster 
immune_sct_cells_clusters <- data.frame(lung@active.ident)
write.csv(immune_sct_cells_clusters, file="immune_cells_clusters.csv")

Age <- data.frame(lung@meta.data[["age"]])
experiment <- data.frame(lung@meta.data[["orig.ident"]])
cell_name <- row.names(immune_sct_cells_clusters)
cluster_list <- cbind(cell_name, immune_sct_cells_clusters)
cluster_list <- cbind(cluster_list, Age)
cluster_list <- cbind(cluster_list, experiment)
cluster_list <- setNames(cluster_list, c("cell", "cluster", "age", "orig_ident"))
write.csv(cluster_list, file="immune_cells_metadata.csv")

# summarize the abundance of cells from specific time points in clusters 
cluster_list$cluster <- as.factor(cluster_list$cluster)
cluster_list$age <- as.factor(cluster_list$age)
abundance_matrix <- matrix(nrow = nlevels(cluster_list$cluster), ncol = nlevels(cluster_list$age))
rownames(abundance_matrix) <- levels(cluster_list$cluster)
colnames(abundance_matrix) <- levels(cluster_list$age)
for (i in 0:nlevels(cluster_list$cluster)-1) 
  for (j in 1:nlevels(cluster_list$age)) {
    list1 <-cluster_list[cluster_list$cluster == i,]
    abundance_matrix[(i+1),j] <- sum(list1$age==levels(list1$age)[j])
  }
write.csv(abundance_matrix, file="immune_cells_age_abundance_in_clusters.csv")


# calculate the % of cells of each time-point in the clusters
cluster_size <- rowSums(abundance_matrix)
abundance_matrix_percentages <- (abundance_matrix*100)/cluster_size
write.csv(abundance_matrix_percentages, file="immune_cells_age_abundance_percentages_in_clusters.csv")

save.image("immune_163K_after_MAST_and_roc_analysis_data.RData")

#=====================================================================================================
# differential expression analysis of specific clusters
#===============================================================================================================
# run differential expression analyses of distinct clusters
DefaultAssay(lung) <-"RNA"
mf_cl2_vs0_16 <- FindMarkers(lung, ident.1 = "2", ident.2 = c("0", "16"), test.use = "MAST", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.1)
mf_cl2_vs0_16$gene <- row.names(mf_cl2_vs0_16)
mf_cl2_vs0_16 <- join(mf_cl2_vs0_16, m, by="gene")
write.csv(mf_cl2_vs0_16, file= "mf_cl2_vs0_16_dif_expression_analysis.csv")

mf_cl0_vs2_16 <- FindMarkers(lung, ident.1 = "0", ident.2 = c("2", "16"), test.use = "MAST", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.1)
mf_cl0_vs2_16$gene <- row.names(mf_cl0_vs2_16)
mf_cl0_vs2_16 <- join(mf_cl0_vs2_16, m, by="gene")
write.csv(mf_cl0_vs2_16, file= "mf_cl0_vs2_16_dif_expression_analysis.csv")

mf_cl16_vs0_2 <- FindMarkers(lung, ident.1 = "16", ident.2 = c("0", "2"), test.use = "MAST", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.1)
mf_cl16_vs0_2$gene <- row.names(mf_cl16_vs0_2)
mf_cl16_vs0_2 <- join(mf_cl16_vs0_2, m, by="gene")
write.csv(mf_cl16_vs0_2, file= "mf_cl16_vs0_2_dif_expression_analysis.csv")
#==========================================================================================

levels(lung) <- c("17","20", "18", "15", "1",  "12", "21", "7",  "8",  "9",  "19", "11", "5",  "13", "14", "4",  "10", "16", "0",  "2",  "3", "6" )
new_order <- rev(levels(lung))

# analysis based of the top markers according to the avg_log2fc(10 genes) and then the difference in the percent of positives
top2 <- lung_cluster.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top2 <- top2 %>% group_by(cluster) %>% top_n(n = 5, wt = Dpct)
top2 <-top2[, c(7,6)]
top2 <- top2[!duplicated(top2[ , "gene"]),]


e <- top2[top2$cluster %in% new_order[1],]
for (i in 2:length(new_order)){
  e1 <- top2[top2$cluster %in% new_order[i],]
  e <- rbind(e, e1)
}

#      general,     STEM      ,    MF , monocyte,        DC     ,  neutrophil  ,     prolif       ,lymphoid,T-cells,             ILC               ,      NK                ,B-cells , MAST/BAS, megakaryocyte, platelets
e <- c("PTPRC", "CD34", "CD38", "C1QB", "S100A8", "CD1C", "FLT3","MPO", "ELANE",   "PCNA", "MKI67",  "LEF1", "CD3D", "CD2", "CD5", "GATA3", "IL1R1", "NCR1",  "CD8A", "NKG7", "CD79B", "CPA3"  , "ITGA2B"     ,   "TUBB1", e$gene)
e <- e[!duplicated(e)]
#6.5x24 inches 
DotPlot(lung, assay = "RNA", features = e, scale = TRUE) + RotatedAxis()

# analysis based of the top markers according to the avg_log2fc(10 genes) and then the difference in the percent of positives
top3 <- lung_cluster.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top3 <- top3 %>% group_by(cluster) %>% top_n(n = 5, wt = Dpct)
top3 <-top3[, c(7,6)]
top3 <- top3[!duplicated(top3[ , "gene"]),]

f <- top3[top3$cluster %in% new_order[1],]
for (i in 2:length(new_order)){
  f1 <- top3[top3$cluster %in% new_order[i],]
  f <- rbind(f, f1)
}
write.csv(f, file="Immune_Fig_B_gene_list.csv") 
#
f <- c("PTPRC", # general
       "CD34", "CD38", # STEM
       "CD163", "LYVE1", "C1QB", # macrophage
       "S100A8", "VCAN", "FCGR3B", # monocyte
       "HLA-DQA1", "CLEC9A", "CLEC4C", # DCs
       "MPO", "ELANE", # neutrophil
       "PCNA", "MKI67", 
       "LEF1", # lymphoid
       "CD3D", "CD3E", "CD3G" ,"CD5", # T-cells
       "GATA3", "IL1RL1", "KLRB1", # ILC2
       "IL1R1", "KIT", "IL23R", # ILC3
       "NKG7", "NCR1",  "CD8A", "KLRD1", "KLRF1", #NK 
       "CD19", "MS4A1", "CD79B", #B-cells
       "CPA3", "ENPP3", "HDC", #MAST/BAS
       "GP9", "ITGA2B", "TUBB1", # megakaryocytes
       f$gene)

known_genes <- c("PTPRC", # general
                 "CD34", "CD38", # STEM
                 "CD163", "LYVE1", "C1QB", # macrophage
                 "S100A8", "VCAN", "FCGR3B", # monocyte
                 "HLA-DQA1", "CLEC9A", "CLEC4C", # DCs
                 "MPO", "ELANE", # neutrophil
                 "PCNA", "MKI67", 
                 "LEF1", # lymphoid
                 "CD3D", "CD3E", "CD3G" ,"CD5", # T-cells
                 "GATA3", "IL1RL1", "KLRB1", # ILC2
                 "IL1R1", "KIT", "IL23R", # ILC3
                 "NKG7", "NCR1",  "CD8A", "KLRD1", "KLRF1", #NK 
                 "CD19", "MS4A1", "CD79B", #B-cells
                 "CPA3", "ENPP3", "HDC", #MAST/BAS
                 "GP9", "ITGA2B", "TUBB1") # megakaryocytes
  
  write.csv(known_genes, file="Immune_Fig_B_known_genes_list.csv")


f <- f[!duplicated(f)]
#6x26 inches 
DotPlot(lung, assay = "RNA", features = f, scale = TRUE) + RotatedAxis()

