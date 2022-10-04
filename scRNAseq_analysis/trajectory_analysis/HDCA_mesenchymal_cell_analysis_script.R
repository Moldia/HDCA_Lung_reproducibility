# this script has beeb tested in Ubuntu 20.04.2 LTS
setwd()

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
library(loomR)
library(harmony)
library(pheatmap)
library(RColorBrewer)
library(sctransform)
library(SeuratWrappers)
library(velocyto.R)
library(monocle3)
library(slingshot)
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
#remotes::install_github("czarnewski/niceRplots")
library(niceRplots)
library(igraph)
library(rafalib)

# increase the maximum of allowed RAM
options(future.globals.maxSize = 188000 * 1024^2)
#plan("multicore", workers = 4)

# load the mesenchymal dataset
load("mes_163k_org_data_with.velo.RData")
#=============================================

# the analusis is based on the SCTranform tutorial of https://satijalab.org/seurat/v3.0/integration.html
# create an object list based on the original datasets
lung <- mes
rm(mes)

DefaultAssay(lung) <-"RNA"

lung <- SetIdent(lung, value = lung@meta.data[["donor"]]) 
all.genes <- rownames(lung)

# human orthologues of stress induced genes because of tissue diggestion, according to https://pubmed.ncbi.nlm.nih.gov/32487174/.
# We have tested the analysis without regressing out stress induced genes with similar results but the following approach gave better clustering
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

lung@assays[["spliced"]] <- NULL
lung@assays[["unspliced"]] <- NULL
lung@assays[["ambiguous"]] <- NULL
VlnPlot(lung, features = "XIST")

# add the sex in metadata
lung.list <- SplitObject(lung, split.by = "donor")

rm(lung)

for (i in 1:length(lung.list)) { 
  lung.list[[i]] <- SCTransform(lung.list[[i]], vars.to.regress = c("stress_Features1"), variable.features.n = 5000, verbose = TRUE)
}

save(lung.list, file= "lung_list.RData")
load("lung_list.RData")
#=============================================
# Next, select features for downstream integration, and run PrepSCTIntegration, which ensures that all necessary Pearson residuals have been calculated.
lung.features <- SelectIntegrationFeatures(object.list = lung.list, nfeatures = 5000)
lung.list <- PrepSCTIntegration(object.list = lung.list, anchor.features = lung.features,
                                verbose = TRUE)

reference_datasets <- which(names(lung.list) == "XDD:385") 
lung.anchors <- FindIntegrationAnchors(object.list = lung.list, normalization.method = "SCT", 
                                       anchor.features = lung.features, reference = reference_datasets)

rm(lung.list)
#save(lung.anchors, file= "lung_anchors_reference.RData")

lung <- IntegrateData(anchorset = lung.anchors, normalization.method = "SCT", 
                      verbose = TRUE)
rm(lung.anchors, lung.features)

ages <- data.frame(lung@meta.data[["age"]])
row.names(ages) <- colnames(lung)
names(ages)[] <- "age_n"
ages$age_n<- as.numeric(str_replace(ages$age_n, "w", ""))
lung <- AddMetaData(lung, ages$age_n, col.name = "age_n")

save.image("mes_batch_integrated.RData")

#======================================================================
# change the active assay from "integrated" to "RNA"
DefaultAssay(lung) <- "RNA"

# to be able to do it with the correct values, we have to normalize and scale the RNA counts
lung <- NormalizeData(lung, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(lung)
lung <- ScaleData(lung, features = all.genes)

# remove the HBA1 high cells, because they are likely doublets with Red Blood cells 
lung1 <- subset(lung, subset = HBA1 >4)

lung <- lung[, !(colnames(lung) %in% colnames(lung1))]
all.genes <- rownames(lung)
rm(lung1,ages, reference_datasets)
save.image("mes_batch_integrated_filtered.RData")
#-----------------------------------------------------------
load("mes_batch_integrated_filtered.RData")

#Perform linear dimensional reduction
# Next we perform PCA on the scaled data. By default, only the previously determined variable features are used as input, but can be defined using 
# features argument if you wish to choose a different subset.
DefaultAssay(lung) <- "integrated"
lung <- RunPCA(lung, npcs = 50)

# set the number of principle components to use for dimension reduction
pcs=50
dim_n = 1:pcs
lung <- FindNeighbors(lung, dims = dim_n, reduction = "pca", assay = "integrated",k.param = 100, n.trees = 500)
# Run UMAP on a graph
lung <- FindClusters(lung, resolution = 1.2)
# run tSNE
lung <- RunTSNE(lung, reduction = "pca", dims = dim_n)
# Plot on tSNE
DimPlot(lung, reduction = "tsne", label = TRUE, pt.size = 1) + NoLegend()

#run UMAP (original settings)
lung <- RunUMAP(lung, 
                assay = "integrated",
                n.neighbors = 100L, 
                min.dist = 0.5, 
                dims = dim_n, 
                spread = 5,
                metric = "cosine",
                repulsion.strength = 0.2,
                negative.sample.rate = 5,
                n.epochs = NULL,
                seed.use = 42,
                learning.rate = 15,
                n.components = 2L,
                reduction.name = "umap2d",
                reduction.key='umap2d_')

#run UMAP (modified settings)
lung <- RunUMAP(lung, 
                assay = "integrated",
                n.neighbors = 100L, 
                min.dist = 0.5, 
                dims = dim_n, 
                spread = 17,
                metric = "cosine",
                repulsion.strength = 0.01,
                negative.sample.rate = 1,
                n.epochs = NULL,
                seed.use = 42,
                learning.rate = 15,
                n.components = 2L,
                reduction.name = "umap2d",
                reduction.key='umap2d_')

DimPlot(lung, reduction = "umap2d", label = TRUE, pt.size = 1)

save.image("mes_163K_after_umap2d_analysis_data.RData")
#=================================================================================
# run 3D UMAP 
lung <- RunUMAP(lung, 
                assay = "integrated",
                n.neighbors = 100L, 
                min.dist = 0.5, 
                dims = dim_n, 
                spread = 25L,
                metric = "cosine",
                repulsion.strength = 0.05,
                negative.sample.rate = 1,
                n.epochs = NULL,
                seed.use = 42,
                learning.rate = 15,
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

# Plot data
plot_ly(data = plot.data, 
        x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, 
        color = ~seurat_clusters, 
        colors = colors,
        type = "scatter3d", 
        mode = "markers",
        marker = list(size = 1.5, width=1.5), # controls size of points
        text=~label, #This is that extra column we made earlier for which we will use for cell ID
        hoverinfo="text") #When you visualize your plotly object, hovering your mouse pointer over a point shows cell names


# plot the 3D UMAP using plot3D package, theta and phi change the orinetation of the graph 
scatter3D(x = plot.data[,1], y = plot.data[,2], z = plot.data[,3], 
          xlab= "UMAP-1", ylab= "UMAP-2", zlab= "UMAP-3",
          colvar= as.integer(plot.data[,4]),
          col = colors,
          bty = "b2", 
          pch = 20, cex = 0.2, ticktype = "detailed",
          theta = 125, phi = 15,
          clab= "Clusters"
)


pdf("mes_163k_umap3d_clusters.pdf", width = 11.69, height = 8.27, paper = "special")
pdf <-scatter3D(x = plot.data[,1], y = plot.data[,2], z = plot.data[,3], 
                xlab= "UMAP-1", ylab= "UMAP-2", zlab= "UMAP-3",
                colvar= as.integer(plot.data[,4]),
                col = colors,
                bty = "b2", 
                pch = 20, cex = 0.2, ticktype = "detailed",
                theta = 125, phi = 15,
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
        marker = list(size = 1.5, width=1.5), # controls size of points
        text=~label, #This is that extra column we made earlier for which we will use for cell ID
        hoverinfo="text") #When you visualize your plotly object, hovering your mouse pointer over a point shows cell names


# plot the 3D UMAP using plot3D package, theta and phi change the orinetation of the graph 
scatter3D(x = plot.data_meta[,1], y = plot.data_meta[,2], z = plot.data_meta[,3], 
          xlab= "UMAP-1", ylab= "UMAP-2", zlab= "UMAP-3",
          colvar= as.integer(as.factor(plot.data_meta[,4])),
          col = colors1,
          bty = "b2", 
          pch = 20, cex = 0.2, ticktype = "detailed",
          theta = 125, phi = 15,
          clab= metadata_to_plot
)

pdf(paste0("mes_163k_umap3d_",metadata_to_plot,".pdf"), width = 11.69, height = 8.27, paper = "special")
pdf <-scatter3D(x = plot.data_meta[,1], y = plot.data_meta[,2], z = plot.data_meta[,3], 
                xlab= "UMAP-1", ylab= "UMAP-2", zlab= "UMAP-3",
                colvar= as.integer(as.factor(plot.data_meta[,4])),
                col = colors1,
                bty = "b2", 
                pch = 20, cex = 0.2, ticktype = "detailed",
                theta = 125, phi = 15,
                clab= metadata_to_plot
)
plot(pdf)
dev.off()
rm(pdf)

#=======================================================================================================
# run the function again setting the resolution according to the above results
DefaultAssay(lung) <- "integrated"

# export Umap plots with the identified clusters, the age and sample name
pdf("mes_cells_umap2d_clusters.pdf",width=11.69,height=8.27,paper="special")
mes_cells_umap2d_clusters <- DimPlot(lung, reduction = "umap2d", cols =  "alphabet2", label = TRUE, pt.size = 0.5)
plot(mes_cells_umap2d_clusters)
rm(mes_cells_umap2d_clusters)
dev.off()

pdf("mes_cells_umap2d_age.pdf",width=11.69,height=8.27,paper="special")
mes_cells_umap2d_age <- DimPlot(lung, reduction = "umap2d", cols =  "alphabet2", group.by = "age", pt.size = 0.5)
plot(mes_cells_umap2d_age)
rm(mes_cells_umap2d_age)
dev.off()

pdf("mes_cells_umap2d_SampleName.pdf",width=11.69,height=8.27,paper="special")
mes_cells_umap2d_SampleName <- DimPlot(lung, reduction = "umap2d", cols =  "alphabet2", group.by = "SampleName", pt.size = 0.5)
plot(mes_cells_umap2d_SampleName)
rm(mes_cells_umap2d_SampleName)
dev.off()

pdf("mes_cells_umap2d_donor.pdf",width=11.69,height=8.27,paper="special")
mes_cells_umap2d_donor <- DimPlot(lung, reduction = "umap2d", cols =  "alphabet2", group.by = "donor", pt.size = 0.5)
plot(mes_cells_umap2d_donor)
rm(mes_cells_umap2d_donor)
dev.off()

pdf("mes_cells_umap2d_163k_anno.pdf",width=11.69,height=8.27,paper="special")
mes_cells_umap2d_163k_anno <- DimPlot(lung, reduction = "umap2d", cols =  "alphabet2", group.by = "anno_163k", pt.size = 0.5, label = TRUE)
plot(mes_cells_umap2d_163k_anno)
rm(mes_cells_umap2d_163k_anno)
dev.off()

pdf("mes_cells_umap2d_chemistry.pdf",width=11.69,height=8.27,paper="special")
mes_cells_umap2d_chemistry <- DimPlot(lung, reduction = "umap2d", group.by = "chemistry", pt.size = 0.5)
plot(mes_cells_umap2d_chemistry)
rm(mes_cells_umap2d_chemistry)
dev.off()

pdf("mes_cells_umap2d_features.pdf",width=11.69,height=8.27,paper="special")
mes_cells_umap2d_features <- FeaturePlot(lung, reduction = "umap2d",  features = "nFeature_RNA", pt.size = 0.5)
plot(mes_cells_umap2d_features)
rm(mes_cells_umap2d_features)
dev.off()

pdf("mes_cells_umap2d_percent_mt.pdf",width=11.69,height=8.27,paper="special")
mes_cells_umap2d_percent_mt <- FeaturePlot(lung, reduction = "umap2d",  features = "MT_ratio", pt.size = 0.5)
plot(mes_cells_umap2d_percent_mt)
rm(mes_cells_umap2d_percent_mt)
dev.off()

pdf("mes_cells_tsne_clusters.pdf",width=11.69,height=8.27,paper="special")
s10X140_cells_tsne_clusters <- DimPlot(lung, reduction = "tsne", cols =  "alphabet2", label = TRUE, pt.size = 0.5)
plot(s10X140_cells_tsne_clusters)
rm(s10X140_cells_tsne_clusters)
dev.off()

pdf("mes_cells_tsne_age.pdf",width=11.69,height=8.27,paper="special")
mes_cells_tsne_age <- DimPlot(lung, reduction = "tsne", cols =  "alphabet2", group.by = "age", pt.size = 0.5)
plot(mes_cells_tsne_age)
rm(mes_cells_tsne_age)
dev.off()

pdf("mes_cells_tsne_SampleName.pdf",width=11.69,height=8.27,paper="special")
mes_cells_tsne_SampleName <- DimPlot(lung, reduction = "tsne", cols =  "alphabet2", group.by = "SampleName", pt.size = 0.5)
plot(mes_cells_tsne_SampleName)
rm(mes_cells_tsne_SampleName)
dev.off()

pdf("mes_cells_tsne_donor.pdf",width=11.69,height=8.27,paper="special")
mes_cells_tsne_donor <- DimPlot(lung, reduction = "tsne", cols =  "alphabet2", group.by = "donor", pt.size = 0.5)
plot(mes_cells_tsne_donor)
rm(mes_cells_tsne_donor)
dev.off()

pdf("mes_cells_tsne_chemistry.pdf",width=11.69,height=8.27,paper="special")
mes_cells_tsne_chemistry <- DimPlot(lung, reduction = "tsne", group.by = "chemistry", pt.size = 0.5)
plot(mes_cells_tsne_chemistry)
rm(mes_cells_tsne_chemistry)
dev.off()

pdf("mes_cells_tsne_percent_mt.pdf",width=11.69,height=8.27,paper="special")
mes_cells_tsne_percent_mt <- FeaturePlot(lung, reduction = "tsne", features = "MT_ratio", pt.size = 0.5)
plot(mes_cells_tsne_percent_mt)
rm(mes_cells_tsne_percent_mt)
dev.off()

pdf("mes_cells_tsne_features.pdf",width=11.69,height=8.27,paper="special")
mes_cells_tsne_features <- FeaturePlot(lung, reduction = "tsne", features = "nFeature_RNA", pt.size = 0.5)
plot(mes_cells_tsne_features)
rm(mes_cells_tsne_features)
dev.off()

#=============================================================================================
# create *.tiff files with individual clusters being colored with red

for (i in 0:(nlevels(lung@active.ident)-1)) {
  cluster <- WhichCells(lung, ident = i)
  name = sprintf("mes_umap2d_cluster%s.tiff", i)
  tiff(name, width = 707, height = 500, compression = "lzw")
  image <- DimPlot(lung, reduction = "umap2d",  cells.highlight = cluster, sizes.highlight= 0.5, pt.size = 0.5) + NoLegend()
  plot(image)
  rm(image)
  dev.off()
}

for (i in 0:(nlevels(lung@active.ident)-1)) {
  cluster <- WhichCells(lung, ident = i)
  name = sprintf("mes_tsne_cluster%s.tiff", i)
  tiff(name, width = 707, height = 500, compression = "lzw")
  image <- DimPlot(lung, reduction = "tsne", cells.highlight = cluster, sizes.highlight= 3, pt.size = 0.5) + NoLegend()
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
                  pch = 20, cex = 0.2, ticktype = "detailed",
                  theta = 125, phi = 15,
                  clab= paste0("Cluster", i)
  )
  plot(pdf)
  dev.off()
  rm(pdf)
}
#==============================================================================================
# differential expression analysis to identify specific or selective genes for each cluster

# change the active assay from "integrated" to "RNA"
DefaultAssay(lung) <- "RNA"

# to be able to do it with the correct values, we have to normalize and scale the RNA counts
lung <- NormalizeData(lung, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(lung)
lung <- ScaleData(lung, features = all.genes)

# Run differential expression analysis
min_cluster <- min(table(lung@active.ident))
lung_cluster.markers <- FindAllMarkers(lung, test.use = "MAST", min.diff.pct = 0.1, max.cells.per.ident = min_cluster, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.1)
# calculate the percentage difference of the positive cells in the indicated cluster, in comparison to the others
lung_cluster.markers$Dpct <- (lung_cluster.markers$pct.1 - lung_cluster.markers$pct.2)

c <-subset(lung, idents =0 )
set.seed(123)
sampled.c <- sample(x = colnames(c), size = min_cluster, replace = F)
c <- c[, colnames(c) %in% sampled.c]
counts <- as.matrix(c[["RNA"]]@data)
m <- data.frame(rowMeans(counts))
gene <- c@assays[["RNA"]]@data@Dimnames[[1]]
m <- cbind(gene, m)
names(m)[2] <- paste0("cluster", "0") 

for (i in 1:(nlevels(lung@active.ident)-1)) {
  c1 <-subset(lung, idents =i )
  set.seed(123)
  sampled.c1 <- sample(x = colnames(c1), size = min_cluster, replace = F)
  c1 <- c1[, colnames(c1) %in% sampled.c1]
  counts1 <- as.matrix(c1[["RNA"]]@data)
  m1 <- data.frame(rowMeans(counts1))
  names(m1)[1] <- paste0("cluster", (i))
  m <- cbind(m, m1)
}
rm(c, counts,sampled.c, counts1, c1, m1, sampled.c1)


lung_cluster.markers <-join(lung_cluster.markers, m, by="gene")
write.csv(lung_cluster.markers, file="mes_sct_lung_markers.csv")


# filter out insignificant genes
lung_cluster.markers <- lung_cluster.markers[lung_cluster.markers$p_val_adj<0.001, ]

# obtain the top 25 genes of each cluster, based on the average logFC
top25 <- lung_cluster.markers %>% group_by(cluster) %>% top_n(n = 25, wt = avg_log2FC)
write.csv(top25, file="mes_sct_lung_cluster_top25_markers.csv")


top5 <- lung_cluster.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
write.csv(top5, file="mes_sct_lung_cluster_top5_markers.csv")


#-------------------------------------------------------------------------------------------
# The transcription factor and co-factor lists were downloaded from the http://bioinfo.life.hust.edu.cn/AnimalTFDB/#!/download
# The kinase list "protein_class_kinases.tsv" was downloaded from https://www.proteinatlas.org/search/protein_class:kinases
# The cd marker list "protein_class_CD.tsv" was downloaded from https://www.proteinatlas.org/search/protein_class%3ACD+markers

# import the list of transcription factors 
TFs <- read.csv("TFs.csv")

# create a new file by merging the TF-list and the enriched genes in each cluster, to identify the TFs, only
lung_cluster_tfs <- join(lung_cluster.markers, TFs, by="gene", type="inner")
top10_cluster_tfs <- lung_cluster_tfs %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
# export the list of the top 10 TFs based on average expression difference
write.csv(top10_cluster_tfs, file="mes_sct__top10_cluster_tfs.csv")

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
write.csv(top10_cluster_cotfs, file="mes_sct__top10_cluster_cotfs.csv")


# import the list of kinases
kinases <- data.frame(read.delim("protein_class_kinases.tsv"))
kinases <- kinases[1:10]

# create a new file by merging the kinase-list and the enriched genes in each cluster, to identify the kinases, only
lung_cluster_kinases <- join(lung_cluster.markers, kinases, by="gene", type="inner")
top10_cluster_kinases <- lung_cluster_kinases %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
# export the list of the top 10 kinases based on average expression difference
write.csv(top10_cluster_kinases, file="mes_sct__top10_cluster_kinases.csv")


# import the list of CD markers
cds <- data.frame(read.delim("protein_class_CD.tsv"))
cds <- cds[1:10]

# create a new file by merging the CD-list and the enriched genes in each cluster, to identify the surface markers, only
lung_cluster_cds <- join(lung_cluster.markers, cds, by="gene", type="inner")
top10_cluster_cds <- lung_cluster_cds %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
# export the list of the top 10 CDs based on average expression difference
write.csv(top10_cluster_cds, file="top10_cluster_cds.csv")

#----------------------------------------------------------------------------------------------------------------------
# Run differential expression analysis with ROC
DefaultAssay(lung) <-"RNA"
lung_roc_cluster.markers <- FindAllMarkers(lung, test.use = "roc", max.cells.per.ident = min_cluster, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.1)
# calculate the percentage difference of the positive cells in the indicated cluster, in comparison to the others
lung_roc_cluster.markers$Dpct <- (lung_roc_cluster.markers$pct.1 - lung_roc_cluster.markers$pct.2)

lung_roc_cluster.markers <-join(lung_roc_cluster.markers, m, by="gene")
write.csv(lung_roc_cluster.markers, file="mes_163k_roc_markers.csv")

lung_roc_cluster.markers <- lung_roc_cluster.markers[lung_roc_cluster.markers$Dpct >0.1, ]

# obtain the top 25 genes of each cluster, based on the power
top25_roc <- lung_roc_cluster.markers %>% group_by(cluster) %>% top_n(n = 25, wt = power)
write.csv(top25_roc, file="mes_163k_roc_cluster_top25_markers.csv")


top5_roc <- lung_roc_cluster.markers %>% group_by(cluster) %>% top_n(n = 5, wt = power)
write.csv(top5_roc, file="mes_163k_roc_cluster_top5_markers.csv")


top10_roc <- lung_roc_cluster.markers %>% group_by(cluster) %>% top_n(n = 10, wt = power)
write.csv(top10_roc, file="mes_163k_roc_cluster_top10_markers.csv")

#-------------------------------------------------------------------------------------------
# The transcription factor and co-factor lists were downloaded from the http://bioinfo.life.hust.edu.cn/AnimalTFDB/#!/download
# The kinase list "protein_class_kinases.tsv", cd marker list "protein_class_CD.tsv" and "protein_class_secreted.tsv" were downloaded from https://www.proteinatlas.org/


# import the list of transcription factors 
TFs <- read.csv("TFs.csv")

# create a new file by merging the TF-list and the enriched genes in each cluster, to identify the TFs, only
lung_roc_cluster_tfs <- join(lung_roc_cluster.markers, TFs, by="gene", type="inner")
top10_roc_cluster_tfs <- lung_roc_cluster_tfs %>% group_by(cluster) %>% top_n(n = 10, wt = power)
# export the list of the top 10 TFs based on average expression difference
write.csv(top10_roc_cluster_tfs, file="mes_163k_roc_top10_cluster_tfs.csv")


# import the list of transcription cofactors 
coTFs <- read.csv("TF_cofactors.csv")

# create a new file by merging the coTF-list and the enriched genes in each cluster, to identify the coTFs, only
lung_roc_cluster_cotfs <- join(lung_roc_cluster.markers, coTFs, by="gene", type="inner")
top10_roc_cluster_cotfs <- lung_roc_cluster_cotfs %>% group_by(cluster) %>% top_n(n = 10, wt = power)
# export the list of the top 10 coTFs based on average expression difference
write.csv(top10_roc_cluster_cotfs, file="mes_163k_roc_top10_cluster_cotfs.csv")


# import the list of kinases
kinases <- data.frame(read.delim("protein_class_kinases.tsv"))
kinases <- kinases[1:10]

# create a new file by merging the kinase-list and the enriched genes in each cluster, to identify the kinases, only
lung_roc_cluster_kinases <- join(lung_roc_cluster.markers, kinases, by="gene", type="inner")
top10_roc_cluster_kinases <- lung_roc_cluster_kinases %>% group_by(cluster) %>% top_n(n = 10, wt = power)
# export the list of the top 10 kinases based on average expression difference
write.csv(top10_roc_cluster_kinases, file="mes_163k_roc_top10_cluster_kinases.csv")


# import the list of CD markers
cds <- data.frame(read.delim("protein_class_CD.tsv"))
cds <- cds[1:10]

# create a new file by merging the CD-list and the enriched genes in each cluster, to identify the surface markers, only
lung_roc_cluster_cds <- join(lung_roc_cluster.markers, cds, by="gene", type="inner")
top10_roc_cluster_cds <- lung_roc_cluster_cds %>% group_by(cluster) %>% top_n(n = 10, wt = power)
# export the list of the top 10 CDs based on average expression difference
write.csv(top10_roc_cluster_cds, file="top10_roc_cluster_cds.csv")

#============================================
# create and export a table with the cell names and their corresponding cluster 
mes_sct_cells_clusters <- data.frame(lung@active.ident)
write.csv(mes_sct_cells_clusters, file="mes_sct_cells_clusters_sct.csv")

Age <- data.frame(lung@meta.data[["age"]])
experiment <- data.frame(lung@meta.data[["orig.ident"]])
cell_name <- row.names(mes_sct_cells_clusters)
cluster_list <- cbind(cell_name, mes_sct_cells_clusters)
cluster_list <- cbind(cluster_list, Age)
cluster_list <- cbind(cluster_list, experiment)
cluster_list <- setNames(cluster_list, c("cell", "cluster", "age", "orig_ident"))
write.csv(cluster_list, file="mes_sct_metadata.csv")

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
write.csv(abundance_matrix, file="age_abundance_in_clusters.csv")


# calculate the % of cells of each time-point in the clusters
cluster_size <- rowSums(abundance_matrix)
abundance_matrix_percentages <- (abundance_matrix*100)/cluster_size
write.csv(abundance_matrix_percentages, file="age_abundance_percentages_in_clusters.csv")

save.image("mes_163K_after_MAST_and_roc_analysis_data.RData")

#=====================================================================================================
# differential expression analysis of specific clusters

chondro_clusters <- c("20", "12")
for (i in 1:length(chondro_clusters)){
  a <- FindMarkers(lung, ident.1 = chondro_clusters[[i]], ident.2 = chondro_clusters[- i], test.use = "MAST",max.cells.per.ident = min_cluster, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.1)
  gene <-row.names(a)
  a <-cbind(gene, a)
  a <-join(a, m, by="gene")
  write.csv(a, file = paste0("chondro_cluster", chondro_clusters[[i]], "vs_other_chondro.csv"))
  rm(a, gene)          
}

#===============================================================================================================
#===============================================================================================
# Run PAGA-plot analysis according to Paulo script
#Create graph abstraction
DefaultAssay(lung) <- "integrated"
g2 <- graph_abstraction( lung , red = "umap2d" , clustering = "seurat_clusters" , 
                         graph = "integrated_snn", cutoff = 0)

# set colors
n_clusters <- nlevels(lung$seurat_clusters)
colors <- DiscretePalette(n_clusters, palette = "alphabet2")

n_vars <- nlevels(as.factor(lung@meta.data[["donor"]]))
colors1 <- DiscretePalette(n_vars, palette = "alphabet2")

n_vars <- nlevels(as.factor(lung@meta.data[["age"]]))
colors2 <- DiscretePalette(n_vars, palette = "alphabet2")

#Plot graph abstraction onto "umap_harmony" and visualize some metadata too
png(filename = "Graph_abstraction.png", width = 1600*4,height = 1600*4,res = 400)
mypar(4,4)
plot_meta(x = lung,red = "umap2d",feat = "age",frame=F,col = colors2,label = T)
plot_meta(x = lung,red = "umap2d",feat = "donor",frame=F,col = colors1,label = T)
plot_meta(x = lung,red = "umap2d",feat = "seurat_clusters",frame=F,col = colors,label = T)

for(i in 10^-(0:10) ){
  centroids <- g2[g2$s==g2$p.Var,]
  g <- g2[ as.numeric(g2$p.Freq) > i,]
  g <- g[ as.numeric(g$p.Freq) > 0,]
  g <- g[ g$s != g$p.Var,]
  g$norm.Freq <- g$p.Freq / max(g$p.Freq)
  plot_meta(x = lung,red = "umap2d",feat = "seurat_clusters",frame=T,col = colors,label = F,main = paste0("cutoff=",i))
  gpal <- colorRampPalette(c("grey20","black"))(20)
  apply(g,1,function(x){
    lines(as.numeric(x[4:5]),as.numeric(x[6:7]),
          lwd=5*as.numeric(x[8]),
          col=gpal[round(as.numeric(x[8])*18)+1]) })
  points(centroids[,4],centroids[,6],
         bg='white',
         pch=21,cex=1)
  text(centroids[,4],centroids[,6],labels = centroids[,1],cex=.3)
}
dev.off()

#Transform graph into probabilities
probs <- matrix( g2$p.Freq , sqrt( nrow(g2) ) , sqrt( nrow(g2) ))
colnames(probs) <- unique(g2$p.Var)
rownames(probs) <- unique(g2$p.Var)
g3 <- graph_from_adjacency_matrix( probs , weighted = T , diag = T , mode="undirected")

# Then you can use this function to find the shortest distance from cluster 0 to cluster 6 (for example)
igraph::shortest_paths( g3, from = "4", to = "10" , weights = 1/E(g3)$weight )


#
for(i in 2*10^-(5) ){
  centroids <- g2[g2$s==g2$p.Var,]
  g <- g2[ as.numeric(g2$p.Freq) > i,]
  g <- g[ as.numeric(g$p.Freq) > 0,]
  g <- g[ g$s != g$p.Var,]
  g$norm.Freq <- g$p.Freq / max(g$p.Freq)
  plot_meta(x = lung,red = "umap2d",feat = "seurat_clusters",frame=T,col = colors,label = F,main = paste0("cutoff=",i))
  gpal <- colorRampPalette(c("grey20","black"))(20)
  apply(g,1,function(x){
    lines(as.numeric(x[4:5]),as.numeric(x[6:7]),
          lwd=5*as.numeric(x[8]),
          col=gpal[round(as.numeric(x[8])*18)+1]) })
  points(centroids[,4],centroids[,6],
         bg='white',
         pch=21,cex=1)
  text(centroids[,4],centroids[,6],labels = centroids[,1],cex=.3)
}

dev.off()
