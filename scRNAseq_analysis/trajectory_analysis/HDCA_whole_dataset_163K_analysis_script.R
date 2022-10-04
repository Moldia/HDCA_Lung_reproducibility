# this script has beeb tested in Ubuntu 20.04.2 LTS
setwd()
library(plyr)
library(dplyr)
library(Seurat)
library(future)
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
library(future)
library(glmGamPoi)

# increase the maximum of allowed RAM
options(future.globals.maxSize = 250000 * 1024^2)

#check the plan for multiprocess
plan()

#change the plan to allow parallel
plan("multicore", workers = 4)
plan()
#=========================================================================================
load("/home/samakovlis/Desktop/163K_cells_100pcs/raw/s10X140w85_filtered_doublets.RData")
DefaultAssay(s10X140w85) <-"RNA"
s10X140w85@assays[["SCT"]] <- NULL
load("/home/samakovlis/Desktop/163K_cells_100pcs/raw/s10X256w14_filtered_doublets.RData")
DefaultAssay(s10X256w14) <-"RNA"
s10X256w14@assays[["SCT"]] <- NULL
load("/home/samakovlis/Desktop/163K_cells_100pcs/raw/s10X263w12_filtered_doublets.RData")
DefaultAssay(s10X263w12) <-"RNA"
s10X263w12@assays[["SCT"]] <- NULL
load("/home/samakovlis/Desktop/163K_cells_100pcs/raw/s10X151w6_filtered_doublets.RData")
DefaultAssay(s10X151w6) <-"RNA"
s10X151w6@assays[["SCT"]] <- NULL
load("/home/samakovlis/Desktop/163K_cells_100pcs/raw/s10X272w115_filtered_doublets.RData")
DefaultAssay(s10X272w115) <-"RNA"
s10X272w115@assays[["SCT"]] <- NULL
load("/home/samakovlis/Desktop/163K_cells_100pcs/raw/s10X158w8_filtered_doublets.RData")
DefaultAssay(s10X158w8) <-"RNA"
s10X158w8@assays[["SCT"]] <- NULL
load("/home/samakovlis/Desktop/163K_cells_100pcs/raw/s10X274w115_filtered_doublets.RData")
DefaultAssay(s10X274w115) <-"RNA"
s10X274w115@assays[["SCT"]] <- NULL
load("/home/samakovlis/Desktop/163K_cells_100pcs/raw/s10X166w85_filtered_doublets.RData")
DefaultAssay(s10X166w85) <-"RNA"
s10X166w85@assays[["SCT"]] <- NULL
load("/home/samakovlis/Desktop/163K_cells_100pcs/raw/s10X285w10_filtered_doublets.RData")
DefaultAssay(s10X285w10) <-"RNA"
s10X285w10@assays[["SCT"]] <- NULL
load("/home/samakovlis/Desktop/163K_cells_100pcs/raw/s10X180w5_filtered_doublets.RData")
DefaultAssay(s10X180w5) <-"RNA"
s10X180w5@assays[["SCT"]] <- NULL
load("/home/samakovlis/Desktop/163K_cells_100pcs/raw/s10X289w6_filtered_doublets.RData")
DefaultAssay(s10X289w6) <-"RNA"
s10X289w6@assays[["SCT"]] <- NULL
load("/home/samakovlis/Desktop/163K_cells_100pcs/raw/s10X189w12_filtered_doublets.RData")
DefaultAssay(s10X189w12) <-"RNA"
s10X189w12@assays[["SCT"]] <- NULL
load("/home/samakovlis/Desktop/163K_cells_100pcs/raw/s10X299_1w7_filtered_doublets.RData")
DefaultAssay(s10X299_1w7) <-"RNA"
s10X299_1w7@assays[["SCT"]] <- NULL
load("/home/samakovlis/Desktop/163K_cells_100pcs/raw/s10X214w13_filtered_doublets.RData")
DefaultAssay(s10X214w13) <-"RNA"
s10X214w13@assays[["SCT"]] <- NULL
load("/home/samakovlis/Desktop/163K_cells_100pcs/raw/s10X299_2w7_filtered_doublets.RData")
DefaultAssay(s10X299_2w7) <-"RNA"
s10X299_2w7@assays[["SCT"]] <- NULL
load("/home/samakovlis/Desktop/163K_cells_100pcs/raw/s10X247w13_filtered_doublets.RData")
DefaultAssay(s10X247w13) <-"RNA"
s10X247w13@assays[["SCT"]] <- NULL
load("/home/samakovlis/Desktop/163K_cells_100pcs/raw/s10X303w55_filtered_doublets.RData")
DefaultAssay(s10X303w55) <-"RNA"
s10X303w55@assays[["SCT"]] <- NULL
load("/home/samakovlis/Desktop/163K_cells_100pcs/raw/s10X253w14_filtered_doublets.RData")
DefaultAssay(s10X253w14) <-"RNA"
s10X253w14@assays[["SCT"]] <- NULL

lung <- merge(
  x = s10X140w85,
  y = c(s10X263w12, s10X151w6, s10X272w115, s10X158w8, s10X274w115, s10X166w85, s10X285w10, s10X180w5, s10X289w6, s10X189w12,
        s10X299_1w7, s10X214w13, s10X299_2w7, s10X247w13, s10X303w55, s10X253w14, s10X256w14),
  merge.data = TRUE,
  project = "SeuratProject"
)

rm(s10X140w85, s10X263w12, s10X151w6, s10X272w115, s10X158w8, s10X274w115, s10X166w85, s10X285w10, s10X180w5, s10X289w6, s10X189w12,
   s10X299_1w7, s10X214w13, s10X299_2w7, s10X247w13, s10X303w55, s10X253w14, s10X256w14)

save(lung, file = "lung_all_filtered_cells_org_velo.RData")

lung
#163236 cells

# filter out the genes with less than 100 UMIs in all analyzed cells
genes_non0 <- data.frame(rowSums(lung@assays[["RNA"]]@counts))
genes <- rownames(genes_non0)
genes_non0 <- cbind(genes, genes_non0)
genes_non01 <- genes_non0[genes_non0$rowSums.lung.assays...RNA....counts. >100,]

lung <- subset(lung, features = genes_non01$genes)

rm(genes, genes_non0, genes_non01)

#====================================================================================
# create an object without velocyto values

lung@assays[["spliced"]] <-NULL
lung@assays[["unspliced"]] <-NULL
lung@assays[["ambiguous"]] <-NULL
save(lung, file = "lung_all_filtered_cells_org_without_velo.RData")
#=============================================
load("lung_all_filtered_cells_org_without_velo.RData")

# create an object list based on the original datasets
lung.list <- SplitObject(lung, split.by = "orig_ident")
# remove the "lung"object to save memory
rm(lung)

# run the SCTransform according to https://satijalab.org/seurat/articles/integration_rpca.html
for (i in 1:length(lung.list)) {
  lung.list[[i]] <- SCTransform(lung.list[[i]], verbose = TRUE, variable.features.n= 5000, return.only.var.genes= FALSE)
}

features <- SelectIntegrationFeatures(object.list = lung.list, nfeatures = 5000)
save.image("lung_list_after_sctransform_5000genes.RData")

lung.list <- PrepSCTIntegration(object.list = lung.list, anchor.features = features)
save.image("lung_list_after_sct_integration_5000genes.RData")

# run the step with sequential plan, because it gives error with multicore
#load("lung_list_after_sct_integration_5000genes.RData")
reference_datasets <- which(naall_cells(lung.list) == "10X253") 
lung.anchors <- FindIntegrationAnchors(object.list = lung.list, normalization.method = "SCT", 
                                       anchor.features = features, reference = reference_datasets)
rm(lung.list)
save(lung.anchors, file= "lung_anchors_reference.RData")
lung <- IntegrateData(anchorset = lung.anchors, normalization.method = "SCT")
rm(lung.anchors)

ages <- data.frame(lung@meta.data[["age"]])
rownames(ages) <- colnames(lung)
colnames(ages)[] <- "age_n"
ages$age_n<- as.numeric(str_replace(ages$age_n, "w", ""))
lung <- AddMetaData(lung, ages$age_n, col.name = "age_n")

save.image("lung_integrated_reference.RData")

#Perform linear dimensional reduction
# Next we perform PCA on the scaled data. By default, only the previously determined variable features are used as input, but can be defined using 
# features argument if you wish to choose a different subset.
lung <- RunPCA(lung, npcs = 100)

pcs=100
dim_n = 1:pcs
lung <- FindNeighbors(lung, dims = dim_n, k.param = 20)

# The FindClusters function iteratively groups cells together and contains a resolution parameter that sets the granularity of the downstream clustering, 
lung <- FindClusters(lung, resolution = 0.8)

# obtain cluster colors, similar to 2D default Seurat plot
n_clusters <- nlevels(lung$seurat_clusters)
#colors <- DiscretePalette(n_clusters, palette = "alphabet2")
cols1 <- DiscretePalette(26, palette = "alphabet2")
cols2 <- DiscretePalette((n_clusters-26), palette = "alphabet")
colors <- c(cols1, cols2)

#Run non-linear dimensional reduction (UMAP/tSNE)
lung <- RunTSNE(lung, dims = dim_n)
# Plot tSNE
DimPlot(lung, reduction = "tsne", label = TRUE, pt.size = 0.5, cols = colors) + NoLegend()
#DimPlot(lung, reduction = "tsne", label = FALSE, pt.size = 0.5)

#run UMAP
lung <- RunUMAP(lung, 
                assay = "integrated",
                n.neighbors = 30L, 
                min.dist = 0.3, 
                dims = dim_n, 
                spread = 10L,
                metric = "cosine",
                repulsion.strength = 0.25,
                negative.sample.rate = 5,
                n.epochs = NULL,
                seed.use = 42,
                learning.rate = 5,
                n.components = 2L,
                reduction.name = "umap2d",
                reduction.key='umap2d_')
#plot UMAP
DimPlot(lung, reduction = "umap2d", label = TRUE, cols = colors, pt.size = 1)
#================================================================================
# run 3d UMAP
DefaultAssay(lung) <- "integrated"

lung <- RunUMAP(lung, 
                assay = "integrated",
                n.neighbors = 30L, 
                min.dist = 0.3, 
                dims = dim_n, 
                spread = 5L,
                metric = "cosine",
                repulsion.strength = 0.1,
                negative.sample.rate = 5,
                n.epochs = NULL,
                seed.use = 42,
                learning.rate = 5,
                n.components = 3L)

# plot of 3D UMAP-plots
umap_1 <- lung[["umap"]]@cell.embeddings[,1]
umap_2 <- lung[["umap"]]@cell.embeddings[,2]
umap_3 <- lung[["umap"]]@cell.embeddings[,3]

# Visualize what headings are called so that you can extract them to form a dataframe
Embeddings(object = lung, reduction = "umap")

# Prepare a dataframe for cell plotting
plot.data <- FetchData(object = lung, vars = c("UMAP_1", "UMAP_2", "UMAP_3", "seurat_clusters"))

# Make a column of row name identities (these will be your cell/barcode naall_cells)
plot.data$label <- paste(rownames(plot.data))

dev.off()

# Plot your data, in this example my Seurat object had 21 clusters (0-20)
plot_ly(data = plot.data, 
        x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, 
        color = ~seurat_clusters, 
        colors = colors,
        type = "scatter3d", 
        mode = "markers",
        marker = list(size = 0.5, width=0.5), # controls size of points
        text=~label, #This is that extra column we made earlier for which we will use for cell ID
        hoverinfo="text") #When you visualize your plotly object, hovering your mouse pointer over a point shows cell naall_cells


# plot the 3D UMAP using plot3D package, theta and phi change the orinetation of the graph 
scatter3D(x = plot.data[,1], y = plot.data[,2], z = plot.data[,3], 
          xlab= "UMAP-1", ylab= "UMAP-2", zlab= "UMAP-3",
          colvar= as.integer(plot.data[,4]),
          col = colors,
          bty = "b2", 
          pch = 20, cex = 0.2, ticktype = "detailed",
          theta = 60, phi = 5,
          clab= "Clusters"
)


pdf("all_cells_163k_umap3d_clusters.pdf", width = 11.69, height = 8.27, paper = "special")
pdf <-scatter3D(x = plot.data[,1], y = plot.data[,2], z = plot.data[,3], 
                xlab= "UMAP-1", ylab= "UMAP-2", zlab= "UMAP-3",
                colvar= as.integer(plot.data[,4]),
                col = colors,
                bty = "b2", 
                pch = 20, cex = 0.2, ticktype = "detailed",
                theta = 60, phi = 5,
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

# Make a column of row name identities (these will be your cell/barcode naall_cells)
plot.data_meta$label <- paste(rownames(plot.data_meta))

# Plot your data, in this example my Seurat object had 21 clusters (0-20)
plot_ly(data = plot.data_meta, 
        x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, 
        color = plot.data_meta[,4], 
        colors = colors1,
        type = "scatter3d", 
        mode = "markers",
        marker = list(size = 0.5, width=0.5), # controls size of points
        text=~label, #This is that extra column we made earlier for which we will use for cell ID
        hoverinfo="text") #When you visualize your plotly object, hovering your mouse pointer over a point shows cell naall_cells


# plot the 3D UMAP using plot3D package, theta and phi change the orinetation of the graph 
scatter3D(x = plot.data_meta[,1], y = plot.data_meta[,2], z = plot.data_meta[,3], 
          xlab= "UMAP-1", ylab= "UMAP-2", zlab= "UMAP-3",
          colvar= as.integer(as.factor(plot.data_meta[,4])),
          col = colors1,
          bty = "b2", 
          pch = 20, cex = 0.2, ticktype = "detailed",
          theta = 60, phi = 5,
          clab= metadata_to_plot
)

pdf(paste0("all_cells_163k_umap3d_",metadata_to_plot,".pdf"), width = 11.69, height = 8.27, paper = "special")
pdf <-scatter3D(x = plot.data_meta[,1], y = plot.data_meta[,2], z = plot.data_meta[,3], 
                xlab= "UMAP-1", ylab= "UMAP-2", zlab= "UMAP-3",
                colvar= as.integer(as.factor(plot.data_meta[,4])),
                col = colors1,
                bty = "b2", 
                pch = 20, cex = 0.2, ticktype = "detailed",
                theta = 60, phi = 5,
                clab= metadata_to_plot
)
plot(pdf)
dev.off()
rm(pdf)

#=======================================================================================================
# run the function again setting the resolution according to the above results
DefaultAssay(lung) <- "integrated"

# export Umap plots with the identified clusters, the age and sample name
pdf("all_cells_umap2d_clusters.pdf",width=11.69,height=8.27,paper="special")
all_cells_umap2d_clusters <- DimPlot(lung, reduction = "umap2d", cols =  colors, label = TRUE, pt.size = 0.5)
plot(all_cells_umap2d_clusters)
rm(all_cells_umap2d_clusters)
dev.off()

pdf("all_cells_umap2d_clusters1.pdf",width=11.69,height=8.27,paper="special")
all_cells_umap2d_clusters1 <- DimPlot(lung, reduction = "umap2d", cols =  colors, label = TRUE, pt.size = 0.1)
plot(all_cells_umap2d_clusters1)
rm(all_cells_umap2d_clusters1)
dev.off()

pdf("all_cells_umap2d_age.pdf",width=11.69,height=8.27,paper="special")
all_cells_umap2d_age <- DimPlot(lung, reduction = "umap2d", cols =  "alphabet2", group.by = "age", pt.size = 0.5)
plot(all_cells_umap2d_age)
rm(all_cells_umap2d_age)
dev.off()

pdf("all_cells_umap2d_age1.pdf",width=11.69,height=8.27,paper="special")
all_cells_umap2d_age1 <- DimPlot(lung, reduction = "umap2d", cols =  "alphabet2", group.by = "age", shuffle = TRUE, pt.size = 0.1)
plot(all_cells_umap2d_age1)
rm(all_cells_umap2d_age1)
dev.off()

pdf("all_cells_umap2d_orig_ident.pdf",width=11.69,height=8.27,paper="special")
all_cells_umap2d_orig_ident <- DimPlot(lung, reduction = "umap2d", cols =  "alphabet2", group.by = "orig.ident", pt.size = 0.5)
plot(all_cells_umap2d_orig_ident)
rm(all_cells_umap2d_orig_ident)
dev.off()

pdf("all_cells_umap2d_orig_ident1.pdf",width=11.69,height=8.27,paper="special")
all_cells_umap2d_orig_ident1 <- DimPlot(lung, reduction = "umap2d", cols =  "alphabet2", group.by = "orig.ident", shuffle = TRUE,pt.size = 0.1)
plot(all_cells_umap2d_orig_ident1)
rm(all_cells_umap2d_orig_ident1)
dev.off()

pdf("all_cells_umap2d_SampleName.pdf",width=11.69,height=8.27,paper="special")
all_cells_umap2d_SampleName <- DimPlot(lung, reduction = "umap2d", cols =  "alphabet2", group.by = "SampleName", pt.size = 0.5)
plot(all_cells_umap2d_SampleName)
rm(all_cells_umap2d_SampleName)
dev.off()

pdf("all_cells_umap2d_SampleName1.pdf",width=11.69,height=8.27,paper="special")
all_cells_umap2d_SampleName1 <- DimPlot(lung, reduction = "umap2d", cols =  "alphabet2", group.by = "SampleName", shuffle = TRUE, pt.size = 0.1)
plot(all_cells_umap2d_SampleName1)
rm(all_cells_umap2d_SampleName1)
dev.off()

pdf("all_cells_umap2d_donor.pdf",width=11.69,height=8.27,paper="special")
all_cells_umap2d_donor <- DimPlot(lung, reduction = "umap2d", cols =  "alphabet2", group.by = "donor", pt.size = 0.5)
plot(all_cells_umap2d_donor)
rm(all_cells_umap2d_donor)
dev.off()

pdf("all_cells_umap2d_donor1.pdf",width=11.69,height=8.27,paper="special")
all_cells_umap2d_donor1 <- DimPlot(lung, reduction = "umap2d", cols =  "alphabet2", group.by = "donor",shuffle = TRUE, pt.size = 0.1)
plot(all_cells_umap2d_donor1)
rm(all_cells_umap2d_donor1)
dev.off()

pdf("all_cells_umap2d_chemistry.pdf",width=11.69,height=8.27,paper="special")
all_cells_umap2d_chemistry <- DimPlot(lung, reduction = "umap2d", group.by = "chemistry", pt.size = 0.5)
plot(all_cells_umap2d_chemistry)
rm(all_cells_umap2d_chemistry)
dev.off()

pdf("all_cells_umap2d_chemistry1.pdf",width=11.69,height=8.27,paper="special")
all_cells_umap2d_chemistry1 <- DimPlot(lung, reduction = "umap2d", group.by = "chemistry",shuffle = TRUE, pt.size = 0.1)
plot(all_cells_umap2d_chemistry1)
rm(all_cells_umap2d_chemistry1)
dev.off()

pdf("all_cells_umap2d_features.pdf",width=11.69,height=8.27,paper="special")
all_cells_umap2d_features <- FeaturePlot(lung, reduction = "umap2d",  features = "nFeature_RNA", pt.size = 0.5)
plot(all_cells_umap2d_features)
rm(all_cells_umap2d_features)
dev.off()

pdf("all_cells_umap2d_percent_mt.pdf",width=11.69,height=8.27,paper="special")
all_cells_umap2d_percent_mt <- FeaturePlot(lung, reduction = "umap2d",  features = "percent.MT", pt.size = 0.5)
plot(all_cells_umap2d_percent_mt)
rm(all_cells_umap2d_percent_mt)
dev.off()

pdf("all_cells_tsne_clusters.pdf",width=11.69,height=8.27,paper="special")
s10X140_cells_tsne_clusters <- DimPlot(lung, reduction = "tsne", cols =  colors, label = TRUE, pt.size = 0.5)
plot(s10X140_cells_tsne_clusters)
rm(s10X140_cells_tsne_clusters)
dev.off()

pdf("all_cells_tsne_age.pdf",width=11.69,height=8.27,paper="special")
all_cells_tsne_age <- DimPlot(lung, reduction = "tsne", cols =  "alphabet2", group.by = "age", pt.size = 0.5)
plot(all_cells_tsne_age)
rm(all_cells_tsne_age)
dev.off()

pdf("all_cells_tsne_SampleName.pdf",width=11.69,height=8.27,paper="special")
all_cells_tsne_SampleName <- DimPlot(lung, reduction = "tsne", cols =  "alphabet2", group.by = "SampleName", pt.size = 0.5)
plot(all_cells_tsne_SampleName)
rm(all_cells_tsne_SampleName)
dev.off()

pdf("all_cells_tsne_donor.pdf",width=11.69,height=8.27,paper="special")
all_cells_tsne_donor <- DimPlot(lung, reduction = "tsne", cols =  "alphabet2", group.by = "donor", pt.size = 0.5)
plot(all_cells_tsne_donor)
rm(all_cells_tsne_donor)
dev.off()

pdf("all_cells_tsne_chemistry.pdf",width=11.69,height=8.27,paper="special")
all_cells_tsne_chemistry <- DimPlot(lung, reduction = "tsne", group.by = "chemistry", pt.size = 0.5)
plot(all_cells_tsne_chemistry)
rm(all_cells_tsne_chemistry)
dev.off()

pdf("all_cells_tsne_percent_mt.pdf",width=11.69,height=8.27,paper="special")
all_cells_tsne_percent_mt <- FeaturePlot(lung, reduction = "tsne", features = "percent.MT", pt.size = 0.5)
plot(all_cells_tsne_percent_mt)
rm(all_cells_tsne_percent_mt)
dev.off()

pdf("all_cells_tsne_features.pdf",width=11.69,height=8.27,paper="special")
all_cells_tsne_features <- FeaturePlot(lung, reduction = "tsne", features = "nFeature_RNA", pt.size = 0.5)
plot(all_cells_tsne_features)
rm(all_cells_tsne_features)
dev.off()

#=============================================================================================
# create *.tiff files with individual clusters being colored with red

for (i in 0:(nlevels(lung@active.ident)-1)) {
  cluster <- WhichCells(lung, ident = i)
  name = sprintf("all_cells_umap2d_cluster%s.tiff", i)
  tiff(name, width = 707, height = 500, compression = "lzw")
  image <- DimPlot(lung, reduction = "umap2d",  cells.highlight = cluster, sizes.highlight= 0.5, pt.size = 0.5) + NoLegend()
  plot(image)
  rm(image)
  dev.off()
}

for (i in 0:(nlevels(lung@active.ident)-1)) {
  cluster <- WhichCells(lung, ident = i)
  name = sprintf("all_cells_tsne_cluster%s.tiff", i)
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
  
  # Make a column of row name identities (these will be your cell/barcode naall_cells)
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
                  theta = 60, phi = 5,
                  clab= paste0("Cluster", i)
  )
  plot(pdf)
  dev.off()
  rm(pdf)
}

save.image(file = "all_cells_sct_analysis1.RData")
#=========================================================================

#=========================================================================================
# plot genes in 3D
goi <- "FOXJ1"
DefaultAssay(lung) <- "RNA"
plotting.data <- FetchData(object = lung, vars = c("UMAP_1", "UMAP_2", "UMAP_3", "Expression"=goi), slot = 'data')

# Add the label column, so that now the column has 'cellname-its expression value'
plotting.data$label <- paste(rownames(plotting.data)," - ", plotting.data[,goi], sep="")

# Plot your data, in this example my Seurat object had 21 clusters (0-20), and cells express a gene called ACTB
plot_ly(data = plotting.data,
        # name = goi,
        x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, 
        color = plotting.data[,4], #~TOP2A,
        opacity = .5,
        colors = c('darkgrey', 'red'), 
        type = "scatter3d", 
        mode = "markers",
        marker = list(size = 1), 
        text=~label,
        hoverinfo="text"
) %>%layout(title=goi)


# plot the 3D UMAP using plot3D package, theta and phi change the orinetation of the graph 
scatter3D(x = plotting.data[,1], y = plotting.data[,2], z = plotting.data[,3], 
          xlab= "UMAP-1", ylab= "UMAP-2", zlab= "UMAP-3",
          colvar= plotting.data[,4],
          col = ramp.col(c('darkgrey', 'red')),
          bty = "b2", 
          pch = 20, cex = 0.5, ticktype = "detailed",
          theta = 15, phi = 20,
          clab = c("FOXJ1")
)


#=====================================================================================================================================
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
min_cluster = 126
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
write.csv(lung_cluster.markers, file="all_cells_sct_lung_markers.csv")


# filter out insignificant genes
lung_cluster.markers <- lung_cluster.markers[lung_cluster.markers$p_val_adj<0.001, ]

# obtain the top 25 genes of each cluster, based on the average logFC
top25 <- lung_cluster.markers %>% group_by(cluster) %>% top_n(n = 25, wt = avg_log2FC)
write.csv(top25, file="all_cells_sct_lung_cluster_top25_markers.csv")


top5 <- lung_cluster.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
write.csv(top5, file="all_cells_sct_lung_cluster_top5_markers.csv")


#-------------------------------------------------------------------------------------------
# The transcription factor and co-factor lists were downloaded from the http://bioinfo.life.hust.edu.cn/AnimalTFDB/#!/download
# The kinase list "protein_class_kinases.tsv", cd marker list "protein_class_CD.tsv" and "protein_class_secreted.tsv" were downloaded from https://www.proteinatlas.org/


# import the list of transcription factors 
TFs <- read.csv("TFs.csv")

# create a new file by merging the TF-list and the enriched genes in each cluster, to identify the TFs, only
lung_cluster_tfs <- join(lung_cluster.markers, TFs, by="gene", type="inner")
top10_cluster_tfs <- lung_cluster_tfs %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
# export the list of the top 10 TFs based on average expression difference
write.csv(top10_cluster_tfs, file="all_cells_sct__top10_cluster_tfs.csv")

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
write.csv(top10_cluster_cotfs, file="all_cells_sct__top10_cluster_cotfs.csv")


# import the list of kinases
kinases <- data.frame(read.delim("protein_class_kinases.tsv"))
kinases <- kinases[1:10]

# create a new file by merging the kinase-list and the enriched genes in each cluster, to identify the kinases, only
lung_cluster_kinases <- join(lung_cluster.markers, kinases, by="gene", type="inner")
top10_cluster_kinases <- lung_cluster_kinases %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
# export the list of the top 10 kinases based on average expression difference
write.csv(top10_cluster_kinases, file="all_cells_sct__top10_cluster_kinases.csv")


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
write.csv(lung_roc_cluster.markers, file="all_cells_163k_roc_markers.csv")

lung_roc_cluster.markers <- lung_roc_cluster.markers[lung_roc_cluster.markers$Dpct >0.1, ]

# obtain the top 25 genes of each cluster, based on the power
top25_roc <- lung_roc_cluster.markers %>% group_by(cluster) %>% top_n(n = 25, wt = power)
write.csv(top25_roc, file="all_cells_163k_roc_cluster_top25_markers.csv")


top5_roc <- lung_roc_cluster.markers %>% group_by(cluster) %>% top_n(n = 5, wt = power)
write.csv(top5_roc, file="all_cells_163k_roc_cluster_top5_markers.csv")


top10_roc <- lung_roc_cluster.markers %>% group_by(cluster) %>% top_n(n = 10, wt = power)
write.csv(top10_roc, file="all_cells_163k_roc_cluster_top10_markers.csv")

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
write.csv(top10_roc_cluster_tfs, file="all_cells_163k_roc_top10_cluster_tfs.csv")


# import the list of transcription cofactors 
coTFs <- read.csv("TF_cofactors.csv")

# create a new file by merging the coTF-list and the enriched genes in each cluster, to identify the coTFs, only
lung_roc_cluster_cotfs <- join(lung_roc_cluster.markers, coTFs, by="gene", type="inner")
top10_roc_cluster_cotfs <- lung_roc_cluster_cotfs %>% group_by(cluster) %>% top_n(n = 10, wt = power)
# export the list of the top 10 coTFs based on average expression difference
write.csv(top10_roc_cluster_cotfs, file="all_cells_163k_roc_top10_cluster_cotfs.csv")


# import the list of kinases
kinases <- data.frame(read.delim("protein_class_kinases.tsv"))
kinases <- kinases[1:10]

# create a new file by merging the kinase-list and the enriched genes in each cluster, to identify the kinases, only
lung_roc_cluster_kinases <- join(lung_roc_cluster.markers, kinases, by="gene", type="inner")
top10_roc_cluster_kinases <- lung_roc_cluster_kinases %>% group_by(cluster) %>% top_n(n = 10, wt = power)
# export the list of the top 10 kinases based on average expression difference
write.csv(top10_roc_cluster_kinases, file="all_cells_163k_roc_top10_cluster_kinases.csv")


# import the list of CD markers
cds <- data.frame(read.delim("protein_class_CD.tsv"))
cds <- cds[1:10]

# create a new file by merging the CD-list and the enriched genes in each cluster, to identify the surface markers, only
lung_roc_cluster_cds <- join(lung_roc_cluster.markers, cds, by="gene", type="inner")
top10_roc_cluster_cds <- lung_roc_cluster_cds %>% group_by(cluster) %>% top_n(n = 10, wt = power)
# export the list of the top 10 CDs based on average expression difference
write.csv(top10_roc_cluster_cds, file="top10_roc_cluster_cds.csv")

#============================================
# create and export a table with the cell age and their corresponding cluster 
all_cells_sct_cells_clusters <- data.frame(lung@active.ident)
write.csv(all_cells_sct_cells_clusters, file="all_cells_sct_cells_clusters_sct.csv")

Age <- data.frame(lung@meta.data[["age"]])
experiment <- data.frame(lung@meta.data[["orig.ident"]])
cell_name <- rownames(all_cells_sct_cells_clusters)
cluster_list <- cbind(cell_name, all_cells_sct_cells_clusters)
cluster_list <- cbind(cluster_list, Age)
cluster_list <- cbind(cluster_list, experiment)
cluster_list <- setNames(cluster_list, c("cell", "cluster", "age", "orig_ident"))
write.csv(cluster_list, file="all_cells_sct_metadata.csv")

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
#============================================
donor <- data.frame(lung@meta.data[["donor"]])
experiment <- data.frame(lung@meta.data[["age"]])
cell_name <- rownames(all_cells_sct_cells_clusters)
cluster_list <- cbind(cell_name, all_cells_sct_cells_clusters)
cluster_list <- cbind(cluster_list, donor)
cluster_list <- cbind(cluster_list, experiment)
cluster_list <- setNames(cluster_list, c("cell", "cluster", "donor", "age"))
write.csv(cluster_list, file="all_cells_sct_metadata_donor.csv")

# summarize the abundance of cells from specific time points in clusters 
cluster_list$cluster <- as.factor(cluster_list$cluster)
cluster_list$donor <- as.factor(cluster_list$donor)
abundance_matrix <- matrix(nrow = nlevels(cluster_list$cluster), ncol = nlevels(cluster_list$donor))
rownames(abundance_matrix) <- levels(cluster_list$cluster)
colnames(abundance_matrix) <- levels(cluster_list$donor)
for (i in 0:nlevels(cluster_list$cluster)-1) 
  for (j in 1:nlevels(cluster_list$donor)) {
    list1 <-cluster_list[cluster_list$cluster == i,]
    abundance_matrix[(i+1),j] <- sum(list1$donor==levels(list1$donor)[j])
  }
write.csv(abundance_matrix, file="donor_abundance_in_clusters.csv")


# calculate the % of cells of each time-point in the clusters
donor_size <- colSums(abundance_matrix)
abundance_matrix_percentdonors <- (abundance_matrix*100)/donor_size
write.csv(abundance_matrix_percentdonors, file="cluster_percent_in_donors.csv")

save.image(file = "all_cells_analysis_data_after_statistics.RData")
#=============================================================================================
# create subsets of the dataset for individual cluster analyses
epi <- subset(lung, idents = c("12", "13", "24", "27"))
save(epi, file = "epi.RData")
rm(epi)

endo <- subset(lung, idents = c("10", "26"))
save(endo, file = "endo.RData")
rm(endo)

neuro_ne <- subset(lung, idents = c("23", "24"))
save(neuro_ne, file = "neuronal_ne.RData")
rm(neuro_ne)

neuro <- subset(lung, idents = c("23"))
save(neuro, file = "neuronal.RData")
rm(neuro)

immune <- subset(lung, idents = c("17", "22", "28"))
save(immune, file = "immune.RData")
rm(immune)

immune_rbc <- subset(lung, idents = c("17", "21", "22", "28"))
save(immune_rbc, file = "immune_rbc.RData")
rm(immune_rbc)

endo_pericytes <- subset(lung, idents = c("10", "16", "26"))
save(endo_pericytes, file = "endo_pericytes.RData")
rm(endo_pericytes)
#==============================================================================================
#=============================================================================================
