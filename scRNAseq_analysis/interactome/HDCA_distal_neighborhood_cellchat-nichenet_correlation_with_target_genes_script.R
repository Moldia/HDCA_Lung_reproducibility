setwd()

wd <- "C:/Users/alex/Desktop/HDCA_lung_project/cellchat/cellchat_distal03_dif/"

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
library(harmony)
library(pheatmap)
library(RColorBrewer)
library(sctransform)
library(SeuratWrappers)
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
library(nichenetr)
library(circlize)
library(CellChat)
library(patchwork)
library(NMF)
library(ggalluvial)
options(stringsAsFactors = FALSE)

# download the ligand-target matrix of nichenet
#ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
# or load it, if you have stored it locally
load("ligand_target_matrix_nichenet.RData")

# obtain the significant ligand-receptor pairs from CellChat analysis
sign_LR <- cellchat@LR[["LRsig"]]

#===================================================================
# create a subset seurat object with only the receiver clusters
receivers_senders <- subset(lung, ident=c("epi_cl2", "epi_cl10", "epi_cl3", "epi_cl9", "epi_cl8", "epi_cl13", "epi_cl5", 
                                          "mes_cl0", "mes_cl2","mes_cl6", "mes_cl8", "mes_cl20", "mes_cl12"))
# select the ligands that are found in our dataset
sign_LR <- sign_LR[sign_LR$ligand %in% row.names(receivers_senders),]

data.input <- data.input[row.names(data.input) %in% row.names(max_value),]

# keep only the ligand-target pairs with corresponding ligands in the ligand-receptor pairs
filtered_LT <- data.frame(ligand_target_matrix[,colnames(ligand_target_matrix) %in% sign_LR$ligand])

# keep only the target genes that are found in their mean log2 expression is 0.5, in at least 1 cluster of the whole dataset. In that way, we remove genes with negligible expression 
filtered_LT <- data.frame(filtered_LT[rownames(filtered_LT) %in% max_value$gene,])

# to avoid ubiquitously expressed genes we run differential expression analysis and select only the target genes being in the differentially expressed genes 
min_cluster <- min(table(lung@active.ident))
table(lung@active.ident)
#min_cluster <- 250
lung_cluster.markers <- FindAllMarkers(lung, test.use = "MAST", min.diff.pct = 0.1, max.cells.per.ident = min_cluster, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.1)
# calculate the percentage difference of the positive cells in the indicated cluster, in comparison to the others
lung_cluster.markers$Dpct <- (lung_cluster.markers$pct.1 - lung_cluster.markers$pct.2)

filtered_LT <- data.frame(filtered_LT[rownames(filtered_LT) %in% lung_cluster.markers$gene,])

# export pdf files for all the ligands with their top10 target genes, according to Nichenet
for (i in 1:ncol(filtered_LT)){
  a<- data.frame(filtered_LT[,i])
  row.names(a)<- row.names(filtered_LT)
  names(a) <- colnames(filtered_LT)[i]
  a$target <- row.names(a)
  a <- a[order(a[,1],decreasing=TRUE),]
  
  #take the 1st 10 genes
  a <- top_n(a, 10, a[,1])
  a<- a[1:10,]
  
  #remove duplicates
 f<-c(colnames(a)[1], rownames(a))
 f <- f[!duplicated(f)]
 
  pdf(paste0(wd,"ligand_targets_strict/",colnames(a)[1], "_top10_target_genes_dotplot.pdf"),width=8,height=4,paper="special")
  b <- DotPlot(receivers_senders, assay = "RNA", features = f, scale = TRUE) + RotatedAxis()
  plot(b)
  rm(b)
  dev.off()
}

# export pdf files for all the ligands with their top10 target genes, according to Nichenet
for (i in 1:ncol(filtered_LT)){
  a<- data.frame(filtered_LT[,i])
  row.names(a)<- row.names(filtered_LT)
  names(a) <- colnames(filtered_LT)[i]
  a$target <- row.names(a)
  a <- a[order(a[,1],decreasing=TRUE),]
  
  #take the 1st 10 genes
  a <- top_n(a, 20, a[,1])
  a<- a[1:20,]
  
  #remove duplicates
  f<-c(colnames(a)[1], rownames(a))
  f <- f[!duplicated(f)]
  
  pdf(paste0(wd,"ligand_targets_strict/",colnames(a)[1], "_top20_target_genes_dotplot.pdf"),width=8,height=4,paper="special")
  b <- DotPlot(receivers_senders, assay = "RNA", features = f, scale = TRUE) + RotatedAxis()
  plot(b)
  rm(b)
  dev.off()
}

# export pdf files for all the ligands with their top10 target genes, according to Nichenet
for (i in 1:ncol(filtered_LT)){
  a<- data.frame(filtered_LT[,i])
  row.names(a)<- row.names(filtered_LT)
  names(a) <- colnames(filtered_LT)[i]
  a$target <- row.names(a)
  a <- a[order(a[,1],decreasing=TRUE),]
  
  #take the 1st 10 genes
  a <- top_n(a, 25, a[,1])
  a<- a[1:25,]
  
  #remove duplicates
  f<-c(colnames(a)[1], rownames(a))
  f <- f[!duplicated(f)]
  
  pdf(paste0(wd,"ligand_targets_strict/",colnames(a)[1], "_top25_target_genes_dotplot.pdf"),width=9,height=4,paper="special")
  b <- DotPlot(receivers_senders, assay = "RNA", features = f, scale = TRUE) + RotatedAxis()
  plot(b)
  rm(b)
  dev.off()
}

# export pdf files for all the ligands with their top10 target genes, according to Nichenet
for (i in 1:ncol(filtered_LT)){
  a<- data.frame(filtered_LT[,i])
  row.names(a)<- row.names(filtered_LT)
  names(a) <- colnames(filtered_LT)[i]
  a$target <- row.names(a)
  a <- a[order(a[,1],decreasing=TRUE),]
  
  #take the 1st 10 genes
  a <- top_n(a, 50, a[,1])
  a<- a[1:50,]
  
  #remove duplicates
  f<-c(colnames(a)[1], rownames(a))
  f <- f[!duplicated(f)]
  
  pdf(paste0(wd,"ligand_targets_strict/",colnames(a)[1], "_top50_target_genes_dotplot.pdf"),width=20,height=4,paper="special")
  b <- DotPlot(receivers_senders, assay = "RNA", features = f, scale = TRUE) + RotatedAxis()
  plot(b)
  rm(b)
  dev.off()
}
save.image("final_analysis_data.RData")
#==============================================
rownames <- row.names(filtered_LT)
filtered_LT1 <- data.frame(filtered_LT[,colnames(filtered_LT)=="FGF9"])
row.names(filtered_LT1) <- rownames
names(filtered_LT1) <- "FGF9"


#take the 1st 10 genes
filtered_LT1$target <-row.names(filtered_LT1) 
filtered_LT1 <- filtered_LT1[order(filtered_LT1[,1],decreasing=TRUE),]
a <- top_n(filtered_LT1, 20, filtered_LT1[,1])
a <- row.names(a)
a<- a[1:20]

#remove duplicates
f<-c("FGF9", "FGFR2", a)
f <- f[!duplicated(f)]
write.csv(f, file = "Ext_Data_Fig9H_FGF9_target_genes.csv")
pdf("FGF9_LR_top20_target_genes_dotplot_strict.pdf",width=8,height=4,paper="special")
b <- DotPlot(receivers_senders, assay = "RNA", features = f, scale = TRUE) + RotatedAxis()
plot(b)
rm(b)
dev.off()

#==============================================
rownames <- row.names(filtered_LT)
filtered_LT1 <- data.frame(filtered_LT[,colnames(filtered_LT)=="WNT7B"])
row.names(filtered_LT1) <- rownames
names(filtered_LT1) <- "WNT7B"


#take the 1st 10 genes
filtered_LT1$target <-row.names(filtered_LT1) 
filtered_LT1 <- filtered_LT1[order(filtered_LT1[,1],decreasing=TRUE),]
a <- top_n(filtered_LT1, 20, filtered_LT1[,1])
a <- row.names(a)
a<- a[1:20]

#remove duplicates
f<-c("WNT7B", "FZD1", "FZD2", "FZD7", "LRP5", "LRP6", a)
f <- f[!duplicated(f)]

pdf("WNT7B_LR_top20_target_genes_dotplot_strict.pdf",width=9,height=4,paper="special")
b <- DotPlot(receivers_senders, assay = "RNA", features = f, scale = TRUE) + RotatedAxis()
plot(b)
rm(b)
dev.off()

pdf("WNT7B_LR_top20_target_genes_dotplot_strict_unscalled.pdf",width=9,height=4,paper="special")
b <- DotPlot(receivers_senders, assay = "RNA", features = f, scale = F) + RotatedAxis()
plot(b)
rm(b)
dev.off()

