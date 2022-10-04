setwd()
wd <-"C:/Users/alex/Desktop/cellchat/cellchat_proximal03_dif/"


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
receivers_senders <- subset(lung, ident=c("epi_cl1", "epi_cl4", "epi_cl0", "epi_cl6", "epi_cl7", "epi_cl11", "epi_cl12", "epi_cl14", "neuronal", "mes_cl13", "mes_cl5", "mes_cl16", "mes_cl18", "mes_cl10"))


levels(receivers_senders) <- c("epi_cl14", "epi_cl11", "epi_cl12", "epi_cl7","epi_cl6", "epi_cl0", "epi_cl4", "epi_cl1","neuronal", "mes_cl10", "mes_cl13", "mes_cl16", "mes_cl5", "mes_cl18")

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

#==============================================
rownames <- row.names(filtered_LT)
filtered_LT1 <- data.frame(filtered_LT[,colnames(filtered_LT)=="WNT2"])
row.names(filtered_LT1) <- rownames
names(filtered_LT1) <- "WNT2"

#take the 1st 10 genes
filtered_LT1$target <-row.names(filtered_LT1) 
filtered_LT1 <- filtered_LT1[order(filtered_LT1[,1],decreasing=TRUE),]
a <- top_n(filtered_LT1, 20, filtered_LT1[,1])
a <- row.names(a)
a<- a[1:20]

#remove duplicates
f<-c("WNT2", "FZD1", "FZD2", "FZD7", "LRP6", a)
f <- f[!duplicated(f)]

pdf("WNT2_LR_top20_target_genes_dotplot_strict.pdf",width=8,height=4,paper="special")
b <- DotPlot(receivers_senders, assay = "RNA", features = f, scale = TRUE) + RotatedAxis()
plot(b)
rm(b)
dev.off()

rownames <- row.names(filtered_LT)
filtered_LT1 <- data.frame(filtered_LT[,colnames(filtered_LT)=="WNT5A"])
row.names(filtered_LT1) <- rownames
names(filtered_LT1) <- "WNT5A"

#take the 1st 10 genes
filtered_LT1$target <-row.names(filtered_LT1) 
filtered_LT1 <- filtered_LT1[order(filtered_LT1[,1],decreasing=TRUE),]
a <- top_n(filtered_LT1, 20, filtered_LT1[,1])
a <- row.names(a)
a<- a[1:20]

#remove duplicates
f<-c("WNT5A", "FZD1", "FZD2", "FZD7", "LRP6", "MCAM", a)
f <- f[!duplicated(f)]

pdf("WNT5A_LR_top20_target_genes_dotplot_strict.pdf",width=8,height=4,paper="special")
b <- DotPlot(receivers_senders, assay = "RNA", features = f, scale = TRUE) + RotatedAxis()
plot(b)
rm(b)
dev.off()


rownames <- row.names(filtered_LT)
filtered_LT1 <- data.frame(filtered_LT[,colnames(filtered_LT)=="WNT11"])
row.names(filtered_LT1) <- rownames
names(filtered_LT1) <- "WNT11"

#take the 1st 10 genes
filtered_LT1$target <-row.names(filtered_LT1) 
filtered_LT1 <- filtered_LT1[order(filtered_LT1[,1],decreasing=TRUE),]
a <- top_n(filtered_LT1, 20, filtered_LT1[,1])
a <- row.names(a)
a<- a[1:20]

#remove duplicates
f<-c("WNT11", "FZD1", "FZD2", "FZD7", "LRP6", a)
f <- f[!duplicated(f)]

pdf("WNT11_LR_top20_target_genes_dotplot_strict.pdf",width=8,height=4,paper="special")
b <- DotPlot(receivers_senders, assay = "RNA", features = f, scale = TRUE) + RotatedAxis()
plot(b)
rm(b)
dev.off()

#--------------------------------------------------------------------------
rownames <- row.names(filtered_LT)
filtered_LT1 <- data.frame(filtered_LT[,colnames(filtered_LT)=="BMP4"])
row.names(filtered_LT1) <- rownames
names(filtered_LT1) <- "BMP4"

#take the 1st 10 genes
filtered_LT1$target <-row.names(filtered_LT1) 
filtered_LT1 <- filtered_LT1[order(filtered_LT1[,1],decreasing=TRUE),]
a <- top_n(filtered_LT1, 20, filtered_LT1[,1])
a <- row.names(a)
a<- a[1:20]

#remove duplicates

f<-c("BMP4", "BMPR1A","BMPR1A", "ACVR1", "BMPR2", "ACVR2A", a)
f <- f[!duplicated(f)]

pdf("BMP4_LR_top20_target_genes_dotplot_strict.pdf",width=8,height=4,paper="special")
b <- DotPlot(receivers_senders, assay = "RNA", features = f, scale = TRUE) + RotatedAxis()
plot(b)
rm(b)
dev.off()

rownames <- row.names(filtered_LT)
filtered_LT1 <- data.frame(filtered_LT[,colnames(filtered_LT)=="BMP5"])
row.names(filtered_LT1) <- rownames
names(filtered_LT1) <- "BMP5"

#take the 1st 10 genes
filtered_LT1$target <-row.names(filtered_LT1) 
filtered_LT1 <- filtered_LT1[order(filtered_LT1[,1],decreasing=TRUE),]
a <- top_n(filtered_LT1, 20, filtered_LT1[,1])
a <- row.names(a)
a<- a[1:20]

#remove duplicates

f<-c("BMP5", "BMPR1A","BMPR1A", "ACVR1", "BMPR2", "ACVR2A", a)
f <- f[!duplicated(f)]

pdf("BMP5_LR_top20_target_genes_dotplot_strict.pdf",width=8,height=4,paper="special")
b <- DotPlot(receivers_senders, assay = "RNA", features = f, scale = TRUE) + RotatedAxis()
plot(b)
rm(b)
dev.off()

#============================================================
rownames <- row.names(filtered_LT)
filtered_LT1 <- data.frame(filtered_LT[,colnames(filtered_LT)=="DLL3"])
row.names(filtered_LT1) <- rownames
names(filtered_LT1) <- "DLL3"

#take the 1st 10 genes
filtered_LT1$target <-row.names(filtered_LT1) 
filtered_LT1 <- filtered_LT1[order(filtered_LT1[,1],decreasing=TRUE),]
a <- top_n(filtered_LT1, 20, filtered_LT1[,1])
a <- row.names(a)
a<- a[1:20]

#remove duplicates

f<-c("DLL3", "NOTCH2", a)
f <- f[!duplicated(f)]

pdf("DLL3_LR_top20_target_genes_dotplot_strict.pdf",width=8,height=4,paper="special")
b <- DotPlot(receivers_senders, assay = "RNA", features = f, scale = TRUE) + RotatedAxis()
plot(b)
rm(b)
dev.off()

rownames <- row.names(filtered_LT)
filtered_LT1 <- data.frame(filtered_LT[,colnames(filtered_LT)=="JAG1"])
row.names(filtered_LT1) <- rownames
names(filtered_LT1) <- "JAG1"

#take the 1st 10 genes
filtered_LT1$target <-row.names(filtered_LT1) 
filtered_LT1 <- filtered_LT1[order(filtered_LT1[,1],decreasing=TRUE),]
a <- top_n(filtered_LT1, 20, filtered_LT1[,1])
a <- row.names(a)
a<- a[1:20]

#remove duplicates

f<-c("JAG1", "NOTCH2", a)
f <- f[!duplicated(f)]

pdf("JAG1_LR_top20_target_genes_dotplot_strict.pdf",width=8,height=4,paper="special")
b <- DotPlot(receivers_senders, assay = "RNA", features = f, scale = TRUE) + RotatedAxis()
plot(b)
rm(b)
dev.off()

#============================================================
rownames <- row.names(filtered_LT)
filtered_LT1 <- data.frame(filtered_LT[,colnames(filtered_LT)=="FGF7"])
row.names(filtered_LT1) <- rownames
names(filtered_LT1) <- "FGF7"

#take the 1st 20 genes
filtered_LT1$target <-row.names(filtered_LT1) 
filtered_LT1 <- filtered_LT1[order(filtered_LT1[,1],decreasing=TRUE),]
a <- top_n(filtered_LT1, 20, filtered_LT1[,1])
a <- row.names(a)
a<- a[1:20]

#remove duplicates

f<-c("FGF7", "FGFR1", "FGFR2", a)
f <- f[!duplicated(f)]

pdf("FGF7_LR_top20_target_genes_dotplot_strict.pdf",width=8,height=4,paper="special")
b <- DotPlot(receivers_senders, assay = "RNA", features = f, scale = TRUE) + RotatedAxis()
plot(b)
rm(b)
dev.off()

#============================================================
rownames <- row.names(filtered_LT)
filtered_LT1 <- data.frame(filtered_LT[,colnames(filtered_LT)=="FGF18"])
row.names(filtered_LT1) <- rownames
names(filtered_LT1) <- "FGF18"

#take the 1st 20 genes
filtered_LT1$target <-row.names(filtered_LT1) 
filtered_LT1 <- filtered_LT1[order(filtered_LT1[,1],decreasing=TRUE),]
a <- top_n(filtered_LT1, 20, filtered_LT1[,1])
a <- row.names(a)
a<- a[1:20]

#remove duplicates

f<-c("FGF18", "FGFR1", "FGFR2", a)
f <- f[!duplicated(f)]

pdf("FGF18_LR_top20_target_genes_dotplot_strict.pdf",width=8,height=4,paper="special")
b <- DotPlot(receivers_senders, assay = "RNA", features = f, scale = TRUE) + RotatedAxis()
plot(b)
rm(b)
dev.off()

#============================================================
rownames <- row.names(filtered_LT)
filtered_LT1 <- data.frame(filtered_LT[,colnames(filtered_LT)=="NRXN1"])
row.names(filtered_LT1) <- rownames
names(filtered_LT1) <- "NRXN1"

#take the 1st 20 genes
filtered_LT1$target <-row.names(filtered_LT1) 
filtered_LT1 <- filtered_LT1[order(filtered_LT1[,1],decreasing=TRUE),]
a <- top_n(filtered_LT1, 20, filtered_LT1[,1])
a <- row.names(a)
a<- a[1:20]

#remove duplicates

f<-c("NRXN1", "NLGN1", a)
f <- f[!duplicated(f)]

pdf("NRXN1_LR_top20_target_genes_dotplot_strict.pdf",width=8,height=4,paper="special")
b <- DotPlot(receivers_senders, assay = "RNA", features = f, scale = TRUE) + RotatedAxis()
plot(b)
rm(b)
dev.off()

pdf("NRXN1_LR_top20_target_genes_dotplot_strict_not_scaled.pdf",width=8,height=4,paper="special")
b <- DotPlot(receivers_senders, assay = "RNA", features = f, scale = F) + RotatedAxis()
plot(b)
rm(b)
dev.off()
#============================================================
rownames <- row.names(filtered_LT)
filtered_LT1 <- data.frame(filtered_LT[,colnames(filtered_LT)=="NRXN3"])
row.names(filtered_LT1) <- rownames
names(filtered_LT1) <- "NRXN3"

#take the 1st 20 genes
filtered_LT1$target <-row.names(filtered_LT1) 
filtered_LT1 <- filtered_LT1[order(filtered_LT1[,1],decreasing=TRUE),]
a <- top_n(filtered_LT1, 20, filtered_LT1[,1])
a <- row.names(a)
a<- a[1:20]

#remove duplicates

f<-c("NRXN3", "NLGN1", a)
f <- f[!duplicated(f)]

pdf("NRXN3_LR_top20_target_genes_dotplot_strict.pdf",width=8,height=4,paper="special")
b <- DotPlot(receivers_senders, assay = "RNA", features = f, scale = TRUE) + RotatedAxis()
plot(b)
rm(b)
dev.off()

pdf("NRXN3_LR_top20_target_genes_dotplot_strict_not_scaled.pdf",width=8,height=4,paper="special")
b <- DotPlot(receivers_senders, assay = "RNA", features = f, scale = F) + RotatedAxis()
plot(b)
rm(b)
dev.off()

#============================================================
rownames <- row.names(filtered_LT)
filtered_LT1 <- data.frame(filtered_LT[,colnames(filtered_LT)=="HGF"])
row.names(filtered_LT1) <- rownames
names(filtered_LT1) <- "HGF"

#take the 1st 20 genes
filtered_LT1$target <-row.names(filtered_LT1) 
filtered_LT1 <- filtered_LT1[order(filtered_LT1[,1],decreasing=TRUE),]
a <- top_n(filtered_LT1, 20, filtered_LT1[,1])
a <- row.names(a)
a<- a[1:20]

#remove duplicates

f<-c("HGF", "MET", a)
f <- f[!duplicated(f)]

pdf("HGF_LR_top20_target_genes_dotplot_strict.pdf",width=8,height=4,paper="special")
b <- DotPlot(receivers_senders, assay = "RNA", features = f, scale = TRUE) + RotatedAxis()
plot(b)
rm(b)
dev.off()

pdf("HGF_LR_top20_target_genes_dotplot_strict_not_scaled.pdf",width=8,height=4,paper="special")
b <- DotPlot(receivers_senders, assay = "RNA", features = f, scale = F) + RotatedAxis()
plot(b)
rm(b)
dev.off()

#============================================================
rownames <- row.names(filtered_LT)
filtered_LT1 <- data.frame(filtered_LT[,colnames(filtered_LT)=="SST"])
row.names(filtered_LT1) <- rownames
names(filtered_LT1) <- "SST"

#take the 1st 20 genes
filtered_LT1$target <-row.names(filtered_LT1) 
filtered_LT1 <- filtered_LT1[order(filtered_LT1[,1],decreasing=TRUE),]
a <- top_n(filtered_LT1, 20, filtered_LT1[,1])
a <- row.names(a)
a<- a[1:20]

#remove duplicates

f<-c("SST", "SSTR2", a)
f <- f[!duplicated(f)]

pdf("SST_LR_top20_target_genes_dotplot_strict.pdf",width=8,height=4,paper="special")
b <- DotPlot(receivers_senders, assay = "RNA", features = f, scale = TRUE) + RotatedAxis()
plot(b)
rm(b)
dev.off()

pdf("SST_LR_top20_target_genes_dotplot_strict_not_scaled.pdf",width=8,height=4,paper="special")
b <- DotPlot(receivers_senders, assay = "RNA", features = f, scale = F) + RotatedAxis()
plot(b)
rm(b)
dev.off()
#============================================================
rownames <- row.names(filtered_LT)
filtered_LT1 <- data.frame(filtered_LT[,colnames(filtered_LT)=="SST"])
row.names(filtered_LT1) <- rownames
names(filtered_LT1) <- "SST"

#take the 1st 100 genes
filtered_LT1$target <-row.names(filtered_LT1) 
filtered_LT1 <- filtered_LT1[order(filtered_LT1[,1],decreasing=TRUE),]
a <- top_n(filtered_LT1, 100, filtered_LT1[,1])
a <- row.names(a)
a<- a[1:100]

#remove duplicates

f<-c("SST", "SSTR2", a)
f <- f[!duplicated(f)]

pdf("SST_LR_top100_target_genes_dotplot_strict.pdf",width=25,height=4,paper="special")
b <- DotPlot(receivers_senders, assay = "RNA", features = f, scale = TRUE) + RotatedAxis()
plot(b)
rm(b)
dev.off()

pdf("SST_LR_top100_target_genes_dotplot_strict_not_scaled.pdf",width=25,height=4,paper="special")
b <- DotPlot(receivers_senders, assay = "RNA", features = f, scale = F) + RotatedAxis()
plot(b)
rm(b)
dev.off()

#============================
f<-c("SST", "SSTR1","SSTR2","SSTR3","SSTR4", a)
f <- f[!duplicated(f)]

pdf("SST_LR_top100_target_genes_dotplot_strict1.pdf",width=25,height=4,paper="special")
b <- DotPlot(receivers_senders, assay = "RNA", features = f, scale = TRUE) + RotatedAxis()
plot(b)
rm(b)
dev.off()

pdf("SST_LR_top100_target_genes_dotplot_strict_not_scaled1.pdf",width=25,height=4,paper="special")
b <- DotPlot(receivers_senders, assay = "RNA", features = f, scale = F) + RotatedAxis()
plot(b)
rm(b)
dev.off()


#==============================================
rownames <- row.names(filtered_LT)
filtered_LT1 <- data.frame(filtered_LT[,colnames(filtered_LT)=="THBS1"])
row.names(filtered_LT1) <- rownames
names(filtered_LT1) <- "THBS1"


#take the 1st 10 genes
filtered_LT1$target <-row.names(filtered_LT1) 
filtered_LT1 <- filtered_LT1[order(filtered_LT1[,1],decreasing=TRUE),]
a <- top_n(filtered_LT1, 20, filtered_LT1[,1])
a <- row.names(a)
a<- a[1:20]

#remove duplicates
f<-c("THBS1", "SDC4","CD47", a)
f <- f[!duplicated(f)]

pdf("THBS1_LR_top20_target_genes_dotplot_strict.pdf",width=9,height=4,paper="special")
b <- DotPlot(receivers_senders, assay = "RNA", features = f, scale = TRUE) + RotatedAxis()
plot(b)
rm(b)
dev.off()

pdf("THBS1_LR_top20_target_genes_dotplot_strict_unscalled.pdf",width=9,height=4,paper="special")
b <- DotPlot(receivers_senders, assay = "RNA", features = f, scale = F) + RotatedAxis()
plot(b)
rm(b)
dev.off()

#==============================================
rownames <- row.names(filtered_LT)
filtered_LT1 <- data.frame(filtered_LT[,colnames(filtered_LT)=="HBEGF"])
row.names(filtered_LT1) <- rownames
names(filtered_LT1) <- "HBEGF"


#take the 1st 10 genes
filtered_LT1$target <-row.names(filtered_LT1) 
filtered_LT1 <- filtered_LT1[order(filtered_LT1[,1],decreasing=TRUE),]
a <- top_n(filtered_LT1, 20, filtered_LT1[,1])
a <- row.names(a)
a<- a[1:20]

#remove duplicates
f<-c("HBEGF", "EGFR", a)
f <- f[!duplicated(f)]

pdf("HBEGF_LR_top20_target_genes_dotplot_strict.pdf",width=9,height=4,paper="special")
b <- DotPlot(receivers_senders, assay = "RNA", features = f, scale = TRUE) + RotatedAxis()
plot(b)
rm(b)
dev.off()

pdf("HBEGF_LR_top20_target_genes_dotplot_strict_unscalled.pdf",width=9,height=4,paper="special")
b <- DotPlot(receivers_senders, assay = "RNA", features = f, scale = F) + RotatedAxis()
plot(b)
rm(b)
dev.off()

#==============================================
rownames <- row.names(filtered_LT)
filtered_LT1 <- data.frame(filtered_LT[,colnames(filtered_LT)=="HBEGF"])
row.names(filtered_LT1) <- rownames
names(filtered_LT1) <- "HBEGF"


#take the 1st 10 genes
filtered_LT1$target <-row.names(filtered_LT1) 
filtered_LT1 <- filtered_LT1[order(filtered_LT1[,1],decreasing=TRUE),]
a <- top_n(filtered_LT1, 20, filtered_LT1[,1])
a <- row.names(a)
a<- a[1:20]

#remove duplicates
f<-c("HBEGF", "EGFR", "ERBB4", a)
f <- f[!duplicated(f)]

pdf("HBEGF_LR_top20_target_genes_dotplot_strict1.pdf",width=9,height=4,paper="special")
b <- DotPlot(receivers_senders, assay = "RNA", features = f, scale = TRUE) + RotatedAxis()
plot(b)
rm(b)
dev.off()

pdf("HBEGF_LR_top20_target_genes_dotplot_strict_unscalled1.pdf",width=9,height=4,paper="special")
b <- DotPlot(receivers_senders, assay = "RNA", features = f, scale = F) + RotatedAxis()
plot(b)
rm(b)
dev.off()
#==============================================
rownames <- row.names(filtered_LT)
filtered_LT1 <- data.frame(filtered_LT[,colnames(filtered_LT)=="EDN1"])
row.names(filtered_LT1) <- rownames
names(filtered_LT1) <- "EDN1"


#take the 1st 10 genes
filtered_LT1$target <-row.names(filtered_LT1) 
filtered_LT1 <- filtered_LT1[order(filtered_LT1[,1],decreasing=TRUE),]
a <- top_n(filtered_LT1, 20, filtered_LT1[,1])
a <- row.names(a)
a<- a[1:20]

#remove duplicates
f<-c("EDN1", "EDNRA", a)
f <- f[!duplicated(f)]

pdf("EDN1_LR_top20_target_genes_dotplot_strict.pdf",width=9,height=4,paper="special")
b <- DotPlot(receivers_senders, assay = "RNA", features = f, scale = TRUE) + RotatedAxis()
plot(b)
rm(b)
dev.off()

pdf("EDN1_LR_top20_target_genes_dotplot_strict_unscalled.pdf",width=9,height=4,paper="special")
b <- DotPlot(receivers_senders, assay = "RNA", features = f, scale = F) + RotatedAxis()
plot(b)
rm(b)
dev.off()

#==============================================
rownames <- row.names(filtered_LT)
filtered_LT1 <- data.frame(filtered_LT[,colnames(filtered_LT)=="EDN1"])
row.names(filtered_LT1) <- rownames
names(filtered_LT1) <- "EDN1"


#take the 1st 10 genes
filtered_LT1$target <-row.names(filtered_LT1) 
filtered_LT1 <- filtered_LT1[order(filtered_LT1[,1],decreasing=TRUE),]
a <- top_n(filtered_LT1, 20, filtered_LT1[,1])
a <- row.names(a)
a<- a[1:20]

#remove duplicates
f<-c("EDN1", "EDNRA","EDNRB", a)
f <- f[!duplicated(f)]

pdf("EDN1_LR_top20_target_genes_dotplot_strict1.pdf",width=9,height=4,paper="special")
b <- DotPlot(receivers_senders, assay = "RNA", features = f, scale = TRUE) + RotatedAxis()
plot(b)
rm(b)
dev.off()

pdf("EDN1_LR_top20_target_genes_dotplot_strict1_unscalled.pdf",width=9,height=4,paper="special")
b <- DotPlot(receivers_senders, assay = "RNA", features = f, scale = F) + RotatedAxis()
plot(b)
rm(b)
dev.off()

#===================================================
f<-c("SST", "SSTR2","BCL2", "CITED2", "IRF2BP2")
levels(receivers_senders) <- rev(levels(receivers_senders))

pdf("SST_LR_top_target_genes_dotplot_for_text.pdf",width=4.5,height=4,paper="special")
b <- DotPlot(receivers_senders, assay = "RNA", features = f, scale = TRUE) + RotatedAxis()
plot(b)
rm(b)
dev.off()

pdf("SST_LR_top_target_genes_dotplot_text_unscalled.pdf",width=4.5,height=4,paper="special")
b <- DotPlot(receivers_senders, assay = "RNA", features = f, scale = F) + RotatedAxis()
plot(b)
rm(b)
dev.off()

#================================
#==============================================
rownames <- row.names(filtered_LT)
filtered_LT1 <- data.frame(filtered_LT[,colnames(filtered_LT)=="FGF7"])
row.names(filtered_LT1) <- rownames
names(filtered_LT1) <- "FGF7"


#take the 1st 10 genes
filtered_LT1$target <-row.names(filtered_LT1) 
filtered_LT1 <- filtered_LT1[order(filtered_LT1[,1],decreasing=TRUE),]
a <- top_n(filtered_LT1, 100, filtered_LT1[,1])
a <- row.names(a)
a<- a[1:100]

b <- lung_cluster.markers[lung_cluster.markers$cluster=="mes_cl18",]

c <- b[b$gene %in% a,]
DotPlot(receivers_senders, assay = "RNA", features = c$gene, scale = TRUE) + RotatedAxis()

