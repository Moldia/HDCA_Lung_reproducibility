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
library(tradeSeq)
library(slingshot)
library(niceRplots)
library(igraph)
library(scales)
library(rafalib)
library(fpc)
library(stats)
library(dendextend)


# set the name of the analyzed dataset
cell_type <- "prox"
load("C:/Users/alex/Desktop/163k_epi_final_strict/epi_163K_after_MAST_and_roc_analysis_data.RData")

DefaultAssay(lung) <- "integrated"

selected_clusters <- c(0, 6, 7, 11, 12)

sel <- lung$seurat_clusters %in% selected_clusters

# subset the original dataset
temp <- lung[,sel]

temp <- RunUMAP(temp, 
                assay = "integrated",
                n.neighbors = 10L, 
                min.dist = 0.5, 
                dims = 1:50, 
                spread = 5,
                metric = "cosine",
                repulsion.strength = 0.01,
                negative.sample.rate = 5,
                n.epochs = NULL,
                seed.use = 42L,
                learning.rate = 15,
                n.components = 2L,
                reduction.name = "umap_traj",
                reduction.key="umaptraj_")

DimPlot(temp, reduction = "umap_traj", cols = "alphabet2", label = TRUE, pt.size = 2)

# retrieve the coordinates
dimred <- temp@reductions$umap_traj@cell.embeddings
clustering <- factor(lung$seurat_clusters[sel])

dimvis <- temp@reductions$umap_traj@cell.embeddings


set.seed(1)
lineages <- getLineages(data = dimred,
                        clusterLabels = clustering,
                        start.clus = "6",
                        end.clus=c("0", "11"))


lineages
lineages@elementMetadata@listData[["reducedDim"]]<- dimvis
curves <- getCurves(lineages, thresh = 0.001, stretch = .000001, allow.breaks = T, approx_points = 30)
pseudotime <- slingPseudotime(curves, na = FALSE)
cellWeights <- slingCurveWeights(curves)


n_clusters <- nlevels(lung$seurat_clusters)
colors <- DiscretePalette(n_clusters, palette = "alphabet2")

cl_colors <- data.frame(c(0:(nlevels(lung$seurat_clusters)-1)), colors)
colnames(cl_colors) <- c("cluster", "color")

colors3 <- cl_colors[cl_colors$cluster %in% selected_clusters,]


png(filename = paste0(cell_type, "_trajectory_umap.png"), width = 800*2, height = 800*2, res = 300)
#mypar(4,4)
plot_meta(x = temp,red = "umap_traj", feat = "seurat_clusters",frame=F,label = T, col = colors3$color)
lines(curves@metadata[["curves"]][["Lineage1"]][["s"]], lwd = 5, col = c('black',"firebrick","navy","darkgreen"))
lines(curves@metadata[["curves"]][["Lineage2"]][["s"]], lwd = 5, col = c('black',"firebrick","navy","darkgreen"))
dev.off()

pdf(file = paste0(cell_type, "_trajectory_umap.pdf"), width = 4, height = 4)
#mypar(4,4)
plot_meta(x = temp,red = "umap_traj", feat = "seurat_clusters",frame=F,label = T, col = colors3$color)
lines(curves@metadata[["curves"]][["Lineage1"]][["s"]], lwd = 5, col = c('black',"firebrick","navy","darkgreen"))
lines(curves@metadata[["curves"]][["Lineage2"]][["s"]], lwd = 5, col = c('black',"firebrick","navy","darkgreen"))
dev.off()

DimPlot(temp, reduction = "umap_traj", cols = colors3$color, label = TRUE, pt.size = 2)

save.image(file="prox2_analysis_data.RData")

# because the datase is big, we randomly select 1000 cells
set.seed(1)
sel2 <- sample( colnames(lung[,sel]) , size = min(1000,ncol(temp)))
counts <- temp@assays$RNA@counts[Matrix::rowSums(temp@assays$RNA@counts[,sel2] > 1) > 50 ,sel2]
dim(counts)

#===========================================================================
#run TradeSeq

opt_k <- evaluateK(counts = counts,
                   pseudotime = pseudotime[sel2,],
                   cellWeights = cellWeights[sel2,],
                   verbose = T,parallel=T)

# create a singecellexperiment object
sce <- fitGAM(counts = counts,
              pseudotime = pseudotime[sel2,],
              cellWeights = cellWeights[sel2,],
              nknots = 7, verbose = T,parallel=T)

# run the analysis for gene identification
patternRes <- patternTest(sce)

# duplicate the result 
patternRes1 <- patternRes
method <- "ward.D2"

#===========================================================================
# all cell markers
# filter the results
patternRes1 <- patternRes
patternRes1$pvalue[is.na(patternRes1$pvalue)] <- 1
patternRes1 <- patternRes1[order(patternRes1$waldStat,decreasing = T),]
patternRes1$FDR <- p.adjust(patternRes1$pvalue,method = "BH")
patternRes1 <- patternRes1[ patternRes1$FDR < 0.0001 , ]
patternRes1 <- patternRes1[ patternRes1$fcMedian > 1 , ]
head(patternRes1,20)

# order all genes
gene_list <- rownames(patternRes1)[1:min(569,nrow(patternRes1))]
lims <- quantile(pseudotime,c(0.05,.95) )
res <- t(sapply(gene_list,
                pseudotime = pseudotime,
                cellWeights = cellWeights,
                lims=lims,
                mycounts=temp@assays$RNA@data[gene_list,],
                function(gene,mycounts,pseudotime,cellWeights,lims) {
                  ll <- lapply( 1:ncol(pseudotime),function(i){
                    l1 <- (cellWeights[,i] == 1 ) & (pseudotime[,i] > lims[1]) & (pseudotime[,i] < lims[2])
                    l1 <- colnames(temp)[l1] #temp or temp2, depend on what was used above
                    sm <- spline(smooth.spline( pseudotime[ l1, i ], mycounts[gene,l1], nknots = 20,spar = .8),n = 100)
                  })
                  return( c( rev( ll[[1]]$y ) ,rep(NA,1), ll[[2]]$y ) )
                }))


#SMOOTHED HEATMAP
png(filename = paste0(cell_type, "_trajectory_branch_markers.png"),width = 800*3,height = 800*9,res = 300)
x <- t(apply(res,1,function(x){scale(x,T,T)}))
to_order <- apply(x,1,function(x){ mean( (1:length(x))[order(x,decreasing = T)][1] )})
to_order <- apply( (res[,1:20] - res[,22:41]),1,function(x){ mean( (1:length(x))[order(x,decreasing = F)][1] )})
to_order <- apply(x,1,function(x){ which.max(x) })
o <- order(to_order)
o <- hclust(as.dist((1-cor(t(res),use = "pairwise.complete.obs"))/2),method = method )$order
par(mar=c(2,2,3,6))
image( t(x[o,][nrow(res):1,]) ,axes=F,frame=F , col=c("grey95",colorRampPalette(c("grey95","grey70","firebrick","firebrick4"))(90)))
text( par("usr")[2], seq(nrow(res),0,length.out = nrow(res)) / nrow(res) , labels = rownames(x)[o] , pos=4 ,xpd=T , cex=.25)
text( par("usr")[1], par("usr")[4]-diff(par("usr")[4:3])/30, labels = "lineage 2" , adj= c(0,0),xpd=T , cex=1 , pos=3)
text( par("usr")[2], par("usr")[4]-diff(par("usr")[4:3])/30, labels = "lineage 1" , adj= c(1,0),xpd=T , cex=1, pos=3)
text( mean(par("usr")[1:2]), par("usr")[4]-diff(par("usr")[4:3])/30, labels = "origin" , adj= c(1,0),xpd=T , cex=.8, pos=3)
x <- curves@metadata[["lineages"]][["Lineage1"]]
text( seq(mean(par("usr")[1:2]),par("usr")[1],length.out = length(x)), par("usr")[4], labels = x , adj= c(0,0),xpd=T , cex=.8 , pos=3)
x <- curves@metadata[["lineages"]][["Lineage2"]]
text( seq(mean(par("usr")[1:2]),par("usr")[2],length.out = length(x)), par("usr")[4], labels = x , adj= c(0,0),xpd=T , cex=.8 , pos=3)
dev.off()

#set number of gene modules
gene_dendro <- hclust(as.dist((1-cor(t(res),use = "pairwise.complete.obs"))/2),method = method )
plot(gene_dendro, cex = 0.1, hang = -1)
k=9
rect.hclust(gene_dendro, k = k, border = "red")

cb_TF <- clusterboot(as.dist((1-cor(t(res),use = "pairwise.complete.obs"))/2), distances=TRUE, B=1000, bootmethod=
                       "boot",clustermethod=disthclustCBI,
                     k=k, cut="number",method=method, showplots=FALSE, seed=1)

write.csv(cb_TF$bootmean, file = paste0(cell_type,"_stability_values_markers1.csv"))
gene_dendro1 <- cutree(gene_dendro, k=k)[gene_dendro$order]
write.csv(gene_dendro1, file = paste0(cell_type,"_marker_modules1.csv"))

#===============================================================================
plot(gene_dendro1)

# reorder the leaves of the dendrogram
a <- which(gene_dendro1==c(9))
b <- which(gene_dendro1==c(8))
c <- which(gene_dendro1==c(3))
d <- rev(which(gene_dendro1==c(1)))
e <- which(gene_dendro1==c(7))
f <- which(gene_dendro1==c(6))
g <- which(gene_dendro1==c(4))
h <- which(gene_dendro1==c(5))
i <- which(gene_dendro1==c(2))

plot(rotate(gene_dendro, c(a,b, c, d, e, f, g, h, i)),cex = 0.1, hang = -1)
rect.hclust(rotate(gene_dendro, c(a,b, c, d, e, f, g, h, i)), k = k, border = "red")
a <-rotate(gene_dendro, c(a,b, c, d, e, f, g, h, i))

png(filename = paste0(cell_type, "_trajectory_branch_markers_reordered.png"),width = 800*3,height = 800*9,res = 300)
x <- t(apply(res,1,function(x){scale(x,T,T)}))
to_order <- apply(x,1,function(x){ mean( (1:length(x))[order(x,decreasing = T)][1] )})
to_order <- apply( (res[,1:20] - res[,22:41]),1,function(x){ mean( (1:length(x))[order(x,decreasing = F)][1] )})
to_order <- apply(x,1,function(x){ which.max(x) })
o <- a[["order"]]
par(mar=c(2,2,3,6))
image( t(x[o,][nrow(res):1,]) ,axes=F,frame=F , col=c("grey95",colorRampPalette(c("grey95","grey70","firebrick","firebrick4"))(90)))
text( par("usr")[2], seq(nrow(res),0,length.out = nrow(res)) / nrow(res) , labels = rownames(x)[o] , pos=4 ,xpd=T , cex=.25)
text( par("usr")[1], par("usr")[4]-diff(par("usr")[4:3])/30, labels = "lineage 2" , adj= c(0,0),xpd=T , cex=1 , pos=3)
text( par("usr")[2], par("usr")[4]-diff(par("usr")[4:3])/30, labels = "lineage 1" , adj= c(1,0),xpd=T , cex=1, pos=3)
text( mean(par("usr")[1:2]), par("usr")[4]-diff(par("usr")[4:3])/30, labels = "origin" , adj= c(1,0),xpd=T , cex=.8, pos=3)
x <- curves@metadata[["lineages"]][["Lineage1"]]
text( seq(mean(par("usr")[1:2]),par("usr")[1],length.out = length(x)), par("usr")[4], labels = x , adj= c(0,0),xpd=T , cex=.8 , pos=3)
x <- curves@metadata[["lineages"]][["Lineage2"]]
text( seq(mean(par("usr")[1:2]),par("usr")[2],length.out = length(x)), par("usr")[4], labels = x , adj= c(0,0),xpd=T , cex=.8 , pos=3)
dev.off()

save.image(file="prox2_analysis_data.RData")
#================================================
# create a subset of each gene module to have equal number of genes for each of them
patternRes1 <- patternRes
patternRes1$pvalue[is.na(patternRes1$pvalue)] <- 1
patternRes1 <- patternRes1[order(patternRes1$waldStat,decreasing = T),]
patternRes1$FDR <- p.adjust(patternRes1$pvalue,method = "BH")
patternRes1 <- patternRes1[ patternRes1$FDR < 0.0001 , ]
patternRes1 <- patternRes1[ patternRes1$fcMedian > 1 , ]


gene_dendro2 <- data.frame(gene_dendro1)
gene_dendro2$gene <- row.names(gene_dendro2)
names(gene_dendro2)[1] <- "module"
patternRes1$gene <- row.names(patternRes1)
patternRes1 <- join(patternRes1, gene_dendro2, type="inner", by="gene")

# set the number of wanted genes/module
n=10
#===============================
cluster_order <- c(9, 8, 3, 1, 7, 6, 4, 5, 2)
d <- patternRes1[patternRes1$module=="9",]
d <- d[order(d$waldStat,decreasing = T),]
d <- d[1:min(n,nrow(d))]
gene_list1 <- d

for (i in c(8, 3, 1, 7, 6, 4, 5, 2)) {
  d <- patternRes1[patternRes1$module==i,]
  d <- d[order(d$waldStat,decreasing = T),]
  d <- d[1:min(n,nrow(d)),]
  gene_list1 <- rbind(gene_list1, d)
}

row.names(gene_list1) <- gene_list1$gene


gene_list2 <- data.frame(rownames(gene_list1))
names(gene_list2) <- "gene"
res1 <- data.frame(res)
res1$gene <- rownames(res1)
res1 <- join(gene_list2, res1, by= "gene", type="inner")
row.names(res1) <- res1$gene
res1 <- res1[-1]


png(filename = paste0(cell_type, "_trajectory_branch_markers_short1.png"),width = 800*1.5,height = 800*4.5,res = 300)
x <- t(apply(res1,1,function(x){scale(x,T,T)}))
to_order <- apply(x,1,function(x){ mean( (1:length(x))[order(x,decreasing = T)][1] )})
to_order <- apply( (res1[,1:20] - res1[,22:41]),1,function(x){ mean( (1:length(x))[order(x,decreasing = F)][1] )})
to_order <- apply(x,1,function(x){ which.max(x) })
o <- c(1:nrow(gene_list2))
par(mar=c(2,2,3,6))
image( t(x[o,][nrow(res1):1,]) ,axes=F,frame=F , col=c("grey95",colorRampPalette(c("grey95","grey70","firebrick","firebrick4"))(90)))
text( par("usr")[2], seq(nrow(res1),0,length.out = nrow(res1)) / nrow(res1) , labels = rownames(x)[o] , pos=4 ,xpd=T , cex=.8)
text( par("usr")[1], par("usr")[4]-diff(par("usr")[4:3])/30, labels = "lineage 2" , adj= c(0,0),xpd=T , cex=1 , pos=3)
text( par("usr")[2], par("usr")[4]-diff(par("usr")[4:3])/30, labels = "lineage 1" , adj= c(1,0),xpd=T , cex=1, pos=3)
text( mean(par("usr")[1:2]), par("usr")[4]-diff(par("usr")[4:3])/30, labels = "origin" , adj= c(1,0),xpd=T , cex=.8, pos=3)
x <- curves@metadata[["lineages"]][["Lineage1"]]
text( seq(mean(par("usr")[1:2]),par("usr")[1],length.out = length(x)), par("usr")[4], labels = x , adj= c(0,0),xpd=T , cex=.8 , pos=3)
x <- curves@metadata[["lineages"]][["Lineage2"]]
text( seq(mean(par("usr")[1:2]),par("usr")[2],length.out = length(x)), par("usr")[4], labels = x , adj= c(0,0),xpd=T , cex=.8 , pos=3)
dev.off()


write.csv(gene_list1, file = paste0("top_", n,"module_markers.csv"))

#================================================
#================================================
# run again for TFs
patternRes1 <- patternRes
patternRes1$pvalue[is.na(patternRes1$pvalue)] <- 1
patternRes1 <- patternRes1[order(patternRes1$waldStat,decreasing = T),]
patternRes1$FDR <- p.adjust(patternRes1$pvalue,method = "BH")
patternRes1 <- patternRes1[ patternRes1$FDR < 0.0001 , ]
patternRes1 <- patternRes1[ patternRes1$fcMedian > 1 , ]

gene_dendro2 <- data.frame(gene_dendro1)
gene_dendro2$gene <- row.names(gene_dendro2)
names(gene_dendro2)[1] <- "module"
patternRes1$gene <- row.names(patternRes1)
patternRes1 <- join(patternRes1, gene_dendro2, type="inner", by="gene")
row.names(patternRes1) <- patternRes1$gene

# select only the TFs
TFs <- read.csv("TFs.csv")
patternRes1 <- patternRes1[row.names(patternRes1) %in% TFs$gene,]

gene_dendro2 <- data.frame(gene_dendro1)
gene_dendro2$gene <- row.names(gene_dendro2)
names(gene_dendro2)[1] <- "module"
patternRes1$gene <- row.names(patternRes1)
patternRes1 <- join(patternRes1, gene_dendro2, type="inner", by="gene")
patternRes1 <-patternRes1[,1:7]

# set the number of wanted genes
n=5
#===============================
cluster_order <- c(9, 8, 3, 1, 7, 6, 4, 5, 2)
d <- patternRes1[patternRes1$module=="9",]
d <- d[order(d$waldStat,decreasing = T),]
d <- d[1:min(n,nrow(d)),]
gene_list1 <- d

for (i in c(8, 3, 1, 7, 6, 4, 5, 2)) {
  d <- patternRes1[patternRes1$module==i,]
  d <- d[order(d$waldStat,decreasing = T),]
  d <- d[1:min(n,nrow(d)),]
  gene_list1 <- rbind(gene_list1, d)
}
gene_list1 <- na.omit(gene_list1)

gene_list1 <- na.omit(gene_list1)
row.names(gene_list1) <- gene_list1$gene


gene_list2 <- data.frame(rownames(gene_list1))
names(gene_list2) <- "gene"
res1 <- data.frame(res)
res1$gene <- rownames(res1)
res1 <- join(gene_list2, res1, by= "gene", type="inner")
row.names(res1) <- res1$gene
res1 <- res1[-1]


png(filename = paste0(cell_type, "_trajectory_branch_TFs_short1.png"),width = 800*1.5,height = 800*3,res = 300)
x <- t(apply(res1,1,function(x){scale(x,T,T)}))
to_order <- apply(x,1,function(x){ mean( (1:length(x))[order(x,decreasing = T)][1] )})
to_order <- apply( (res1[,1:20] - res1[,22:41]),1,function(x){ mean( (1:length(x))[order(x,decreasing = F)][1] )})
to_order <- apply(x,1,function(x){ which.max(x) })
o <- c(1:nrow(gene_list2))
par(mar=c(2,2,3,6))
image( t(x[o,][nrow(res1):1,]) ,axes=F,frame=F , col=c("grey95",colorRampPalette(c("grey95","grey70","firebrick","firebrick4"))(90)))
text( par("usr")[2], seq(nrow(res1),0,length.out = nrow(res1)) / nrow(res1) , labels = rownames(x)[o] , pos=4 ,xpd=T , cex=.8)
text( par("usr")[1], par("usr")[4]-diff(par("usr")[4:3])/30, labels = "lineage 2" , adj= c(0,0),xpd=T , cex=1 , pos=3)
text( par("usr")[2], par("usr")[4]-diff(par("usr")[4:3])/30, labels = "lineage 1" , adj= c(1,0),xpd=T , cex=1, pos=3)
text( mean(par("usr")[1:2]), par("usr")[4]-diff(par("usr")[4:3])/30, labels = "origin" , adj= c(1,0),xpd=T , cex=.8, pos=3)
x <- curves@metadata[["lineages"]][["Lineage1"]]
text( seq(mean(par("usr")[1:2]),par("usr")[1],length.out = length(x)), par("usr")[4], labels = x , adj= c(0,0),xpd=T , cex=.8 , pos=3)
x <- curves@metadata[["lineages"]][["Lineage2"]]
text( seq(mean(par("usr")[1:2]),par("usr")[2],length.out = length(x)), par("usr")[4], labels = x , adj= c(0,0),xpd=T , cex=.8 , pos=3)
dev.off()


write.csv(gene_list1, file = paste0("top_", n,"module_TFs.csv"))

levels(temp) <- c( "0","6", "7", "12", "11")
DefaultAssay(temp) <- "RNA"
DotPlot(temp, assay = "RNA", features = gene_list2$gene, scale = TRUE) + RotatedAxis()

#============================================================================
DefaultAssay(lung) <- "RNA"
levels(lung) <- c("14", "11", "12", "7", "6", "0", "4", "1", "10", "2", "9", "3", "8", "13", "5")

for (i in 1:length(gene_list)){
  prefix ="C:/Users/alex/Desktop/163k_epi_final_strict/prox2_trajectory/heatmap_genes/violinplot_"
  name = gene_list[[i]]
  pdf(paste0(prefix,name, ".pdf"), width = 12, height = 3, paper = "special")
  pdf <- VlnPlot(lung, features = gene_list[[i]], pt.size = 0, cols = colors) + NoLegend()
  plot(pdf)
  dev.off()
  rm(pdf)
} 

#===========================================================================
# plot-inspect known markers
DefaultAssay(temp) <- "RNA"
FeaturePlot(temp, reduction = "umap_traj", features = "CALCA", label = F, pt.size = 1)
FeaturePlot(temp, reduction = "umap_traj", features = "ASCL1", label = F, pt.size = 1)
FeaturePlot(temp, reduction = "umap_traj", features = "GRP", label = F, pt.size = 1)
FeaturePlot(temp, reduction = "umap_traj", features = "GHRL", label = F, pt.size = 1)
FeaturePlot(temp, reduction = "umap_traj", features = "ACSL1", label = F, pt.size = 1)
FeaturePlot(temp, reduction = "umap_traj", features = "NKX2-2", label = F, pt.size = 1)
FeaturePlot(temp, reduction = "umap_traj", features = "CALCA", label = F, pt.size = 1)
FeaturePlot(temp, reduction = "umap_traj", features = "TUBB3", label = F, pt.size = 1)
FeaturePlot(temp, reduction = "umap_traj", features = "NEUROD1", label = F, pt.size = 1)

FeaturePlot(temp, reduction = "umap_traj", features = "SALL4", label = F, pt.size = 1)
levels(temp) <- c("6", "0", "7", "12", "11")
colors3 <- colors3[c(2,1,3, 5, 4),]

#3x9 inches
VlnPlot(temp, features = c("ISL1"), cols = colors3$color, pt.size = 0)

DefaultAssay(lung) <- "RNA"
VlnPlot(lung, features = c("DLL3"), pt.size = 0)
gene = "ASCL1"
mv <- as.integer(max(lung@assays[["RNA"]]@data[row.names(lung@assays[["RNA"]]@data)==gene]))
FeaturePlot(lung, reduction = "umap2d", features = gene, label = F, pt.size = 1, order = F)+ scale_colour_gradientn(breaks=c(0, 1, mv), colors = c( "wheat", "gray","blue"))

DefaultAssay(temp) <- "RNA"
min_cluster <- min(table(temp@active.ident))
lung_cluster.markers <- FindAllMarkers(temp, test.use = "MAST", min.diff.pct = 0.1, max.cells.per.ident = min_cluster, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.1)
# calculate the percentage difference of the positive cells in the indicated cluster, in comparison to the others
lung_cluster.markers$Dpct <- (lung_cluster.markers$pct.1 - lung_cluster.markers$pct.2)

c <-subset(temp, idents =6 )
counts <- as.matrix(c[["RNA"]]@data)
m <- data.frame(rowMeans(counts))
gene <- row.names(c)
m <- cbind(gene, m)
names(m)[2] <- paste0("cluster", "6") 

for (i in c("0", "7", "12", "11")) {
  c1 <-subset(lung, idents =i )
  counts1 <- as.matrix(c1[["RNA"]]@data)
  m1 <- data.frame(rowMeans(counts1))
  names(m1)[1] <- paste0("cluster", (i))
  m <- cbind(m, m1)
}
rm(c, counts, counts1, c1, m1)


lung_cluster.markers <-join(lung_cluster.markers, m, by="gene")
write.csv(lung_cluster.markers, file="prox2_markers1.csv")

# run ROC differential expression analysis
DefaultAssay(temp) <-"RNA"
lung_roc_cluster.markers <- FindAllMarkers(temp, test.use = "roc", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.1)
# calculate the percentage difference of the positive cells in the indicated cluster, in comparison to the others
lung_roc_cluster.markers$Dpct <- (lung_roc_cluster.markers$pct.1 - lung_roc_cluster.markers$pct.2)

lung_roc_cluster.markers <-join(lung_roc_cluster.markers, m, by="gene")
write.csv(lung_roc_cluster.markers, file="prox2_roc_markers.csv")

lung_cluster.markers1 <- lung_cluster.markers
lung_cluster.markers <- lung_cluster.markers[lung_cluster.markers$p_val_adj<0.001, ]
# analysis based of the top markers according to roc power (10 genes) and then avg_log2fc(5 genes)
top2 <- lung_cluster.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top2 <- top2 %>% group_by(cluster) %>% top_n(n = 5, wt = Dpct)
top2 <-top2[, c(7,6)]
top2 <- top2[!duplicated(top2[ , "gene"]),]

levels(temp) <- c("11", "12", "7",  "6", "0")
new_order <- c( "0", "6",  "7", "12", "11")

e <- top2[top2$cluster %in% new_order[1],]
for (i in 2:length(new_order)){
  e1 <- top2[top2$cluster %in% new_order[i],]
  e <- rbind(e, e1)
}

e <- c(e$gene)
e <- e[!duplicated(e)]
#3x9 inches
#prox2_top5_roc_markers_dotplot
DotPlot(temp, assay = "RNA", features = e, scale = TRUE) + RotatedAxis()

#===========================================================================
# plot violin plots for the cluster-7 (NE-progenitor)
DefaultAssay(lung) <- "RNA"
levels(lung) <- c("14", "11", "12", "7", "6", "0", "4", "1", "10", "2", "9", "3", "8", "13", "5")

prox2_markers <- lung_cluster.markers1
gene_list <- prox2_markers[prox2_markers$cluster==7,]
gene_list <-gene_list$gene

for (i in 1:length(gene_list)){
  prefix ="C:/Users/alex/Desktop/163k_epi_final_strict/prox2_trajectory_strict/cl7_markers/"
  name = gene_list[[i]]
  pdf(paste0(prefix,name, ".pdf"), width = 12, height = 3, paper = "special")
  pdf <- VlnPlot(lung, features = gene_list[[i]], pt.size = 0, cols = colors) + NoLegend()
  plot(pdf)
  dev.off()
  rm(pdf)
} 

for (i in 1:length(gene_list)){
  prefix ="C:/Users/alex/Desktop/163k_epi_final_strict/prox2_trajectory_strict/cl7_markers/"
  name = gene_list[[i]]
  tiff(paste0(prefix,name, ".tiff"), width = 1200, height = 300, compression = "lzw")
  pdf <- VlnPlot(lung, features = gene_list[[i]], pt.size = 0, cols = colors) + NoLegend()
  plot(pdf)
  dev.off()
  rm(pdf)
} 

#=================================================================
#plot notch genes
DefaultAssay(lung) <- "RNA"
levels(lung) <- c("14", "11", "12", "7", "6", "0", "4", "1", "10", "2", "9", "3", "8", "13", "5")

DotPlot(lung, assay = "RNA", features = c("NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4",
                                           "DLL1", "DLL3", "DLL4", "JAG1", "JAG2", 
                                           "HES1", "HEY1", "HEYL", "HES6", "REST"), scale = FALSE) + RotatedAxis()

#3x11.5
DotPlot(temp, assay = "RNA", features = c("NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4", #receptors
                                          "HES1", "HES2", "HES3", "HES5", "HES6", "HES7", "HEY1", "HEY2", "HEYL", "NRARP", #targets
                                          "DLL1", "DLL4", "JAG1", "JAG2", #ligands
                                          "POFUT1", "POFUT2", "LFNG", "MFNG", "RFNG", "MAML1", "MAML2", "MAML3", "MAMLD1", "RBPJ", #transducers
                                          "DLL3", "DLK1", "DLK2", "NUMB", "NUMBL", # inhibitors
                                          
                                          "REST", "YAP1", "SCGB3A2", "MYCL", "ASCL1", "GRP", "NEUROD1", "GHRL"), scale = T) + RotatedAxis()

#3x10
# notch genes from https://pubmed.ncbi.nlm.nih.gov/31585080/
DotPlot(temp, assay = "RNA", features = c("NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4",
                                          "DLL1", "DLL3", "DLL4", "JAG1", "JAG2", 
                                          "HES1", "HEY1", "HEYL", "HES6", "REST", "SCGB3A2"), scale = T) + RotatedAxis()


DotPlot(temp, assay = "RNA", features = c("ATP4A", "ATP6V0D1", "ATP4B", "ATP6V1C1", "ATP6V1D", "ATP6V1E1", "ATP6V1F", "ATP6V1G1", 
                                          "ATP6V0E1", "ATP6V1G2", "ATP6V0C", "ATP6V0B", "CACNA1A", "STX1A", "CLTA", "CLTB", "CLTC", 
                                          "CACNA1B", "SLC6A1", "SLC6A2", "SLC6A3", "SLC6A4", "SLC6A5", "SLC6A9", "SLC6A7", "SLC6A11", 
                                          "SLC6A12", "SLC6A13", "SLC1A1", "SLC1A2", "SLC1A3", "SLC1A6", "SLC1A7", "NSF", "RAB3A", "SLC18A1", 
                                          "SLC18A2", "STX1B", "TFAP2A", "AP2B1", "AP2M1", "AP2S1", "SLC17A6", "SLC17A7", "SLC17A8", "VAMP2", 
                                          "SLC18A3", "SLC32A1", "SYT1", "RIMS1", "STXBP1", "UNC13A", "UNC13B", "UNC13C", "CPLX1", "CPLX2", 
                                          "CPLX3", "CPLX4", "NAPA", "SNAP25", "DNM2"), scale = F) + RotatedAxis()

#========================================================
# analysis based of the top markers according to avg_log2fc(20 genes) and then difference in the percent of positive cells (10 genes)
top2 <- lung_cluster.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
top2 <- top2 %>% group_by(cluster) %>% top_n(n = 10, wt = Dpct)
top2 <-top2[, c(7,6)]
top2 <- top2[!duplicated(top2[ , "gene"]),]

levels(temp) <- c("11", "12", "7",  "6", "0")
new_order <- c( "0", "6",  "7", "12", "11")

e <- top2[top2$cluster %in% new_order[1],]
for (i in 2:length(new_order)){
  e1 <- top2[top2$cluster %in% new_order[i],]
  e <- rbind(e, e1)
}

e <- c(e$gene)
e <- e[!duplicated(e)]
#3x12 inches
#prox2_top10_markers_dotplot
DotPlot(temp, assay = "RNA", features = e, scale = TRUE) + RotatedAxis()

e <- lung_roc_cluster.markers %>% group_by(cluster) %>% top_n(n = 5, wt = power)
e <- c(e$gene)
e <- e[!duplicated(e)]
DotPlot(temp, assay = "RNA", features = e, scale = TRUE) + RotatedAxis()

#========================================================
save.image(file="prox2_analysis_data.RData")
#========================================================
# plot top5 TFs
TFs <- read.csv("TFs.csv")
prox2_markers <- read.csv("C:/Users/alex/Desktop/163k_epi_final_strict/prox2_trajectory/prox2_markers1.csv", row.names=1)
#lung_cluster.markers1 <- prox2_markers[prox2_markers$p_val_adj<0.001, ]
lung_cluster.markers1 <- prox2_markers

# create a new file by merging the TF-list and the enriched genes in each cluster, to identify the TFs, only
lung_cluster_tfs <- join(lung_cluster.markers1, TFs, by="gene", type="inner")
write.csv(lung_cluster_tfs, file="prox2_TFs.csv")

# analysis based of the top markers according to avg_log2fc(10 TFs) and then difference in the percent of positive cells (5 TFs)
top3 <- lung_cluster_tfs %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top3 <- top3 %>% group_by(cluster) %>% top_n(n = 5, wt = Dpct)
top3 <-top3[, c(7,6)]
top3 <- top3[!duplicated(top3[ , "gene"]),]

levels(temp) <- c("11", "12", "7", "6", "0")
new_order <- c(  "0","6", "7", "12", "11")

f <- top3[top3$cluster %in% new_order[1],]
for (i in 2:length(new_order)){
  f1 <- top3[top3$cluster %in% new_order[i],]
  f <- rbind(f, f1)
}

f <- c(f$gene)
f <- f[!duplicated(f)]
#3x9 inches 
DotPlot(temp, assay = "RNA", features = f, scale = TRUE) + RotatedAxis()
DotPlot(temp, assay = "RNA", features = f, scale = FALSE) + RotatedAxis()


# create a new file by merging the TF-list and the enriched genes in each cluster, to identify the TFs, only
lung_roc_cluster_tfs <- join(lung_roc_cluster.markers, TFs, by="gene", type="inner")
top5_roc_cluster_tfs <- lung_roc_cluster_tfs %>% group_by(cluster) %>% top_n(n = 5, wt = power)

f <- c(top5_roc_cluster_tfs$gene)
f <- f[!duplicated(f)]
DotPlot(temp, assay = "RNA", features = f, scale = TRUE) + RotatedAxis()

#==============================================================================
#=============================================================================
# plot violin plots for all transcription factors
DefaultAssay(lung) <- "RNA"
levels(lung) <- c("14", "11", "12", "7", "6", "0", "4", "1", "10", "2", "9", "3", "8", "13", "5")

for (i in 1:nrow(lung_cluster_tfs)){
  prefix ="C:/Users/alex/Desktop/163k_epi_final_strict/prox2_trajectory_strict/TF_genes/violinplot_"
  name = lung_cluster_tfs$gene[i]
  pdf(paste0(prefix,name, ".pdf"), width = 12, height = 3, paper = "special")
  pdf <- VlnPlot(lung, features = lung_cluster_tfs$gene[i], pt.size = 0, cols = colors) + NoLegend()
  plot(pdf)
  dev.off()
  rm(pdf)
} 

for (i in 1:nrow(lung_cluster_tfs)){
  prefix ="C:/Users/alex/Desktop/163k_epi_final_strict/prox2_trajectory_strict/TF_genes/violinplot_"
  name = lung_cluster_tfs$gene[i]
  tiff(paste0(prefix,name, ".tiff"), width = 1200, height = 300, compression = "lzw")
  pdf <- VlnPlot(lung, features = lung_cluster_tfs$gene[i], pt.size = 0, cols = colors) + NoLegend()
  plot(pdf)
  dev.off()
  rm(pdf)
} 

VlnPlot(lung, features = "MYC", pt.size = 0, cols = colors) + NoLegend()


for (i in 1:nrow(lung_roc_cluster_tfs)){
  prefix ="C:/Users/alex/Desktop/163k_epi_final_strict/prox2_trajectory_strict/TF_genes/violinplot_"
  name = lung_roc_cluster_tfs$gene[i]
  pdf(paste0(prefix,name, ".pdf"), width = 12, height = 3, paper = "special")
  pdf <- VlnPlot(lung, features = lung_roc_cluster_tfs$gene[i], pt.size = 0, cols = colors) + NoLegend()
  plot(pdf)
  dev.off()
  rm(pdf)
} 

for (i in 1:nrow(lung_roc_cluster_tfs)){
  prefix ="C:/Users/alex/Desktop/163k_epi_final_strict/prox2_trajectory_strict/TF_genes/violinplot_"
  name = lung_roc_cluster_tfs$gene[i]
  tiff(paste0(prefix,name, ".tiff"), width = 1200, height = 300, compression = "lzw")
  pdf <- VlnPlot(lung, features = lung_roc_cluster_tfs$gene[i], pt.size = 0, cols = colors) + NoLegend()
  plot(pdf)
  dev.off()
  rm(pdf)
} 


#========================================================
#========================================================
# plot top5 coTFs
coTFs <- read.csv("TF_cofactors.csv")
prox2_markers <- read.csv("C:/Users/alex/Desktop/163k_epi_final_strict/prox2_trajectory/prox2_markers.csv", row.names=1)
lung_cluster.markers1 <- prox2_markers[prox2_markers$p_val_adj<0.01, ]

# create a new file by merging the coTF-list and the enriched genes in each cluster, to identify the coTFs, only
lung_cluster_cotfs <- join(lung_cluster.markers1, coTFs, by="gene", type="inner")
write.csv(lung_cluster_cotfs, file="prox2_coTFs.csv")

# analysis based of the top markers according to avg_log2fc(10 genes) and then difference in the percent of positive cells (5 genes)
top3 <- lung_cluster_cotfs %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top3 <- top3 %>% group_by(cluster) %>% top_n(n = 5, wt = Dpct)
top3 <-top3[, c(7,6)]
top3 <- top3[!duplicated(top3[ , "gene"]),]

levels(temp) <- c("11", "12", "7", "6", "0")
new_order <- c(  "0","6", "7", "12", "11")

f <- top3[top3$cluster %in% new_order[1],]
for (i in 2:length(new_order)){
  f1 <- top3[top3$cluster %in% new_order[i],]
  f <- rbind(f, f1)
}

f <- c(f$gene)
f <- f[!duplicated(f)]
#3x8 inches 
DotPlot(temp, assay = "RNA", features = f, scale = TRUE) + RotatedAxis()
DotPlot(temp, assay = "RNA", features = f, scale = FALSE) + RotatedAxis()

#========================================================
# plot top5 secr
secr <- data.frame(read.delim("protein_class_secreted.tsv"))
secr <- secr[1:10]

# create a new file by merging the secreted protein list and the enriched genes in each cluster, to identify the secreted proteins, only
lung_cluster_secr <- join(lung_cluster.markers, secr, by="gene", type="inner")

# analysis based of the top markers according to avg_log2fc(10 genes) and then difference in the percent of positive cells (5 genes)
top4 <- lung_cluster_secr %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top4 <- top4 %>% group_by(cluster) %>% top_n(n = 5, wt = Dpct)
top4 <-top4[, c(7,6)]
top4 <- top4[!duplicated(top4[ , "gene"]),]

levels(temp) <- c("11", "12", "7", "0", "6")
new_order <- c( "6", "0", "7", "12", "11")

h <- top4[top4$cluster %in% new_order[1],]
for (i in 2:length(new_order)){
  h1 <- top4[top4$cluster %in% new_order[i],]
  h <- rbind(h, h1)
}

h <- c(h$gene)
h <- h[!duplicated(h)]
#3x8 inches 
DotPlot(temp, assay = "RNA", features = h, scale = TRUE) + RotatedAxis()

#========================================================
# plot top5 CDs
cds <- data.frame(read.delim("protein_class_CD.tsv"))
cds <- cds[1:10]

# create a new file by merging the CD-list and the enriched genes in each cluster, to identify the CDs, only
lung_cluster_cds <- join(lung_cluster.markers, cds, by="gene", type="inner")

# analysis based of the top markers according to avg_log2fc(20 genes) and then difference in the percent of positive cells (10 genes)
top5 <- lung_cluster_cds %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top5 <- top5 %>% group_by(cluster) %>% top_n(n = 5, wt = Dpct)
top5 <-top5[, c(7,6)]
top5 <- top5[!duplicated(top5[ , "gene"]),]

levels(temp) <- c("11", "12", "7", "0", "6")
new_order <- c( "6", "0", "7", "12", "11")

g <- top5[top5$cluster %in% new_order[1],]
for (i in 2:length(new_order)){
  g1 <- top5[top5$cluster %in% new_order[i],]
  g <- rbind(g, g1)
}

g <- c(g$gene)
g <- g[!duplicated(g)]
#5x20 inches 
DotPlot(temp, assay = "RNA", features = g, scale = TRUE) + RotatedAxis()


#========================================================
# plot top5 kinases
kinases <- data.frame(read.delim("protein_class_kinases.tsv"))
kinases <- kinases[1:10]

# create a new file by merging the kinases and the enriched genes in each cluster, to identify the kiases, only
lung_cluster_kinases <- join(lung_cluster.markers1, kinases, by="gene", type="inner")

# analysis based of the top markers according to avg_log2fc(10 genes) and then difference in the percent of positive cells (5 genes)
top5 <- lung_cluster_kinases %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top5 <- top5 %>% group_by(cluster) %>% top_n(n = 5, wt = Dpct)
top5 <-top5[, c(7,6)]
top5 <- top5[!duplicated(top5[ , "gene"]),]

levels(temp) <- c("11", "12", "7", "6", "0")
new_order <- c(  "0","6", "7", "12", "11")

g <- top5[top5$cluster %in% new_order[1],]
for (i in 2:length(new_order)){
  g1 <- top5[top5$cluster %in% new_order[i],]
  g <- rbind(g, g1)
}

g <- c(g$gene)
g <- g[!duplicated(g)]
#3x8 inches 
DotPlot(temp, assay = "RNA", features = g, scale = TRUE) + RotatedAxis()

save.image(file="prox2_analysis_data.RData")

DefaultAssay(lung) <- "RNA"
FeaturePlot(lung, reduction = "umap2d", features = "SALL4", label = F, pt.size = 1)


#------------------------------------------------------------------
# analysis of NE cells by pair-wise comparison
ne_clusters <- c(11, 12)

sel1 <- lung$seurat_clusters %in% ne_clusters

# suset the original dataset
temp1 <- lung[,sel1]
DefaultAssay(temp1) <- "RNA"
min_cluster <- min(table(temp1@active.ident))
cl11_vs_cl12_markers <- FindAllMarkers(temp1, test.use = "MAST", min.diff.pct = 0.1, max.cells.per.ident = min_cluster, only.pos = T, min.pct = 0.25, logfc.threshold = 0.1)
# calculate the percentage difference of the positive cells in the indicated cluster, in comparison to the others
cl11_vs_cl12_markers$Dpct <- (cl11_vs_cl12_markers$pct.1 - cl11_vs_cl12_markers$pct.2)
cl11_vs_cl12_markers$gene <- row.names(cl11_vs_cl12_markers)
cl11_vs_cl12_markers <-join(cl11_vs_cl12_markers, m, by="gene")
write.csv(cl11_vs_cl12_markers, file="cl11_vs_cl12_markers.csv")

cl11_vs_cl12_markers <- cl11_vs_cl12_markers[cl11_vs_cl12_markers$p_val_adj<0.001, ]
write.csv(cl11_vs_cl12_markers, file="cl11_vs_cl12_markers_filtered.csv")

cl11_vs_cl12_markers_025logfc <- cl11_vs_cl12_markers[cl11_vs_cl12_markers$avg_log2FC>0.25, ]
write.csv(cl11_vs_cl12_markers_025logfc, file="cl11_vs_cl12_markers_filtered_025logfc.csv")
#--------------------------------------------------------------------------------------------------------------
cl11_vs_cl12_tfs <- join(cl11_vs_cl12_markers, TFs, by="gene", type="inner")
write.csv(cl11_vs_cl12_tfs, file="cl11_vs_cl12_tfs.csv")

cl11_vs_cl12_cotfs <- join(cl11_vs_cl12_markers, coTFs, by="gene", type="inner")
write.csv(cl11_vs_cl12_cotfs, file="cl11_vs_cl12_cotfs.csv")

cl11_vs_cl12_cds <- join(cl11_vs_cl12_markers, cds, by="gene", type="inner")
write.csv(cl11_vs_cl12_cds, file="cl11_vs_cl12_cds.csv")

cl11_vs_cl12_kinases <- join(cl11_vs_cl12_markers, kinases, by="gene", type="inner")
write.csv(cl11_vs_cl12_kinases, file="cl11_vs_cl12_kinases.csv")

secr <- data.frame(read.delim("protein_class_secreted.tsv"))
secr <- secr[1:10]
cl11_vs_cl12_secr <- join(cl11_vs_cl12_markers, secr, by="gene", type="inner")
write.csv(cl11_vs_cl12_secr, file="cl11_vs_cl12_secreted.csv")

# neuropeptides were retrieved from http://proteomics.ucsd.edu/Software/NeuroPedia/#Downloads
neuropeptides <- read_excel("neuropedia_neuropeptides_human.xlsx")
neuro <- levels(as.factor(neuropeptides$"Gene Name"))
cl11_vs_cl12_neuro <- cl11_vs_cl12_markers[cl11_vs_cl12_markers$gene %in% neuro,]
write.csv(cl11_vs_cl12_neuro, file="cl11_vs_cl12_neuropeptides.csv")
#========================================================
save.image(file="prox2_analysis_data.RData")
#-------------------------------------------------------------------------------------------------------------
# analysis based of the top markers according to avg_log2fc(50 genes) and then difference in the percent of positive cells (25 genes)
top6 <- cl11_vs_cl12_markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
top6 <- top6 %>% group_by(cluster) %>% top_n(n = 25, wt = Dpct)
top6 <-top6[, c(7,6)]
top6 <- top6[!duplicated(top6[ , "gene"]),]

levels(temp1) <- c("11", "12")
new_order <- c("12", "11")

h <- top6[top6$cluster %in% new_order[1],]
for (i in 2:length(new_order)){
  h1 <- top6[top6$cluster %in% new_order[i],]
  h <- rbind(h, h1)
}
write.csv(h, file="Ext_Fig5_F_gene_list.csv")
h <- c(h$gene)
h <- h[!duplicated(h)]
#3x13 inches
#prox2_top5_markers_dotplot
DotPlot(temp1, assay = "RNA", features = h, scale = TRUE) + RotatedAxis()

#==============================
ne11_clusters <- c(12)

sel2 <- lung$seurat_clusters %in% ne11_clusters

colors_age <- colors[2:11]

# suset the original dataset
temp2 <- lung[,sel2]
temp2@meta.data[["age"]] <- factor(temp2@meta.data[["age"]])
DefaultAssay(temp2) <- "RNA"
temp2 <-SetIdent(temp2, value = "age")
#4x6
DotPlot(temp2, assay = "RNA", features = c("PROX1", "DPP10", "GRP", "ASCL1", "CALCA", "SST", "GHRL", "ACSL1", "VSTM2L", "CFC1", "RFX6", "NKX2-2", "ARX", "PCSK1"), scale = FALSE) + RotatedAxis()

VlnPlot(temp2, features = "GRP", pt.size = 0, cols = colors_age)
#==============================================
levels(temp) <- c("11", "12", "7", "0", "6")
DotPlot(temp, assay = "RNA", features = c("GHRL", "CRH", "GRP", "CALCA", "NXPH4", "SST", "TAC3", "CALCB"), scale = T) + RotatedAxis()

