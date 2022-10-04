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
cell_type <- "neuro"
load("C:/Users/alex/Desktop/neuronal//neuro_163K_after_MAST_and_roc_analysis_data.RData")

DefaultAssay(lung) <- "integrated"

selected_clusters <- c(4, 0, 3, 2, 6)

sel <- lung$seurat_clusters %in% selected_clusters

# subset the original dataset
temp <- lung[,sel]

temp <- RunUMAP(temp, 
                assay = "integrated",
                n.neighbors = 15L, 
                min.dist = 0.5, 
                dims = 1:50, 
                spread = 1,
                metric = "cosine",
                repulsion.strength = 0.01,
                negative.sample.rate = 20,
                n.epochs = NULL,
                seed.use = 42L,
                learning.rate = 15,
                n.components = 2L,
                reduction.name = "umap_traj",
                reduction.key="umaptraj_")

DimPlot(temp, reduction = "umap_traj", cols = "alphabet2", label = TRUE, pt.size = 2)


dimred <- temp@reductions$umap_traj@cell.embeddings
temp@meta.data[["seurat_clusters"]] <- factor(temp@meta.data[["seurat_clusters"]])
clustering <- factor(temp@meta.data[["seurat_clusters"]])

dimvis <- temp@reductions$umap_traj@cell.embeddings


set.seed(1)
lineages <- getLineages(data = dimred,
                        clusterLabels = clustering,
                        start.clus = "4",
                        end.clus=c("6"))


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
dev.off()

pdf(file = paste0(cell_type, "_trajectory_umap.pdf"), width = 4, height = 4)
#mypar(4,4)
plot_meta(x = temp,red = "umap_traj", feat = "seurat_clusters",frame=F,label = T, col = colors3$color)
lines(curves@metadata[["curves"]][["Lineage1"]][["s"]], lwd = 5, col = c('black',"firebrick","navy","darkgreen"))
dev.off()

DimPlot(temp, reduction = "umap_traj", cols = colors3$color, label = TRUE, pt.size = 2)

FeaturePlot(temp, reduction = "umap_traj", features = "HBA1")

# we run the analysis in all cells because the dataset is small
set.seed(1)
sel2 <- sample( colnames(lung[,sel]) , size = min(1000,ncol(temp)))
counts <- temp@assays$RNA@counts[Matrix::rowSums(temp@assays$RNA@counts[,sel2] > 1) > 50 ,sel2]
dim(counts)


#===========================================================================
#run Tradeseq

opt_k <- evaluateK(counts = counts,
                   pseudotime = pseudotime[sel2,],
                   cellWeights = cellWeights[sel2,],
                   verbose = T,parallel=T)

# create a singecellexperiment object
sce <- fitGAM(counts = counts,
              pseudotime = pseudotime[sel2,],
              cellWeights = cellWeights[sel2,],
              nknots = 4, verbose = T,parallel=T)

# run analysis for gene identification for one trajectory
assoRes <- associationTest(sce)

# duplicate the result 
patternRes1 <- assoRes
method <- "ward.D2"
#method <- "complete"

save.image("neuronal_trajectory.RData")
#==============================================================================
# filter the results
patternRes1 <-assoRes
patternRes1$pvalue[is.na(patternRes1$pvalue)] <- 1
patternRes1 <- patternRes1[order(patternRes1$waldStat,decreasing = T),]
patternRes1$FDR <- p.adjust(patternRes1$pvalue,method = "BH")
patternRes1 <- patternRes1[ patternRes1$FDR < 0.0001 , ]
patternRes1 <- patternRes1[ patternRes1$meanLogFC > 1 , ]
head(patternRes1,20)

# select all genes and order them
gene_list <- rownames(patternRes1)[1:min(434,nrow(patternRes1))]
lims <- quantile(pseudotime,c(0.02,.98) )
res <- t(sapply(gene_list,
                pseudotime = pseudotime,
                cellWeights = cellWeights,
                lims=lims,
                mycounts=temp@assays$RNA@data[gene_list,],
                function(gene,mycounts,pseudotime,cellWeights,lims) {
                  ll <- lapply( 1:ncol(pseudotime),function(i){
                    l1 <- (cellWeights[,i] == 1 ) & (pseudotime[,i] > lims[1])
                    l1 <- colnames(temp)[l1] #temp or temp2, depend on what was used above
                    sm <- spline(smooth.spline( pseudotime[ l1, i ], mycounts[gene,l1], nknots = 50,spar = .8),n = 100)
                  })
                  return( c(rep(NA,1),  ll[[1]]$y ) )
                }))


#SMOOTHED HEATMAP
png(filename = paste0(cell_type, "_trajectory_branch_markers.png"),width = 800*3,height = 800*9,res = 300)
x <- t(apply(res,1,function(x){scale(x,T,T)}))
#to_order <- apply(x,1,function(x){ mean( (1:length(x))[order(x,decreasing = T)][1] )})
#to_order <- apply( (res[,1:40]),1,function(x){ mean( (1:length(x))[order(x,decreasing = F)][1] )})
to_order <- apply(x,1,function(x){ which.max(x) })
o <- order(to_order)
#o <- hclust(as.dist((1-cor(t(res),use = "pairwise.complete.obs"))/2),method = method )$order
par(mar=c(2,2,3,6))
image( t(x[o,][nrow(res):1,]) ,axes=F,frame=F , col=c("grey95",colorRampPalette(c("grey95","grey70","firebrick","firebrick4"))(90)))
text( par("usr")[2], seq(nrow(res),0,length.out = nrow(res)) / nrow(res) , labels = rownames(x)[o] , pos=4 ,xpd=T , cex=.4)
text( par("usr")[1], par("usr")[4]-diff(par("usr")[4:3])/30, labels = "origin" , adj= c(0,0),xpd=T , cex=.8 , pos=3)
text( par("usr")[2], par("usr")[4]-diff(par("usr")[4:3])/30, labels = "lineage1" , adj= c(1,0),xpd=T , cex=.8, pos=3)
x <- rev(curves@metadata[["lineages"]][["Lineage1"]])
text( seq(mean(par("usr")[2]),par("usr")[1],length.out = length(x)), par("usr")[4], labels = x , adj= c(0,0),xpd=T , cex=.8 , pos=3)
dev.off()

#SMOOTHED HEATMAP
png(filename = paste0(cell_type, "_trajectory_branch_markers_hierachical.png"),width = 800*4,height = 800*12,res = 300)
x <- t(apply(res,1,function(x){scale(x,T,T)}))
#to_order <- apply(x,1,function(x){ mean( (1:length(x))[order(x,decreasing = T)][1] )})
#to_order <- apply( (res[,1:40]),1,function(x){ mean( (1:length(x))[order(x,decreasing = F)][1] )})
to_order <- apply(x,1,function(x){ which.max(x) })
o <- order(to_order)
o <- hclust(as.dist((1-cor(t(res),use = "pairwise.complete.obs"))/2),method = method )$order
par(mar=c(2,2,3,6))
image( t(x[o,][nrow(res):1,]) ,axes=F,frame=F , col=c("grey95",colorRampPalette(c("grey95","grey70","firebrick","firebrick4"))(90)))
text( par("usr")[2], seq(nrow(res),0,length.out = nrow(res)) / nrow(res) , labels = rownames(x)[o] , pos=4 ,xpd=T , cex=.4)
text( par("usr")[1], par("usr")[4]-diff(par("usr")[4:3])/30, labels = "origin" , adj= c(0,0),xpd=T , cex=.8 , pos=3)
text( par("usr")[2], par("usr")[4]-diff(par("usr")[4:3])/30, labels = "lineage1" , adj= c(1,0),xpd=T , cex=.8, pos=3)
x <- rev(curves@metadata[["lineages"]][["Lineage1"]])
text( seq(mean(par("usr")[2]),par("usr")[1],length.out = length(x)), par("usr")[4], labels = x , adj= c(0,0),xpd=T , cex=.8 , pos=3)
dev.off()

gene_dendro <- hclust(as.dist((1-cor(t(res),use = "pairwise.complete.obs"))/2),method = method )
plot(gene_dendro, cex = 0.1, hang = -1)
k=11
rect.hclust(gene_dendro, k = k, border = "red")

cb <- clusterboot(as.dist((1-cor(t(res),use = "pairwise.complete.obs"))/2), distances=TRUE, B=1000, bootmethod=
                    "boot",clustermethod=disthclustCBI,
                  k=k, cut="number",method=method, showplots=FALSE, seed=1)

write.csv(cb$bootmean, file = paste0(cell_type,"_stability_values_markers.csv"))
gene_dendro1 <- cutree(gene_dendro, k=k)[gene_dendro$order]
write.csv(gene_dendro1, file = paste0(cell_type,"_marker_modules.csv"))

#------
plot(gene_dendro1)

# reorder the leaves of the dendrogram
a <- which(gene_dendro1==c(9))
b <- which(gene_dendro1==c(11))
c <- which(gene_dendro1==c(8))
d <- which(gene_dendro1==c(4))
e <- which(gene_dendro1==c(3))
f <- which(gene_dendro1==c(7))
g <- rev(which(gene_dendro1==c(10)))
h <- which(gene_dendro1==c(6))
i <- which(gene_dendro1==c(5))
j <- which(gene_dendro1==c(1))
k1 <- which(gene_dendro1==c(2))

plot(rotate(gene_dendro, c(a,b, c, d, e, f, g, h, i, j, k1)),cex = 0.1, hang = -1)
rect.hclust(rotate(gene_dendro, c(a,b, c, d, e, f, g, h, i, j, k1)), k = k, border = "red")
a <-rotate(gene_dendro, c(a,b, c, d, e, f, g, h, i, j, k1))

png(filename = paste0(cell_type, "_trajectory_branch_markers_reordered.png"),width = 800*3,height = 800*9,res = 300)
x <- t(apply(res,1,function(x){scale(x,T,T)}))
to_order <- apply(x,1,function(x){ mean( (1:length(x))[order(x,decreasing = T)][1] )})
to_order <- apply( (res[,1:20] - res[,22:41]),1,function(x){ mean( (1:length(x))[order(x,decreasing = F)][1] )})
to_order <- apply(x,1,function(x){ which.max(x) })
o <- a[["order"]]
par(mar=c(2,2,3,6))
image( t(x[o,][nrow(res):1,]) ,axes=F,frame=F , col=c("grey95",colorRampPalette(c("grey95","grey70","firebrick","firebrick4"))(90)))
text( par("usr")[2], seq(nrow(res),0,length.out = nrow(res)) / nrow(res) , labels = rownames(x)[o] , pos=4 ,xpd=T , cex=.6)
text( par("usr")[1], par("usr")[4]-diff(par("usr")[4:3])/30, labels = "origin" , adj= c(0,0),xpd=T , cex=.8 , pos=3)
text( par("usr")[2], par("usr")[4]-diff(par("usr")[4:3])/30, labels = "lineage1" , adj= c(1,0),xpd=T , cex=.8, pos=3)
x <- rev(curves@metadata[["lineages"]][["Lineage1"]])
text( seq(mean(par("usr")[2]),par("usr")[1],length.out = length(x)), par("usr")[4], labels = x , adj= c(0,0),xpd=T , cex=.8 , pos=3)
dev.off()


#================================================
# create a subset of each gene module to have equal number of genes for each of them
patternRes1 <- assoRes
patternRes1$pvalue[is.na(patternRes1$pvalue)] <- 1
patternRes1 <- patternRes1[order(patternRes1$waldStat,decreasing = T),]
patternRes1$FDR <- p.adjust(patternRes1$pvalue,method = "BH")
patternRes1 <- patternRes1[ patternRes1$FDR < 0.0001 , ]
patternRes1 <- patternRes1[ patternRes1$meanLogFC > 1 , ]


gene_dendro2 <- data.frame(gene_dendro1)
gene_dendro2$gene <- row.names(gene_dendro2)
names(gene_dendro2)[1] <- "module"
patternRes1$gene <- row.names(patternRes1)
patternRes1 <- join(patternRes1, gene_dendro2, type="inner", by="gene")

# set the number of wanted genes
n=10
#===============================
cluster_order <- c(9, 11, 8, 4, 3, 7, 10, 6, 5, 1, 2)
d <- patternRes1[patternRes1$module=="9",]
d <- d[order(d$waldStat,decreasing = T),]
  d <- d[1:min(n,nrow(d)),]
gene_list1 <- d

for (i in c(11, 8, 4, 3, 7, 10, 6, 5, 1, 2)) {
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


png(filename = paste0(cell_type, "_trajectory_branch_markers_short1.png"),width = 800*3,height = 800*9,res = 300)
x <- t(apply(res1,1,function(x){scale(x,T,T)}))
o <- c(1:nrow(gene_list2))
par(mar=c(2,2,3,6))
image( t(x[o,][nrow(res1):1,]) ,axes=F,frame=F , col=c("grey95",colorRampPalette(c("grey95","grey70","firebrick","firebrick4"))(90)))
text( par("usr")[2], seq(nrow(res1),0,length.out = nrow(res1)) / nrow(res1) , labels = rownames(x)[o] , pos=4 ,xpd=T , cex=.6)
text( par("usr")[1], par("usr")[4]-diff(par("usr")[4:3])/30, labels = "origin" , adj= c(0,0),xpd=T , cex=.8 , pos=3)
text( par("usr")[2], par("usr")[4]-diff(par("usr")[4:3])/30, labels = "lineage1" , adj= c(1,0),xpd=T , cex=.8, pos=3)
x <- rev(curves@metadata[["lineages"]][["Lineage1"]])
text( seq(mean(par("usr")[2]),par("usr")[1],length.out = length(x)), par("usr")[4], labels = x , adj= c(0,0),xpd=T , cex=.8 , pos=3)
dev.off()

write.csv(gene_list1, file = paste0("top_", n,"module_markers.csv"))

save.image(file="neuronal_trajectory_analysis_data.RData")

#=====================================================================
#================================================
# create a subset of each gene module to have equal number of TFs
patternRes1 <- assoRes
patternRes1$pvalue[is.na(patternRes1$pvalue)] <- 1
patternRes1 <- patternRes1[order(patternRes1$waldStat,decreasing = T),]
patternRes1$FDR <- p.adjust(patternRes1$pvalue,method = "BH")
patternRes1 <- patternRes1[ patternRes1$FDR < 0.0001 , ]
patternRes1 <- patternRes1[ patternRes1$meanLogFC > 1 , ]


gene_dendro2 <- data.frame(gene_dendro1)
gene_dendro2$gene <- row.names(gene_dendro2)
names(gene_dendro2)[1] <- "module"
patternRes1$gene <- row.names(patternRes1)
patternRes1 <- join(patternRes1, gene_dendro2, type="inner", by="gene")

# select only the TFs
TFs <- read.csv("TFs.csv")
patternRes1 <- patternRes1[patternRes1$gene %in% TFs$gene,]


# set the number of wanted TFs/module
n=20
#===============================
cluster_order <- c(9, 11, 8, 4, 3, 7, 10, 6, 5, 1, 2)
d <- patternRes1[patternRes1$module=="9",]
d <- d[order(d$waldStat,decreasing = T),]
d <- d[1:min(n,nrow(d)),]
gene_list1 <- d

for (i in c(11, 8, 4, 3, 7, 10, 6, 5, 1, 2)) {
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


png(filename = paste0(cell_type, "_trajectory_branch_TFs.png"),width = 800*1.5,height = 800*3,res = 300)
x <- t(apply(res1,1,function(x){scale(x,T,T)}))
o <- c(1:nrow(gene_list2))
par(mar=c(2,2,3,6))
image( t(x[o,][nrow(res1):1,]) ,axes=F,frame=F , col=c("grey95",colorRampPalette(c("grey95","grey70","firebrick","firebrick4"))(90)))
text( par("usr")[2], seq(nrow(res1),0,length.out = nrow(res1)) / nrow(res1) , labels = rownames(x)[o] , pos=4 ,xpd=T , cex=1)
text( par("usr")[1], par("usr")[4]-diff(par("usr")[4:3])/30, labels = "SCP" , adj= c(0,0),xpd=T , cex=1.2 , pos=3)
text( par("usr")[2], par("usr")[4]-diff(par("usr")[4:3])/30, labels = "neuronal" , adj= c(1,0),xpd=T , cex=1.2, pos=3)
x <- rev(curves@metadata[["lineages"]][["Lineage1"]])
text( seq(mean(par("usr")[2]),par("usr")[1],length.out = length(x)), par("usr")[4], labels = x , adj= c(0,0),xpd=T , cex=.8 , pos=3)
dev.off()

write.csv(gene_list1, file = paste0("top_", n,"module_TFs.csv"))