setwd()
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
cell_type <- "mes_cl"

load("C:/Users/alex/Desktop/fibro_trajectory_20220615/mes_subset_fibro_clusters.RData")

fibro <- SetIdent(fibro, value="indiv_clusters_numerical")
selected_clusters <- c(16, 5, 4, 9, 10)

temp <- fibro
DefaultAssay(temp) <- "integrated"
rm(fibro)

temp <- RunUMAP(temp, 
                assay = "integrated",
                n.neighbors = 25L, 
                min.dist = 0.3, 
                dims = 1:50, 
                spread = 5,
                metric = "cosine",
                repulsion.strength = 0.5,
                negative.sample.rate = 20,
                n.epochs = NULL,
                seed.use = 42L,
                learning.rate = 15,
                n.components = 2L,
                reduction.name = "umap_traj",
                reduction.key="umaptraj_")
DimPlot(temp, reduction = "umap_traj", cols = "alphabet2", label = TRUE, pt.size = 2)

dimred <- temp@reductions$umap_traj@cell.embeddings
temp@meta.data[["indiv_clusters_numerical"]] <- factor(temp@meta.data[["indiv_clusters_numerical"]])
clustering <- factor(temp@meta.data[["indiv_clusters_numerical"]])

dimvis <- temp@reductions$umap_traj@cell.embeddings


set.seed(1)
lineages <- getLineages(data = dimred,
                        clusterLabels = clustering,
                        start.clus = "4",
                        end.clus=c("10", "16"))


lineages
lineages@elementMetadata@listData[["reducedDim"]]<- dimvis
curves <- getCurves(lineages, thresh = 0.001, stretch = .000001, allow.breaks = T, approx_points = 30)
pseudotime <- slingPseudotime(curves, na = FALSE)
cellWeights <- slingCurveWeights(curves)


n_clusters <-21
colors <- DiscretePalette(n_clusters, palette = "alphabet2")

cl_colors <- data.frame(c(0:20), colors)
colnames(cl_colors) <- c("cluster", "color")

colors3 <- cl_colors[cl_colors$cluster %in% selected_clusters,]


png(filename = paste0(cell_type, "_trajectory_umap.png"), width = 800*2, height = 800*2, res = 300)
#mypar(4,4)
plot_meta(x = temp,red = "umap_traj", feat = "indiv_clusters_numerical",frame=F,label = T, col = colors3$color)
lines(curves@metadata[["curves"]][["Lineage1"]][["s"]], lwd = 5, col = c('black',"firebrick","navy","darkgreen"))
lines(curves@metadata[["curves"]][["Lineage2"]][["s"]], lwd = 5, col = c('black',"firebrick","navy","darkgreen"))
dev.off()

pdf(file = paste0(cell_type, "_trajectory_umap.pdf"), width = 4, height = 4)
#mypar(4,4)
plot_meta(x = temp,red = "umap_traj", feat = "indiv_clusters_numerical",frame=F,label = T, col = colors3$color)
lines(curves@metadata[["curves"]][["Lineage1"]][["s"]], lwd = 5, col = c('black',"firebrick","navy","darkgreen"))
lines(curves@metadata[["curves"]][["Lineage2"]][["s"]], lwd = 5, col = c('black',"firebrick","navy","darkgreen"))
dev.off()

DimPlot(temp, reduction = "umap_traj", cols = colors3$color, label = TRUE, pt.size = 2)

save.image(file="fibro_analysis_data.RData")

# because the datase is big, we randomly select 1000 cells
set.seed(1)
sel2 <- sample(colnames(temp) , size = ncol(temp))
counts <- temp@assays$RNA@counts[Matrix::rowSums(temp@assays$RNA@counts[,sel2] > 1) > 50 ,sel2]
dim(counts)

#===========================================================================
#DIFFERENTIAL GENE EXPRESSION

opt_k <- evaluateK(counts = counts,
                   pseudotime = pseudotime[sel2,],
                   cellWeights = cellWeights[sel2,], k = 3:10,
                   verbose = T,parallel=T)

# create a singecellexperiment object
sce <- fitGAM(counts = counts,
              pseudotime = pseudotime[sel2,],
              cellWeights = cellWeights[sel2,],
              nknots = 8, verbose = T,parallel=T)

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

# select the first 100 genes
gene_list <- rownames(patternRes1)[1:nrow(patternRes1)]
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
text( par("usr")[2], seq(nrow(res),0,length.out = nrow(res)) / nrow(res) , labels = rownames(x)[o] , pos=4 ,xpd=T , cex=.4)
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
k=5
rect.hclust(gene_dendro, k = k, border = "red")

cb_TF <- clusterboot(as.dist((1-cor(t(res),use = "pairwise.complete.obs"))/2), distances=TRUE, B=1000, bootmethod=
                       "boot",clustermethod=disthclustCBI,
                     k=k, cut="number",method=method, showplots=FALSE, seed=1)

write.csv(cb_TF$bootmean, file = paste0(cell_type,"_stability_values_markers1.csv"))
gene_dendro1 <- cutree(gene_dendro, k=k)[gene_dendro$order]
write.csv(gene_dendro1, file = paste0(cell_type,"_marker_modules1.csv"))

save.image(file="fibro_analysis_data.RData")
#===============================================================================
plot(gene_dendro1)
f <- data.frame(gene_dendro1)
f <- unique(f$gene_dendro1)
f
# 5 1 3 2 4
# 1, 5, 4, 2, 3
# reorder the leaves of the dendrogram
a <- which(gene_dendro1==c(1))
b <- rev(which(gene_dendro1==c(5)))
c <- which(gene_dendro1==c(4))
d <- which(gene_dendro1==c(2))
e <- which(gene_dendro1==c(3))
#f <- which(gene_dendro1==c(5))
#g <- which(gene_dendro1==c(4))
#h <- which(gene_dendro1==c(5))
#i <- which(gene_dendro1==c(2))

plot(rotate(gene_dendro, c(a,b, c, d, e)),cex = 0.1, hang = -1)
rect.hclust(rotate(gene_dendro, c(a,b, c, d, e)), k = k, border = "red")
a <-rotate(gene_dendro, c(a,b, c, d, e))

png(filename = paste0(cell_type, "_trajectory_branch_markers_reordered.png"),width = 800*3,height = 800*9,res = 300)
x <- t(apply(res,1,function(x){scale(x,T,T)}))
to_order <- apply(x,1,function(x){ mean( (1:length(x))[order(x,decreasing = T)][1] )})
to_order <- apply( (res[,1:20] - res[,22:41]),1,function(x){ mean( (1:length(x))[order(x,decreasing = F)][1] )})
to_order <- apply(x,1,function(x){ which.max(x) })
o <- a[["order"]]
par(mar=c(2,2,3,6))
image( t(x[o,][nrow(res):1,]) ,axes=F,frame=F , col=c("grey95",colorRampPalette(c("grey95","grey70","firebrick","firebrick4"))(90)))
text( par("usr")[2], seq(nrow(res),0,length.out = nrow(res)) / nrow(res) , labels = rownames(x)[o] , pos=4 ,xpd=T , cex=.4)
text( par("usr")[1], par("usr")[4]-diff(par("usr")[4:3])/30, labels = "lineage 2" , adj= c(0,0),xpd=T , cex=1 , pos=3)
text( par("usr")[2], par("usr")[4]-diff(par("usr")[4:3])/30, labels = "lineage 1" , adj= c(1,0),xpd=T , cex=1, pos=3)
text( mean(par("usr")[1:2]), par("usr")[4]-diff(par("usr")[4:3])/30, labels = "origin" , adj= c(1,0),xpd=T , cex=.8, pos=3)
x <- curves@metadata[["lineages"]][["Lineage1"]]
text( seq(mean(par("usr")[1:2]),par("usr")[1],length.out = length(x)), par("usr")[4], labels = x , adj= c(0,0),xpd=T , cex=.8 , pos=3)
x <- curves@metadata[["lineages"]][["Lineage2"]]
text( seq(mean(par("usr")[1:2]),par("usr")[2],length.out = length(x)), par("usr")[4], labels = x , adj= c(0,0),xpd=T , cex=.8 , pos=3)
dev.off()

pdf(paste0(cell_type, "_trajectory_branch_markers_reordered.pdf"), width = 6,height = 18,paper = "special")
x <- t(apply(res,1,function(x){scale(x,T,T)}))
to_order <- apply(x,1,function(x){ mean( (1:length(x))[order(x,decreasing = T)][1] )})
to_order <- apply( (res[,1:20] - res[,22:41]),1,function(x){ mean( (1:length(x))[order(x,decreasing = F)][1] )})
to_order <- apply(x,1,function(x){ which.max(x) })
o <- a[["order"]]
par(mar=c(2,2,3,6))
image( t(x[o,][nrow(res):1,]) ,axes=F,frame=F , col=c("grey95",colorRampPalette(c("grey95","grey70","firebrick","firebrick4"))(90)))
text( par("usr")[2], seq(nrow(res),0,length.out = nrow(res)) / nrow(res) , labels = rownames(x)[o] , pos=4 ,xpd=T , cex=.3)
text( par("usr")[1], par("usr")[4]-diff(par("usr")[4:3])/30, labels = "lineage 2" , adj= c(0,0),xpd=T , cex=1 , pos=3)
text( par("usr")[2], par("usr")[4]-diff(par("usr")[4:3])/30, labels = "lineage 1" , adj= c(1,0),xpd=T , cex=1, pos=3)
text( mean(par("usr")[1:2]), par("usr")[4]-diff(par("usr")[4:3])/30, labels = "origin" , adj= c(1,0),xpd=T , cex=.8, pos=3)
x <- curves@metadata[["lineages"]][["Lineage1"]]
text( seq(mean(par("usr")[1:2]),par("usr")[1],length.out = length(x)), par("usr")[4], labels = x , adj= c(0,0),xpd=T , cex=.8 , pos=3)
x <- curves@metadata[["lineages"]][["Lineage2"]]
text( seq(mean(par("usr")[1:2]),par("usr")[2],length.out = length(x)), par("usr")[4], labels = x , adj= c(0,0),xpd=T , cex=.8 , pos=3)
dev.off()

save.image(file="fibro_analysis_data.RData")
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

# set the number of wanted genes
n=10
#===============================
cluster_order <- c(1, 5, 4, 2, 3)
d <- patternRes1[patternRes1$module=="1",]
d <- d[order(d$waldStat,decreasing = T),]
d <- d[1:min(n,nrow(d)),]
gene_list1 <- d

for (i in c(5, 4, 2, 3)) {
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


png(filename = paste0(cell_type, "_trajectory_branch_markers_short1.png"),width = 800*2,height = 800*4,res = 300)
x <- t(apply(res1,1,function(x){scale(x,T,T)}))
to_order <- apply(x,1,function(x){ mean( (1:length(x))[order(x,decreasing = T)][1] )})
to_order <- apply( (res1[,1:20] - res1[,22:41]),1,function(x){ mean( (1:length(x))[order(x,decreasing = F)][1] )})
to_order <- apply(x,1,function(x){ which.max(x) })
o <- c(1:nrow(gene_list2))
par(mar=c(2,2,3,6))
image( t(x[o,][nrow(res1):1,]) ,axes=F,frame=F , col=c("grey95",colorRampPalette(c("grey95","grey70","firebrick","firebrick4"))(90)))
text( par("usr")[2], seq(nrow(res1),0,length.out = nrow(res1)) / nrow(res1) , labels = rownames(x)[o] , pos=4 ,xpd=T , cex=1.25)
text( par("usr")[1], par("usr")[4]-diff(par("usr")[4:3])/30, labels = "lineage 2" , adj= c(0,0),xpd=T , cex=1 , pos=3)
text( par("usr")[2], par("usr")[4]-diff(par("usr")[4:3])/30, labels = "lineage 1" , adj= c(1,0),xpd=T , cex=1, pos=3)
text( mean(par("usr")[1:2]), par("usr")[4]-diff(par("usr")[4:3])/30, labels = "origin" , adj= c(1,0),xpd=T , cex=1.25, pos=3)
x <- curves@metadata[["lineages"]][["Lineage1"]]
text( seq(mean(par("usr")[1:2]),par("usr")[1],length.out = length(x)), par("usr")[4], labels = x , adj= c(0,0),xpd=T , cex=1.25 , pos=3)
x <- curves@metadata[["lineages"]][["Lineage2"]]
text( seq(mean(par("usr")[1:2]),par("usr")[2],length.out = length(x)), par("usr")[4], labels = x , adj= c(0,0),xpd=T , cex=1.25 , pos=3)
dev.off()

pdf(paste0(cell_type, "_trajectory_branch_markers_short1.pdf"),width = 4,height = 9,paper="special")
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
n=20
#===============================
cluster_order <- c(1, 5, 4, 2, 3)
d <- patternRes1[patternRes1$module=="1",]
d <- d[order(d$waldStat,decreasing = T),]
d <- d[1:min(n,nrow(d)),]
gene_list1 <- d

for (i in c(5, 4, 2, 3)) {
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


png(filename = paste0(cell_type, "_trajectory_branch_markers_short2.png"),width = 800*2,height = 800*6,res = 300)
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

pdf(paste0(cell_type, "_trajectory_branch_markers_short2.pdf"),width = 4,height = 12,paper="special")
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
# create a subset of each gene module to have equal number of TFs for each of them
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
n=10
#===============================
cluster_order <- c(1, 5, 4, 2, 3)
d <- patternRes1[patternRes1$module=="1",]
d <- d[order(d$waldStat,decreasing = T),]
d <- d[1:min(n,nrow(d)),]
gene_list1 <- d

for (i in c(5, 4, 2, 3)) {
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

pdf(paste0(cell_type, "_trajectory_branch_TFs_short1.pdf"),width = 3,height = 9,paper="special")
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
#================================================
#================================================
# create a subset of each gene module to have equal number of secreted proteins for each of them
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

# select only the secreted proteins
secreted <-  read.delim("C:/Users/alex/Desktop/fibro_trajectory_20220615/protein_class_secreted.tsv")[,1:4]
patternRes1 <- patternRes1[row.names(patternRes1) %in% secreted$gene,]

gene_dendro2 <- data.frame(gene_dendro1)
gene_dendro2$gene <- row.names(gene_dendro2)
names(gene_dendro2)[1] <- "module"
patternRes1$gene <- row.names(patternRes1)
patternRes1 <- join(patternRes1, gene_dendro2, type="inner", by="gene")
patternRes1 <-patternRes1[,1:7]

# set the number of wanted genes
n=10
#===============================
cluster_order <- c(1, 5, 4, 2, 3)
d <- patternRes1[patternRes1$module=="1",]
d <- d[order(d$waldStat,decreasing = T),]
d <- d[1:min(n,nrow(d)),]
gene_list1 <- d

for (i in c(5, 4, 2, 3)) {
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


png(filename = paste0(cell_type, "_trajectory_branch_secreted_short1.png"),width = 800*1.5,height = 800*3,res = 300)
x <- t(apply(res1,1,function(x){scale(x,T,T)}))
to_order <- apply(x,1,function(x){ mean( (1:length(x))[order(x,decreasing = T)][1] )})
to_order <- apply( (res1[,1:20] - res1[,22:41]),1,function(x){ mean( (1:length(x))[order(x,decreasing = F)][1] )})
to_order <- apply(x,1,function(x){ which.max(x) })
o <- c(1:nrow(gene_list2))
par(mar=c(2,2,3,6))
image( t(x[o,][nrow(res1):1,]) ,axes=F,frame=F , col=c("grey95",colorRampPalette(c("grey95","grey70","firebrick","firebrick4"))(90)))
text( par("usr")[2], seq(nrow(res1),0,length.out = nrow(res1)) / nrow(res1) , labels = rownames(x)[o] , pos=4 ,xpd=T , cex=.6)
text( par("usr")[1], par("usr")[4]-diff(par("usr")[4:3])/30, labels = "lineage 2" , adj= c(0,0),xpd=T , cex=1 , pos=3)
text( par("usr")[2], par("usr")[4]-diff(par("usr")[4:3])/30, labels = "lineage 1" , adj= c(1,0),xpd=T , cex=1, pos=3)
text( mean(par("usr")[1:2]), par("usr")[4]-diff(par("usr")[4:3])/30, labels = "origin" , adj= c(1,0),xpd=T , cex=.8, pos=3)
x <- curves@metadata[["lineages"]][["Lineage1"]]
text( seq(mean(par("usr")[1:2]),par("usr")[1],length.out = length(x)), par("usr")[4], labels = x , adj= c(0,0),xpd=T , cex=.8 , pos=3)
x <- curves@metadata[["lineages"]][["Lineage2"]]
text( seq(mean(par("usr")[1:2]),par("usr")[2],length.out = length(x)), par("usr")[4], labels = x , adj= c(0,0),xpd=T , cex=.8 , pos=3)
dev.off()

pdf(paste0(cell_type, "_trajectory_branch_secreted_short1.pdf"),width = 4,height = 9, paper="special")
x <- t(apply(res1,1,function(x){scale(x,T,T)}))
to_order <- apply(x,1,function(x){ mean( (1:length(x))[order(x,decreasing = T)][1] )})
to_order <- apply( (res1[,1:20] - res1[,22:41]),1,function(x){ mean( (1:length(x))[order(x,decreasing = F)][1] )})
to_order <- apply(x,1,function(x){ which.max(x) })
o <- c(1:nrow(gene_list2))
par(mar=c(2,2,3,6))
image( t(x[o,][nrow(res1):1,]) ,axes=F,frame=F , col=c("grey95",colorRampPalette(c("grey95","grey70","firebrick","firebrick4"))(90)))
text( par("usr")[2], seq(nrow(res1),0,length.out = nrow(res1)) / nrow(res1) , labels = rownames(x)[o] , pos=4 ,xpd=T , cex=.6)
text( par("usr")[1], par("usr")[4]-diff(par("usr")[4:3])/30, labels = "lineage 2" , adj= c(0,0),xpd=T , cex=1 , pos=3)
text( par("usr")[2], par("usr")[4]-diff(par("usr")[4:3])/30, labels = "lineage 1" , adj= c(1,0),xpd=T , cex=1, pos=3)
text( mean(par("usr")[1:2]), par("usr")[4]-diff(par("usr")[4:3])/30, labels = "origin" , adj= c(1,0),xpd=T , cex=.8, pos=3)
x <- curves@metadata[["lineages"]][["Lineage1"]]
text( seq(mean(par("usr")[1:2]),par("usr")[1],length.out = length(x)), par("usr")[4], labels = x , adj= c(0,0),xpd=T , cex=.8 , pos=3)
x <- curves@metadata[["lineages"]][["Lineage2"]]
text( seq(mean(par("usr")[1:2]),par("usr")[2],length.out = length(x)), par("usr")[4], labels = x , adj= c(0,0),xpd=T , cex=.8 , pos=3)
dev.off()

write.csv(gene_list1, file = paste0("top_", n,"module_secreted.csv"))

levels(temp) <- c( "0","6", "7", "12", "11")
DefaultAssay(temp) <- "RNA"
DotPlot(temp, assay = "RNA", features = gene_list2$gene, scale = TRUE) + RotatedAxis()

#==============================================================================
#========================================================================
DefaultAssay(temp) <- "RNA"
temp <- NormalizeData(temp, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(temp)
temp <- ScaleData(temp, features = all.genes)

levels(temp) <- c(10, 9, 4, 5, 16)
#colors3 <- colors3[c(1, 2, 5, 3, 4),]

#3x9 inches
VlnPlot(temp, features = c("NKX3-2"), cols = colors3$color, pt.size = 0)


min_cluster <- min(table(temp@active.ident))
lung_cluster.markers <- FindAllMarkers(temp, test.use = "MAST", min.diff.pct = 0.1, max.cells.per.ident = min_cluster, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.1)
# calculate the percentage difference of the positive cells in the indicated cluster, in comparison to the others
lung_cluster.markers$Dpct <- (lung_cluster.markers$pct.1 - lung_cluster.markers$pct.2)

c <-subset(temp, idents =10 )
counts <- as.matrix(c[["RNA"]]@data)
m <- data.frame(rowMeans(counts))
gene <- row.names(c)
m <- cbind(gene, m)
names(m)[2] <- paste0("cluster", "10") 

for (i in c("9","4", "5", "16")) {
  c1 <-subset(temp, idents =i )
  counts1 <- as.matrix(c1[["RNA"]]@data)
  m1 <- data.frame(rowMeans(counts1))
  names(m1)[1] <- paste0("cluster", (i))
  m <- cbind(m, m1)
}
rm(c, counts, counts1, c1, m1)


lung_cluster.markers <-join(lung_cluster.markers, m, by="gene")
write.csv(lung_cluster.markers, file="fibro_markers.csv")


lung_cluster.markers <- lung_cluster.markers[lung_cluster.markers$p_val_adj<0.001, ]
# analysis based of the top markers according to roc power (10 genes) and then avg_log2fc(5 genes)
top2 <- lung_cluster.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top2 <- top2 %>% group_by(cluster) %>% top_n(n = 5, wt = Dpct)
top2 <-top2[, c(7,6)]
top2 <- top2[!duplicated(top2[ , "gene"]),]

levels(temp) <- c("10","9","4", "5", "16")
new_order <- c( "16",  "5", "4", "9", "10")

e <- top2[top2$cluster %in% new_order[1],]
for (i in 2:length(new_order)){
  e1 <- top2[top2$cluster %in% new_order[i],]
  e <- rbind(e, e1)
}

e <- c(e$gene)
e <- e[!duplicated(e)]
#3x9 inches
#fibro_trajectory_top5_markers_dotplot
DotPlot(temp, assay = "RNA", features = e, scale = TRUE) + RotatedAxis()
DotPlot(temp, assay = "RNA", features = e, scale = F) + RotatedAxis()
#========================================================
# plot top5 TFs
TFs <- read.csv("TFs.csv")

# create a new file by merging the TF-list and the enriched genes in each cluster, to identify the TFs, only
lung_cluster_tfs <- join(lung_cluster.markers, TFs, by="gene", type="inner")

# analysis based of the top markers according to avg_log2fc(10 genes) and then the difference in the percent of positives (5 genes)
top3 <- lung_cluster_tfs %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top3 <- top3 %>% group_by(cluster) %>% top_n(n = 5, wt = Dpct)
top3 <-top3[, c(7,6)]
top3 <- top3[!duplicated(top3[ , "gene"]),]

levels(temp) <- c("10","9","4", "5", "16")
new_order <- c( "16",  "5", "4", "9", "10")

f <- top3[top3$cluster %in% new_order[1],]
for (i in 2:length(new_order)){
  f1 <- top3[top3$cluster %in% new_order[i],]
  f <- rbind(f, f1)
}
write.csv(f, file="Ext_Fig4_Fibro_TFs_list.csv")
f <- c(f$gene)
f <- f[!duplicated(f)]
#fibro_trajectory_top5_TFs_dotplot
#3x9 inches 
DotPlot(temp, assay = "RNA", features = f, scale = TRUE) + RotatedAxis()
DotPlot(temp, assay = "RNA", features = f, scale = F) + RotatedAxis()

#========================================================
# plot top5 secr
secr <- data.frame(read.delim("protein_class_secreted.tsv"))
secr <- secr[1:10]

# create a new file by merging the secreted protein-list and the enriched genes in each cluster, to identify only the secreted proteins
lung_cluster_secr <- join(lung_cluster.markers, secr, by="gene", type="inner")

# analysis based of the top markers according to avg_log2fc(10 genes) and then the difference in the percent of positives (5 genes)
top4 <- lung_cluster_secr %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top4 <- top4 %>% group_by(cluster) %>% top_n(n = 5, wt = Dpct)
top4 <-top4[, c(7,6)]
top4 <- top4[!duplicated(top4[ , "gene"]),]

levels(temp) <- c("10","9","4", "5", "16")
new_order <- c( "16",  "5", "4", "9", "10")

h <- top4[top4$cluster %in% new_order[1],]
for (i in 2:length(new_order)){
  h1 <- top4[top4$cluster %in% new_order[i],]
  h <- rbind(h, h1)
}
write.csv(h, file="Ext_Fig4_Fibro_secreted_list.csv")
h <- c(h$gene)
h <- h[!duplicated(h)]
#3x9 inches 
#fibro_trajectory_top5_secreted_dotplot
DotPlot(temp, assay = "RNA", features = h, scale = TRUE) + RotatedAxis()
DotPlot(temp, assay = "RNA", features = h, scale = F) + RotatedAxis()
#========================================================
# plot top5 CDs
cds <- data.frame(read.delim("protein_class_CD.tsv"))
cds <- cds[1:10]

# create a new file by merging the CD-list and the enriched genes in each cluster, to identify only the CDs
lung_cluster_cds <- join(lung_cluster.markers, cds, by="gene", type="inner")

# analysis based of the top markers according to the avg_log2fc(10 genes) and then the difference in the percent of positives (5 genes)
top5 <- lung_cluster_cds %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top5 <- top5 %>% group_by(cluster) %>% top_n(n = 5, wt = Dpct)
top5 <-top5[, c(7,6)]
top5 <- top5[!duplicated(top5[ , "gene"]),]

levels(temp) <- c("10","9","4", "5", "16")
new_order <- c( "16",  "5", "4", "9", "10")

g <- top5[top5$cluster %in% new_order[1],]
for (i in 2:length(new_order)){
  g1 <- top5[top5$cluster %in% new_order[i],]
  g <- rbind(g, g1)
}

g <- c(g$gene)
g <- g[!duplicated(g)]
#3x9 inches 
DotPlot(temp, assay = "RNA", features = g, scale = TRUE) + RotatedAxis()
DotPlot(temp, assay = "RNA", features = g, scale = F) + RotatedAxis()
#========================================================
# plot top5 Kinases
kinases <- data.frame(read.delim("protein_class_kinases.tsv"))
kinases <- kinases[1:10]

# create a new file by merging the kinase-list and the enriched genes in each cluster, to identify only the kinases
lung_cluster_kinases <- join(lung_cluster.markers, kinases, by="gene", type="inner")

# analysis based of the top markers according to avg_log2fc(10 genes) and then the difference in the percent of positives (5 genes)
top6 <- lung_cluster_kinases %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top6 <- top6 %>% group_by(cluster) %>% top_n(n = 5, wt = Dpct)
top6 <-top6[, c(7,6)]
top6 <- top6[!duplicated(top6[ , "gene"]),]

levels(temp) <- c("10","9","4", "5", "16")
new_order <- c( "16",  "5", "4", "9", "10")

g <- top6[top6$cluster %in% new_order[1],]
for (i in 2:length(new_order)){
  g1 <- top6[top6$cluster %in% new_order[i],]
  g <- rbind(g, g1)
}

g <- c(g$gene)
g <- g[!duplicated(g)]
#3x9 inches 
DotPlot(temp, assay = "RNA", features = g, scale = TRUE) + RotatedAxis()
DotPlot(temp, assay = "RNA", features = g, scale = F) + RotatedAxis()
#========================================================
# plot top5 coTFs
coTFs <- read.csv("TF_cofactors.csv")

# create a new file by merging the coTF-list and the enriched genes in each cluster, to identify the coTFs, only
lung_cluster_cotfs <- join(lung_cluster.markers, coTFs, by="gene", type="inner")

# analysis based of the top markers according to avg_log2fc(10 genes) and then the difference in the percent of positives (5 genes)
top7 <- lung_cluster_cotfs %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top7 <- top7 %>% group_by(cluster) %>% top_n(n = 5, wt = Dpct)
top7 <-top7[, c(7,6)]
top7 <- top7[!duplicated(top7[ , "gene"]),]

levels(temp) <- c("10","9","4", "5", "16")
new_order <- c( "16",  "5", "4", "9", "10")

f <- top7[top7$cluster %in% new_order[1],]
for (i in 2:length(new_order)){
  f1 <- top7[top7$cluster %in% new_order[i],]
  f <- rbind(f, f1)
}

f <- c(f$gene)
f <- f[!duplicated(f)]
#3x9 inches 
DotPlot(temp, assay = "RNA", features = f, scale = TRUE) + RotatedAxis()
DotPlot(temp, assay = "RNA", features = f, scale = F) + RotatedAxis()

#=================================================
# mes_cl16 GO:0031175	neuron projection development genes
g <- c("EPHB2", "TNIK", "BOC", "DOK5", "NCALD", "RPS6KA5", "MAP2", "ATP8A2", "MATN2", "SERPINF1", "NBL1", "NCAM2", 
       "MGLL", "COL25A1", "TENM2", "GFRA2", "BMP4", "BMP5", "SERPINI1", "SEMA6D", "NEGR1", "NRK", "ARHGAP44", 
       "PCDH15", "SYT1", "LRRC7", "FGFR2", "FLRT2", "NOTCH2", "COL4A5", "COL4A6", "TNC", "NTF3", "PPFIA2", 
       "CCDC80", "ALKAL2", "TENM4", "EFNB2", "TENM3", "DNM3", "TPBG", "SEMA3E"
)
mes_cl16_markers <- lung_cluster.markers[lung_cluster.markers$cluster=="16",]
mes_cl16_markers <- mes_cl16_markers[order(mes_cl16_markers$avg_log2FC, decreasing = T),]
g <- mes_cl16_markers[mes_cl16_markers$gene %in% g,] 
g <- g$gene

levels(temp) <- c("10","9","4", "5", "16")

# mes_cl16_GO0031175_genes_dotplot 3x12
DotPlot(temp, assay = "RNA", features = g, scale = TRUE) + RotatedAxis()
DotPlot(temp, assay = "RNA", features = g, scale = F) + RotatedAxis()



save.image(file="fibro_analysis_data.RData")
