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

selected_clusters <- c(0, 2,6, 8, 12, 13)

DefaultAssay(lung) <- "integrated"

sel <- lung$seurat_clusters %in% selected_clusters

temp <- lung[,sel]
levels(temp) <- c(0, 2,6, 8, 12, 13)

save.image(file="arw_sm_org_dataset.RData")

temp <- RunUMAP(temp, 
                assay = "integrated",
                n.neighbors = 10L, 
                min.dist = 0.3, 
                dims = 1:50, 
                spread = 0.5,
                metric = "cosine",
                repulsion.strength = 0.5,
                negative.sample.rate = 5,
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
                        start.clus = "0",
                        end.clus=c("13"))


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

save.image(file="sm_trajectory_no_prol_analysis_data.RData")

DefaultAssay(lung) <- "RNA"
# because the datase is big, we randomly select 1000 cells
set.seed(1)
sel2 <- sample( colnames(temp), size = min(1000,ncol(temp)))
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
              nknots = 7, verbose = T,parallel=T)

# run analysis for gene identification for one trajecroy only
assoRes <- associationTest(sce)

# duplicate the result 
patternRes1 <- assoRes
method <- "ward.D2"
#method <- "complete"

save.image(file="sm_trajectory_no_prol_analysis_data.RData")
#----------------------------------------------------------------------

# filter the results
patternRes1$pvalue[is.na(patternRes1$pvalue)] <- 1
patternRes1 <- patternRes1[order(patternRes1$waldStat,decreasing = T),]
patternRes1$FDR <- p.adjust(patternRes1$pvalue,method = "BH")
patternRes1 <- patternRes1[ patternRes1$FDR < 0.0001 , ]
head(patternRes1,20)

# select the first 100 genes
gene_list <- rownames(patternRes1)[1:min(100,nrow(patternRes1))]
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
png(filename = paste0(cell_type, "_trajectory_branch_markers.png"),width = 800*2,height = 800*6,res = 300)
x <- t(apply(res,1,function(x){scale(x,T,T)}))
#to_order <- apply(x,1,function(x){ mean( (1:length(x))[order(x,decreasing = T)][1] )})
#to_order <- apply( (res[,1:40]),1,function(x){ mean( (1:length(x))[order(x,decreasing = F)][1] )})
to_order <- apply(x,1,function(x){ which.max(x) })
o <- order(to_order)
#o <- hclust(as.dist((1-cor(t(res),use = "pairwise.complete.obs"))/2),method = method )$order
par(mar=c(2,2,3,6))
image( t(x[o,][nrow(res):1,]) ,axes=F,frame=F , col=c("grey95",colorRampPalette(c("grey95","grey70","firebrick","firebrick4"))(90)))
text( par("usr")[2], seq(nrow(res),0,length.out = nrow(res)) / nrow(res) , labels = rownames(x)[o] , pos=4 ,xpd=T , cex=.8)
text( par("usr")[1], par("usr")[4]-diff(par("usr")[4:3])/30, labels = "origin" , adj= c(0,0),xpd=T , cex=.8 , pos=3)
text( par("usr")[2], par("usr")[4]-diff(par("usr")[4:3])/30, labels = "lineage1" , adj= c(1,0),xpd=T , cex=.8, pos=3)
x <- rev(curves@metadata[["lineages"]][["Lineage1"]])
text( seq(mean(par("usr")[2]),par("usr")[1],length.out = length(x)), par("usr")[4], labels = x , adj= c(0,0),xpd=T , cex=.8 , pos=3)
dev.off()

#SMOOTHED HEATMAP
png(filename = paste0(cell_type, "_trajectory_branch_markers_hierachical.png"),width = 800*2,height = 800*6,res = 300)
x <- t(apply(res,1,function(x){scale(x,T,T)}))
#to_order <- apply(x,1,function(x){ mean( (1:length(x))[order(x,decreasing = T)][1] )})
#to_order <- apply( (res[,1:40]),1,function(x){ mean( (1:length(x))[order(x,decreasing = F)][1] )})
to_order <- apply(x,1,function(x){ which.max(x) })
o <- order(to_order)
o <- hclust(as.dist((1-cor(t(res),use = "pairwise.complete.obs"))/2),method = method )$order
par(mar=c(2,2,3,6))
image( t(x[o,][nrow(res):1,]) ,axes=F,frame=F , col=c("grey95",colorRampPalette(c("grey95","grey70","firebrick","firebrick4"))(90)))
text( par("usr")[2], seq(nrow(res),0,length.out = nrow(res)) / nrow(res) , labels = rownames(x)[o] , pos=4 ,xpd=T , cex=.8)
text( par("usr")[1], par("usr")[4]-diff(par("usr")[4:3])/30, labels = "origin" , adj= c(0,0),xpd=T , cex=.8 , pos=3)
text( par("usr")[2], par("usr")[4]-diff(par("usr")[4:3])/30, labels = "lineage1" , adj= c(1,0),xpd=T , cex=.8, pos=3)
x <- rev(curves@metadata[["lineages"]][["Lineage1"]])
text( seq(mean(par("usr")[2]),par("usr")[1],length.out = length(x)), par("usr")[4], labels = x , adj= c(0,0),xpd=T , cex=.8 , pos=3)
dev.off()

gene_dendro <- hclust(as.dist((1-cor(t(res),use = "pairwise.complete.obs"))/2),method = method )
plot(gene_dendro, cex = 0.5, hang = -1)
k=9
rect.hclust(gene_dendro, k = k, border = "red")

cb <- clusterboot(as.dist((1-cor(t(res),use = "pairwise.complete.obs"))/2), distances=TRUE, B=1000, bootmethod=
                       "boot",clustermethod=disthclustCBI,
                     k=k, cut="number",method=method, showplots=FALSE, seed=1)

write.csv(cb$bootmean, file = paste0(cell_type,"_stability_values_markers.csv"))
gene_dendro1 <- cutree(gene_dendro, k=k)[gene_dendro$order]
write.csv(gene_dendro1, file = paste0(cell_type,"_marker_modules.csv"))

#------
plot(gene_dendro1)
f <- data.frame(gene_dendro1)
f <- unique(f$gene_dendro1)
f

# reorder the leaves of the dendrogram
a <- which(gene_dendro1==4)
b <- rev(which(gene_dendro1==c(6)))
c <- which(gene_dendro1==c(8))
d <- which(gene_dendro1==c(9))
e <- which(gene_dendro1==c(7))
f <- rev(which(gene_dendro1==c(3)))
g <- rev(which(gene_dendro1==c(1)))
h <- which(gene_dendro1==c(5))
i <- which(gene_dendro1==c(2))

plot(rotate(gene_dendro, c(a,b, c, d, e, f, g, h, i)),cex = 0.5, hang = -1)
rect.hclust(rotate(gene_dendro, c(a,b, c, d, e, f, g, h, i)), k = k, border = "red")
a <-rotate(gene_dendro, c(a,b, c, d, e, f, g, h, i))


reordered_dendro <- data.frame(cbind(c(a[["labels"]]), (a[["order"]]), gene_dendro1))
names(reordered_dendro) <- c("gene", "order", "module")
reordered_dendro$order <- as.numeric(reordered_dendro$order)
reordered_dendro <- reordered_dendro[order(reordered_dendro$order),] 
  
png(filename = paste0(cell_type, "_trajectory_branch_markers_reordered.png"),width = 800*2,height = 800*6,res = 300)
x <- t(apply(res,1,function(x){scale(x,T,T)}))
to_order <- apply(x,1,function(x){ mean( (1:length(x))[order(x,decreasing = T)][1] )})
to_order <- apply( (res[,1:20] - res[,22:41]),1,function(x){ mean( (1:length(x))[order(x,decreasing = F)][1] )})
to_order <- apply(x,1,function(x){ which.max(x) })
o <- a[["order"]]
par(mar=c(2,2,3,6))
image( t(x[o,][nrow(res):1,]) ,axes=F,frame=F , col=c("grey95",colorRampPalette(c("grey95","grey70","firebrick","firebrick4"))(90)))
text( par("usr")[2], seq(nrow(res),0,length.out = nrow(res)) / nrow(res) , labels = rownames(x)[o] , pos=4 ,xpd=T , cex=.8)
text( par("usr")[1], par("usr")[4]-diff(par("usr")[4:3])/30, labels = "origin" , adj= c(0,0),xpd=T , cex=.8 , pos=3)
text( par("usr")[2], par("usr")[4]-diff(par("usr")[4:3])/30, labels = "lineage1" , adj= c(1,0),xpd=T , cex=.8, pos=3)
x <- rev(curves@metadata[["lineages"]][["Lineage1"]])
text( seq(mean(par("usr")[2]),par("usr")[1],length.out = length(x)), par("usr")[4], labels = x , adj= c(0,0),xpd=T , cex=.8 , pos=3)
dev.off()



#-----------------------------------------------------------
# select only the TFs
TFs <- read.csv("TFs.csv")
patternRes1 <- assoRes
patternRes1 <- patternRes1[row.names(patternRes1) %in% TFs$gene,]

# filter the results
patternRes1$pvalue[is.na(patternRes1$pvalue)] <- 1
patternRes1 <- patternRes1[order(patternRes1$waldStat,decreasing = T),]
patternRes1$FDR <- p.adjust(patternRes1$pvalue,method = "BH")
patternRes1 <- patternRes1[ patternRes1$FDR < 0.001 , ]
head(patternRes1,20)

# select the first 100 genes
gene_list <- rownames(patternRes1)[1:min(100,nrow(patternRes1))]
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
png(filename = paste0(cell_type, "_trajectory_branch_TFs.png"),width = 800*2,height = 800*6,res = 300)
x <- t(apply(res,1,function(x){scale(x,T,T)}))
#to_order <- apply(x,1,function(x){ mean( (1:length(x))[order(x,decreasing = T)][1] )})
#to_order <- apply( (res[,1:40]),1,function(x){ mean( (1:length(x))[order(x,decreasing = F)][1] )})
to_order <- apply(x,1,function(x){ which.max(x) })
o <- order(to_order)
#o <- hclust(as.dist((1-cor(t(res),use = "pairwise.complete.obs"))/2),method = method )$order
par(mar=c(2,2,3,6))
image( t(x[o,][nrow(res):1,]) ,axes=F,frame=F , col=c("grey95",colorRampPalette(c("grey95","grey70","firebrick","firebrick4"))(90)))
text( par("usr")[2], seq(nrow(res),0,length.out = nrow(res)) / nrow(res) , labels = rownames(x)[o] , pos=4 ,xpd=T , cex=.8)
text( par("usr")[1], par("usr")[4]-diff(par("usr")[4:3])/30, labels = "origin" , adj= c(0,0),xpd=T , cex=.8 , pos=3)
text( par("usr")[2], par("usr")[4]-diff(par("usr")[4:3])/30, labels = "lineage1" , adj= c(1,0),xpd=T , cex=.8, pos=3)
x <- rev(curves@metadata[["lineages"]][["Lineage1"]])
text( seq(mean(par("usr")[2]),par("usr")[1],length.out = length(x)), par("usr")[4], labels = x , adj= c(0,0),xpd=T , cex=.8 , pos=3)
dev.off()

#SMOOTHED HEATMAP
png(filename = paste0(cell_type, "_trajectory_branch_TFs_hierachical.png"),width = 800*2,height = 800*6,res = 300)
x <- t(apply(res,1,function(x){scale(x,T,T)}))
#to_order <- apply(x,1,function(x){ mean( (1:length(x))[order(x,decreasing = T)][1] )})
#to_order <- apply( (res[,1:40]),1,function(x){ mean( (1:length(x))[order(x,decreasing = F)][1] )})
to_order <- apply(x,1,function(x){ which.max(x) })
o <- order(to_order)
o <- hclust(as.dist((1-cor(t(res),use = "pairwise.complete.obs"))/2),method = method )$order
par(mar=c(2,2,3,6))
image( t(x[o,][nrow(res):1,]) ,axes=F,frame=F , col=c("grey95",colorRampPalette(c("grey95","grey70","firebrick","firebrick4"))(90)))
text( par("usr")[2], seq(nrow(res),0,length.out = nrow(res)) / nrow(res) , labels = rownames(x)[o] , pos=4 ,xpd=T , cex=.8)
text( par("usr")[1], par("usr")[4]-diff(par("usr")[4:3])/30, labels = "origin" , adj= c(0,0),xpd=T , cex=.8 , pos=3)
text( par("usr")[2], par("usr")[4]-diff(par("usr")[4:3])/30, labels = "lineage1" , adj= c(1,0),xpd=T , cex=.8, pos=3)
x <- rev(curves@metadata[["lineages"]][["Lineage1"]])
text( seq(mean(par("usr")[2]),par("usr")[1],length.out = length(x)), par("usr")[4], labels = x , adj= c(0,0),xpd=T , cex=.8 , pos=3)
dev.off()

gene_dendro <- hclust(as.dist((1-cor(t(res),use = "pairwise.complete.obs"))/2),method = method )
plot(gene_dendro, cex = 0.5, hang = -1)
k=11
rect.hclust(gene_dendro, k = k, border = "red")

cb <- clusterboot(as.dist((1-cor(t(res),use = "pairwise.complete.obs"))/2), distances=TRUE, B=1000, bootmethod=
                       "boot",clustermethod=disthclustCBI,
                     k=k, cut="number",method=method, showplots=FALSE, seed=1)

write.csv(cb$bootmean, file = paste0(cell_type,"_stability_values_TF.csv"))
gene_dendro1 <- cutree(gene_dendro, k=k)[gene_dendro$order]
write.csv(gene_dendro1, file = paste0(cell_type,"_TF_modules.csv"))

#-----------------------------------------------------------
# select only the coTFs
coTFs <- read.csv("TF_cofactors.csv")
patternRes1 <- assoRes
patternRes1 <- patternRes1[row.names(patternRes1) %in% coTFs$gene,]

# filter the results
patternRes1$pvalue[is.na(patternRes1$pvalue)] <- 1
patternRes1 <- patternRes1[order(patternRes1$waldStat,decreasing = T),]
patternRes1$FDR <- p.adjust(patternRes1$pvalue,method = "BH")
patternRes1 <- patternRes1[ patternRes1$FDR < 0.001 , ]
head(patternRes1,20)

# select the first 100 genes
gene_list <- rownames(patternRes1)[1:min(100,nrow(patternRes1))]
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
png(filename = paste0(cell_type, "_trajectory_branch_coTFs.png"),width = 800*2,height = 800*6,res = 300)
x <- t(apply(res,1,function(x){scale(x,T,T)}))
#to_order <- apply(x,1,function(x){ mean( (1:length(x))[order(x,decreasing = T)][1] )})
#to_order <- apply( (res[,1:40]),1,function(x){ mean( (1:length(x))[order(x,decreasing = F)][1] )})
to_order <- apply(x,1,function(x){ which.max(x) })
o <- order(to_order)
#o <- hclust(as.dist((1-cor(t(res),use = "pairwise.complete.obs"))/2),method = method )$order
par(mar=c(2,2,3,6))
image( t(x[o,][nrow(res):1,]) ,axes=F,frame=F , col=c("grey95",colorRampPalette(c("grey95","grey70","firebrick","firebrick4"))(90)))
text( par("usr")[2], seq(nrow(res),0,length.out = nrow(res)) / nrow(res) , labels = rownames(x)[o] , pos=4 ,xpd=T , cex=.8)
text( par("usr")[1], par("usr")[4]-diff(par("usr")[4:3])/30, labels = "origin" , adj= c(0,0),xpd=T , cex=.8 , pos=3)
text( par("usr")[2], par("usr")[4]-diff(par("usr")[4:3])/30, labels = "lineage1" , adj= c(1,0),xpd=T , cex=.8, pos=3)
x <- rev(curves@metadata[["lineages"]][["Lineage1"]])
text( seq(mean(par("usr")[2]),par("usr")[1],length.out = length(x)), par("usr")[4], labels = x , adj= c(0,0),xpd=T , cex=.8 , pos=3)
dev.off()

#SMOOTHED HEATMAP
png(filename = paste0(cell_type, "_trajectory_branch_coTFs_hierachical.png"),width = 800*2,height = 800*6,res = 300)
x <- t(apply(res,1,function(x){scale(x,T,T)}))
#to_order <- apply(x,1,function(x){ mean( (1:length(x))[order(x,decreasing = T)][1] )})
#to_order <- apply( (res[,1:40]),1,function(x){ mean( (1:length(x))[order(x,decreasing = F)][1] )})
to_order <- apply(x,1,function(x){ which.max(x) })
o <- order(to_order)
o <- hclust(as.dist((1-cor(t(res),use = "pairwise.complete.obs"))/2),method = method )$order
par(mar=c(2,2,3,6))
image( t(x[o,][nrow(res):1,]) ,axes=F,frame=F , col=c("grey95",colorRampPalette(c("grey95","grey70","firebrick","firebrick4"))(90)))
text( par("usr")[2], seq(nrow(res),0,length.out = nrow(res)) / nrow(res) , labels = rownames(x)[o] , pos=4 ,xpd=T , cex=.8)
text( par("usr")[1], par("usr")[4]-diff(par("usr")[4:3])/30, labels = "origin" , adj= c(0,0),xpd=T , cex=.8 , pos=3)
text( par("usr")[2], par("usr")[4]-diff(par("usr")[4:3])/30, labels = "lineage1" , adj= c(1,0),xpd=T , cex=.8, pos=3)
x <- rev(curves@metadata[["lineages"]][["Lineage1"]])
text( seq(mean(par("usr")[2]),par("usr")[1],length.out = length(x)), par("usr")[4], labels = x , adj= c(0,0),xpd=T , cex=.8 , pos=3)
dev.off()

gene_dendro <- hclust(as.dist((1-cor(t(res),use = "pairwise.complete.obs"))/2),method = method )
plot(gene_dendro, cex = 0.5, hang = -1)
k=8
rect.hclust(gene_dendro, k = k, border = "red")

cb <- clusterboot(as.dist((1-cor(t(res),use = "pairwise.complete.obs"))/2), distances=TRUE, B=1000, bootmethod=
                    "boot",clustermethod=disthclustCBI,
                  k=k, cut="number",method=method, showplots=FALSE, seed=1)

write.csv(cb$bootmean, file = paste0(cell_type,"_stability_values_coTF.csv"))
gene_dendro1 <- cutree(gene_dendro, k=k)[gene_dendro$order]
write.csv(gene_dendro1, file = paste0(cell_type,"_coTF_modules.csv"))

#-----------------------------------------------------------
# select only the secreted proteins
# import the list of secretory markers
secr <- data.frame(read.delim("protein_class_secreted.tsv"))
secr <- secr[1:10]
patternRes1 <- assoRes
patternRes1 <- patternRes1[row.names(patternRes1) %in% secr$gene,]


# filter the results
patternRes1$pvalue[is.na(patternRes1$pvalue)] <- 1
patternRes1 <- patternRes1[order(patternRes1$waldStat,decreasing = T),]
patternRes1$FDR <- p.adjust(patternRes1$pvalue,method = "BH")
patternRes1 <- patternRes1[ patternRes1$FDR < 0.001 , ]
head(patternRes1,20)

# select the first 100 genes
gene_list <- rownames(patternRes1)[1:min(100,nrow(patternRes1))]
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
png(filename = paste0(cell_type, "_trajectory_branch_secr.png"),width = 800*2,height = 800*6,res = 300)
x <- t(apply(res,1,function(x){scale(x,T,T)}))
#to_order <- apply(x,1,function(x){ mean( (1:length(x))[order(x,decreasing = T)][1] )})
#to_order <- apply( (res[,1:40]),1,function(x){ mean( (1:length(x))[order(x,decreasing = F)][1] )})
to_order <- apply(x,1,function(x){ which.max(x) })
o <- order(to_order)
#o <- hclust(as.dist((1-cor(t(res),use = "pairwise.complete.obs"))/2),method = method )$order
par(mar=c(2,2,3,6))
image( t(x[o,][nrow(res):1,]) ,axes=F,frame=F , col=c("grey95",colorRampPalette(c("grey95","grey70","firebrick","firebrick4"))(90)))
text( par("usr")[2], seq(nrow(res),0,length.out = nrow(res)) / nrow(res) , labels = rownames(x)[o] , pos=4 ,xpd=T , cex=.8)
text( par("usr")[1], par("usr")[4]-diff(par("usr")[4:3])/30, labels = "origin" , adj= c(0,0),xpd=T , cex=.8 , pos=3)
text( par("usr")[2], par("usr")[4]-diff(par("usr")[4:3])/30, labels = "lineage1" , adj= c(1,0),xpd=T , cex=.8, pos=3)
x <- rev(curves@metadata[["lineages"]][["Lineage1"]])
text( seq(mean(par("usr")[2]),par("usr")[1],length.out = length(x)), par("usr")[4], labels = x , adj= c(0,0),xpd=T , cex=.8 , pos=3)
dev.off()

#SMOOTHED HEATMAP
png(filename = paste0(cell_type, "_trajectory_branch_secr_hierachical.png"),width = 800*2,height = 800*6,res = 300)
x <- t(apply(res,1,function(x){scale(x,T,T)}))
#to_order <- apply(x,1,function(x){ mean( (1:length(x))[order(x,decreasing = T)][1] )})
#to_order <- apply( (res[,1:40]),1,function(x){ mean( (1:length(x))[order(x,decreasing = F)][1] )})
to_order <- apply(x,1,function(x){ which.max(x) })
o <- order(to_order)
o <- hclust(as.dist((1-cor(t(res),use = "pairwise.complete.obs"))/2),method = method )$order
par(mar=c(2,2,3,6))
image( t(x[o,][nrow(res):1,]) ,axes=F,frame=F , col=c("grey95",colorRampPalette(c("grey95","grey70","firebrick","firebrick4"))(90)))
text( par("usr")[2], seq(nrow(res),0,length.out = nrow(res)) / nrow(res) , labels = rownames(x)[o] , pos=4 ,xpd=T , cex=.8)
text( par("usr")[1], par("usr")[4]-diff(par("usr")[4:3])/30, labels = "origin" , adj= c(0,0),xpd=T , cex=.8 , pos=3)
text( par("usr")[2], par("usr")[4]-diff(par("usr")[4:3])/30, labels = "lineage1" , adj= c(1,0),xpd=T , cex=.8, pos=3)
x <- rev(curves@metadata[["lineages"]][["Lineage1"]])
text( seq(mean(par("usr")[2]),par("usr")[1],length.out = length(x)), par("usr")[4], labels = x , adj= c(0,0),xpd=T , cex=.8 , pos=3)
dev.off()

gene_dendro <- hclust(as.dist((1-cor(t(res),use = "pairwise.complete.obs"))/2),method = method )
plot(gene_dendro, cex = 0.5, hang = -1)
k=8
rect.hclust(gene_dendro, k = k, border = "red")

cb <- clusterboot(as.dist((1-cor(t(res),use = "pairwise.complete.obs"))/2), distances=TRUE, B=1000, bootmethod=
                    "boot",clustermethod=disthclustCBI,
                  k=k, cut="number",method=method, showplots=FALSE, seed=1)

write.csv(cb$bootmean, file = paste0(cell_type,"_stability_values_secr.csv"))
gene_dendro1 <- cutree(gene_dendro, k=k)[gene_dendro$order]
write.csv(gene_dendro1, file = paste0(cell_type,"_secr_modules.csv"))

plot(gene_dendro1)

# reorder the leaves of the dendrogram
a <- which(gene_dendro1==6)
b <- which(gene_dendro1==c(2))
c <- which(gene_dendro1==c(5))
d <- which(gene_dendro1==c(7))
e <- which(gene_dendro1==c(4))
f <-rev(which(gene_dendro1==c(1)))
g <- which(gene_dendro1==c(3))
h <- which(gene_dendro1==c(8))

plot(rotate(gene_dendro, c(a,b, c, d, e, f, g, h)),cex = 0.5, hang = -1)
rect.hclust(rotate(gene_dendro, c(a,b, c, d, e, f, g, h)), k = k, border = "red")
a <-rotate(gene_dendro, c(a,b, c, d, e, f, g, h))

png(filename = paste0(cell_type, "_trajectory_branch_secr_reordered.png"),width = 800*2,height = 800*6,res = 300)
x <- t(apply(res,1,function(x){scale(x,T,T)}))
to_order <- apply(x,1,function(x){ mean( (1:length(x))[order(x,decreasing = T)][1] )})
to_order <- apply( (res[,1:20] - res[,22:41]),1,function(x){ mean( (1:length(x))[order(x,decreasing = F)][1] )})
to_order <- apply(x,1,function(x){ which.max(x) })
o <- a[["order"]]
par(mar=c(2,2,3,6))
image( t(x[o,][nrow(res):1,]) ,axes=F,frame=F , col=c("grey95",colorRampPalette(c("grey95","grey70","firebrick","firebrick4"))(90)))
text( par("usr")[2], seq(nrow(res),0,length.out = nrow(res)) / nrow(res) , labels = rownames(x)[o] , pos=4 ,xpd=T , cex=.8)
text( par("usr")[1], par("usr")[4]-diff(par("usr")[4:3])/30, labels = "origin" , adj= c(0,0),xpd=T , cex=.8 , pos=3)
text( par("usr")[2], par("usr")[4]-diff(par("usr")[4:3])/30, labels = "lineage1" , adj= c(1,0),xpd=T , cex=.8, pos=3)
x <- rev(curves@metadata[["lineages"]][["Lineage1"]])
text( seq(mean(par("usr")[2]),par("usr")[1],length.out = length(x)), par("usr")[4], labels = x , adj= c(0,0),xpd=T , cex=.8 , pos=3)
dev.off()

#-----------------------------------------------------------
# select only the CD proteins
# import the list of CD markers
CD <- data.frame(read.delim("protein_class_CD.tsv"))
CD <- CD[1:10]
patternRes1 <- assoRes
patternRes1 <- patternRes1[row.names(patternRes1) %in% CD$gene,]


# filter the results
patternRes1$pvalue[is.na(patternRes1$pvalue)] <- 1
patternRes1 <- patternRes1[order(patternRes1$waldStat,decreasing = T),]
patternRes1$FDR <- p.adjust(patternRes1$pvalue,method = "BH")
patternRes1 <- patternRes1[ patternRes1$FDR < 0.001 , ]
head(patternRes1,20)

# select the first 100 genes
gene_list <- rownames(patternRes1)[1:min(100,nrow(patternRes1))]
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
png(filename = paste0(cell_type, "_trajectory_branch_CD.png"),width = 800*2,height = 800*6,res = 300)
x <- t(apply(res,1,function(x){scale(x,T,T)}))
#to_order <- apply(x,1,function(x){ mean( (1:length(x))[order(x,decreasing = T)][1] )})
#to_order <- apply( (res[,1:40]),1,function(x){ mean( (1:length(x))[order(x,decreasing = F)][1] )})
to_order <- apply(x,1,function(x){ which.max(x) })
o <- order(to_order)
#o <- hclust(as.dist((1-cor(t(res),use = "pairwise.complete.obs"))/2),method = method )$order
par(mar=c(2,2,3,6))
image( t(x[o,][nrow(res):1,]) ,axes=F,frame=F , col=c("grey95",colorRampPalette(c("grey95","grey70","firebrick","firebrick4"))(90)))
text( par("usr")[2], seq(nrow(res),0,length.out = nrow(res)) / nrow(res) , labels = rownames(x)[o] , pos=4 ,xpd=T , cex=.8)
text( par("usr")[1], par("usr")[4]-diff(par("usr")[4:3])/30, labels = "origin" , adj= c(0,0),xpd=T , cex=.8 , pos=3)
text( par("usr")[2], par("usr")[4]-diff(par("usr")[4:3])/30, labels = "lineage1" , adj= c(1,0),xpd=T , cex=.8, pos=3)
x <- rev(curves@metadata[["lineages"]][["Lineage1"]])
text( seq(mean(par("usr")[2]),par("usr")[1],length.out = length(x)), par("usr")[4], labels = x , adj= c(0,0),xpd=T , cex=.8 , pos=3)
dev.off()

#SMOOTHED HEATMAP
png(filename = paste0(cell_type, "_trajectory_branch_CD_hierachical.png"),width = 800*2,height = 800*6,res = 300)
x <- t(apply(res,1,function(x){scale(x,T,T)}))
#to_order <- apply(x,1,function(x){ mean( (1:length(x))[order(x,decreasing = T)][1] )})
#to_order <- apply( (res[,1:40]),1,function(x){ mean( (1:length(x))[order(x,decreasing = F)][1] )})
to_order <- apply(x,1,function(x){ which.max(x) })
o <- order(to_order)
o <- hclust(as.dist((1-cor(t(res),use = "pairwise.complete.obs"))/2),method = method )$order
par(mar=c(2,2,3,6))
image( t(x[o,][nrow(res):1,]) ,axes=F,frame=F , col=c("grey95",colorRampPalette(c("grey95","grey70","firebrick","firebrick4"))(90)))
text( par("usr")[2], seq(nrow(res),0,length.out = nrow(res)) / nrow(res) , labels = rownames(x)[o] , pos=4 ,xpd=T , cex=.8)
text( par("usr")[1], par("usr")[4]-diff(par("usr")[4:3])/30, labels = "origin" , adj= c(0,0),xpd=T , cex=.8 , pos=3)
text( par("usr")[2], par("usr")[4]-diff(par("usr")[4:3])/30, labels = "lineage1" , adj= c(1,0),xpd=T , cex=.8, pos=3)
x <- rev(curves@metadata[["lineages"]][["Lineage1"]])
text( seq(mean(par("usr")[2]),par("usr")[1],length.out = length(x)), par("usr")[4], labels = x , adj= c(0,0),xpd=T , cex=.8 , pos=3)
dev.off()

gene_dendro <- hclust(as.dist((1-cor(t(res),use = "pairwise.complete.obs"))/2),method = method )
plot(gene_dendro, cex = 0.5, hang = -1)
k=7
rect.hclust(gene_dendro, k = k, border = "red")

cb <- clusterboot(as.dist((1-cor(t(res),use = "pairwise.complete.obs"))/2), distances=TRUE, B=1000, bootmethod=
                    "boot",clustermethod=disthclustCBI,
                  k=k, cut="number",method=method, showplots=FALSE, seed=1)

write.csv(cb$bootmean, file = paste0(cell_type,"_stability_values_CD.csv"))
gene_dendro1 <- cutree(gene_dendro, k=k)[gene_dendro$order]
write.csv(gene_dendro1, file = paste0(cell_type,"_CD_modules.csv"))

#-----------------------------------------------------------
# select only the kinase proteins
# import the list of kinases markers
kinases <- data.frame(read.delim("protein_class_kinases.tsv"))
kinases <- kinases[1:10]
patternRes1 <- assoRes
patternRes1 <- patternRes1[row.names(patternRes1) %in% kinases$gene,]


# filter the results
patternRes1$pvalue[is.na(patternRes1$pvalue)] <- 1
patternRes1 <- patternRes1[order(patternRes1$waldStat,decreasing = T),]
patternRes1$FDR <- p.adjust(patternRes1$pvalue,method = "BH")
patternRes1 <- patternRes1[ patternRes1$FDR < 0.001 , ]
head(patternRes1,20)

# select the first 100 genes
gene_list <- rownames(patternRes1)[1:min(100,nrow(patternRes1))]
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
png(filename = paste0(cell_type, "_trajectory_branch_kinases.png"),width = 800*2,height = 800*6,res = 300)
x <- t(apply(res,1,function(x){scale(x,T,T)}))
#to_order <- apply(x,1,function(x){ mean( (1:length(x))[order(x,decreasing = T)][1] )})
#to_order <- apply( (res[,1:40]),1,function(x){ mean( (1:length(x))[order(x,decreasing = F)][1] )})
to_order <- apply(x,1,function(x){ which.max(x) })
o <- order(to_order)
#o <- hclust(as.dist((1-cor(t(res),use = "pairwise.complete.obs"))/2),method = method )$order
par(mar=c(2,2,3,6))
image( t(x[o,][nrow(res):1,]) ,axes=F,frame=F , col=c("grey95",colorRampPalette(c("grey95","grey70","firebrick","firebrick4"))(90)))
text( par("usr")[2], seq(nrow(res),0,length.out = nrow(res)) / nrow(res) , labels = rownames(x)[o] , pos=4 ,xpd=T , cex=.8)
text( par("usr")[1], par("usr")[4]-diff(par("usr")[4:3])/30, labels = "origin" , adj= c(0,0),xpd=T , cex=.8 , pos=3)
text( par("usr")[2], par("usr")[4]-diff(par("usr")[4:3])/30, labels = "lineage1" , adj= c(1,0),xpd=T , cex=.8, pos=3)
x <- rev(curves@metadata[["lineages"]][["Lineage1"]])
text( seq(mean(par("usr")[2]),par("usr")[1],length.out = length(x)), par("usr")[4], labels = x , adj= c(0,0),xpd=T , cex=.8 , pos=3)
dev.off()

#SMOOTHED HEATMAP
png(filename = paste0(cell_type, "_trajectory_branch_kinases_hierachical.png"),width = 800*2,height = 800*6,res = 300)
x <- t(apply(res,1,function(x){scale(x,T,T)}))
#to_order <- apply(x,1,function(x){ mean( (1:length(x))[order(x,decreasing = T)][1] )})
#to_order <- apply( (res[,1:40]),1,function(x){ mean( (1:length(x))[order(x,decreasing = F)][1] )})
to_order <- apply(x,1,function(x){ which.max(x) })
o <- order(to_order)
o <- hclust(as.dist((1-cor(t(res),use = "pairwise.complete.obs"))/2),method = method )$order
par(mar=c(2,2,3,6))
image( t(x[o,][nrow(res):1,]) ,axes=F,frame=F , col=c("grey95",colorRampPalette(c("grey95","grey70","firebrick","firebrick4"))(90)))
text( par("usr")[2], seq(nrow(res),0,length.out = nrow(res)) / nrow(res) , labels = rownames(x)[o] , pos=4 ,xpd=T , cex=.8)
text( par("usr")[1], par("usr")[4]-diff(par("usr")[4:3])/30, labels = "origin" , adj= c(0,0),xpd=T , cex=.8 , pos=3)
text( par("usr")[2], par("usr")[4]-diff(par("usr")[4:3])/30, labels = "lineage1" , adj= c(1,0),xpd=T , cex=.8, pos=3)
x <- rev(curves@metadata[["lineages"]][["Lineage1"]])
text( seq(mean(par("usr")[2]),par("usr")[1],length.out = length(x)), par("usr")[4], labels = x , adj= c(0,0),xpd=T , cex=.8 , pos=3)
dev.off()

gene_dendro <- hclust(as.dist((1-cor(t(res),use = "pairwise.complete.obs"))/2),method = method )
plot(gene_dendro, cex = 0.5, hang = -1)
k=9
rect.hclust(gene_dendro, k = k, border = "red")

cb <- clusterboot(as.dist((1-cor(t(res),use = "pairwise.complete.obs"))/2), distances=TRUE, B=1000, bootmethod=
                    "boot",clustermethod=disthclustCBI,
                  k=k, cut="number",method=method, showplots=FALSE, seed=1)

write.csv(cb$bootmean, file = paste0(cell_type,"_stability_values_kinase.csv"))
gene_dendro1 <- cutree(gene_dendro, k=k)[gene_dendro$order]
write.csv(gene_dendro1, file = paste0(cell_type,"_kinase_modules.csv"))

#========================================================================
DefaultAssay(temp) <- "RNA"
levels(temp) <- c(0, 2,6, 8, 12, 13)

# run differential expression analysis with MAST
min_cluster <- min(table(temp@active.ident))
lung_cluster.markers <- FindAllMarkers(temp, test.use = "MAST", min.diff.pct = 0.1, max.cells.per.ident = min_cluster, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.1)
# calculate the percentage difference of the positive cells in the indicated cluster, in comparison to the others
lung_cluster.markers$Dpct <- (lung_cluster.markers$pct.1 - lung_cluster.markers$pct.2)

c <-subset(temp, idents =0 )
counts <- as.matrix(c[["RNA"]]@data)
m <- data.frame(rowMeans(counts))
gene <- row.names(c)
m <- cbind(gene, m)
names(m)[2] <- paste0("cluster", "0") 

for (i in c("2","6", "8", "12", "13")) {
  c1 <-subset(lung, idents =i )
  counts1 <- as.matrix(c1[["RNA"]]@data)
  m1 <- data.frame(rowMeans(counts1))
  names(m1)[1] <- paste0("cluster", (i))
  m <- cbind(m, m1)
}
rm(c, counts, counts1, c1, m1)


lung_cluster.markers <-join(lung_cluster.markers, m, by="gene")
write.csv(lung_cluster.markers, file="sm_markers.csv")


lung_cluster.markers <- lung_cluster.markers[lung_cluster.markers$p_val_adj<0.001, ]
# analysis based of the top markers according to the avg_log2fc(5 genes) and then the difference in the percent of positives (5 genes)
top2 <- lung_cluster.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top2 <- top2 %>% group_by(cluster) %>% top_n(n = 5, wt = Dpct)
top2 <-top2[, c(7,6)]
top2 <- top2[!duplicated(top2[ , "gene"]),]

levels(temp) <- c("13", "12", "8", "6", "2", "0")
new_order <- c( "0",  "2", "6", "8", "12", "13")

e <- top2[top2$cluster %in% new_order[1],]
for (i in 2:length(new_order)){
  e1 <- top2[top2$cluster %in% new_order[i],]
  e <- rbind(e, e1)
}
write.csv(e, "Ext_Fig2_top5_SM-traj_markers.csv")

e <- c(e$gene)
e <- e[!duplicated(e)]
#3x9 inches
#sm_trajectory_top5_markers_dotplot
DotPlot(temp, assay = "RNA", features = e, scale = TRUE) + RotatedAxis()

pdf("sm_trajectory_top5_markers_dotplot1.pdf", width = 10, height = 3)
#mypar(4,4)
a <- DotPlot(temp, assay = "RNA", features = e, scale = TRUE) + RotatedAxis()
plot(a)
rm(a)
dev.off()

DotPlot(temp, assay = "RNA", features = e, scale = F) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) # + RotatedAxis()
#========================================================
# plot top5 TFs
TFs <- read.csv("TFs.csv")

# create a new file by merging the TF-list and the enriched genes in each cluster, to identify the TFs, only
lung_cluster_tfs <- join(lung_cluster.markers, TFs, by="gene", type="inner")

# analysis based of the top markers according to the avg_log2fc(10 genes) and then the difference in the percent of positives (5 genes)
top3 <- lung_cluster_tfs %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top3 <- top3 %>% group_by(cluster) %>% top_n(n = 5, wt = Dpct)
top3 <-top3[, c(7,6)]
top3 <- top3[!duplicated(top3[ , "gene"]),]

levels(temp) <- c("13", "12", "8", "6", "2", "0")
new_order <- c( "0",  "2", "6", "8", "12", "13")

f <- top3[top3$cluster %in% new_order[1],]
for (i in 2:length(new_order)){
  f1 <- top3[top3$cluster %in% new_order[i],]
  f <- rbind(f, f1)
}
write.csv(f, "Ext_Fig2_top5_SM-traj_TFs.csv")
f <- c(f$gene)
f <- f[!duplicated(f)]
#3x9 inches 
DotPlot(temp, assay = "RNA", features = f, scale = TRUE) + RotatedAxis()

pdf("sm_trajectory_top5_TFs_dotplot1.pdf", width = 10, height = 3.5)
#mypar(4,4)
a <- DotPlot(temp, assay = "RNA", features = f, scale = TRUE) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
plot(a)
rm(a)
dev.off()

#========================================================
# plot top5 secr
secr <- data.frame(read.delim("protein_class_secreted.tsv"))
secr <- secr[1:10]

# create a new file by merging the list of secreted proteins and the enriched genes in each cluster, to identify only the secreted proteins
lung_cluster_secr <- join(lung_cluster.markers, secr, by="gene", type="inner")

# analysis based of the top markers according to the avg_log2fc(5 genes) and then the difference in the percent of positives (5 genes)
top4 <- lung_cluster_secr %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top4 <- top4 %>% group_by(cluster) %>% top_n(n = 5, wt = Dpct)
top4 <-top4[, c(7,6)]
top4 <- top4[!duplicated(top4[ , "gene"]),]

levels(temp) <- c("13", "12", "8", "6", "2", "0")
new_order <- c( "0",  "2", "6", "8", "12", "13")

h <- top4[top4$cluster %in% new_order[1],]
for (i in 2:length(new_order)){
  h1 <- top4[top4$cluster %in% new_order[i],]
  h <- rbind(h, h1)
}
write.csv(h, "Ext_Fig2_top5_SM-traj_Secretory.csv")
h <- c(h$gene)
h <- h[!duplicated(h)]
#3x9 inches 
DotPlot(temp, assay = "RNA", features = h, scale = TRUE) + RotatedAxis()

pdf("sm_trajectory_top5_secreted_dotplot1.pdf", width = 10, height = 3.5)
#mypar(4,4)
a <- DotPlot(temp, assay = "RNA", features = h, scale = TRUE) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
plot(a)
rm(a)
dev.off()

#========================================================
# plot top5 CDs
cds <- data.frame(read.delim("protein_class_CD.tsv"))
cds <- cds[1:10]

# create a new file by merging the CD-list and the enriched genes in each cluster, to identify only the CDs
lung_cluster_cds <- join(lung_cluster.markers, cds, by="gene", type="inner")

# analysis based of the top markers according to the avg_log2fc(5 genes) and then the difference in the percent of positives (5 genes)
top5 <- lung_cluster_cds %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top5 <- top5 %>% group_by(cluster) %>% top_n(n = 5, wt = Dpct)
top5 <-top5[, c(7,6)]
top5 <- top5[!duplicated(top5[ , "gene"]),]

levels(temp) <- c("13", "12", "8", "6", "2", "0")
new_order <- c( "0",  "2", "6", "8", "12", "13")

g <- top5[top5$cluster %in% new_order[1],]
for (i in 2:length(new_order)){
  g1 <- top5[top5$cluster %in% new_order[i],]
  g <- rbind(g, g1)
}
write.csv(g, "Ext_Fig2_top5_SM-traj_CDs.csv")
g <- c(g$gene)
g <- g[!duplicated(g)]
#3x9 inches 
DotPlot(temp, assay = "RNA", features = g, scale = TRUE) + RotatedAxis()

pdf("sm_trajectory_top5_CDs_dotplot1.pdf", width = 8, height = 3.5)
#mypar(4,4)
a <- DotPlot(temp, assay = "RNA", features = g, scale = TRUE) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
plot(a)
rm(a)
dev.off()

#========================================================
# plot top5 Kinases
kinases <- data.frame(read.delim("protein_class_kinases.tsv"))
kinases <- kinases[1:10]

# create a new file by merging the kinase-list and the enriched genes in each cluster, to identify only the kinases
lung_cluster_kinases <- join(lung_cluster.markers, kinases, by="gene", type="inner")

# analysis based of the top markers according to the avg_log2fc(5 genes) and then the difference in the percent of positives (5 genes)
top6 <- lung_cluster_kinases %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top6 <- top6 %>% group_by(cluster) %>% top_n(n = 5, wt = Dpct)
top6 <-top6[, c(7,6)]
top6 <- top6[!duplicated(top6[ , "gene"]),]

levels(temp) <- c("13", "12", "8", "6", "2", "0")
new_order <- c( "0",  "2", "6", "8", "12", "13")

g <- top6[top6$cluster %in% new_order[1],]
for (i in 2:length(new_order)){
  g1 <- top6[top6$cluster %in% new_order[i],]
  g <- rbind(g, g1)
}

g <- c(g$gene)
g <- g[!duplicated(g)]
#3x9 inches 
DotPlot(temp, assay = "RNA", features = g, scale = TRUE) + RotatedAxis()

#========================================================
# plot top5 coTFs
coTFs <- read.csv("TF_cofactors.csv")

# create a new file by merging the coTF-list and the enriched genes in each cluster, to identify only the coTFs
lung_cluster_cotfs <- join(lung_cluster.markers, coTFs, by="gene", type="inner")

# analysis based of the top markers according to the avg_log2fc(5 genes) and then the difference in the percent of positives (5 genes)
top7 <- lung_cluster_cotfs %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top7 <- top7 %>% group_by(cluster) %>% top_n(n = 5, wt = Dpct)
top7 <-top7[, c(7,6)]
top7 <- top7[!duplicated(top7[ , "gene"]),]

levels(temp) <- c("13", "12", "8", "6", "2", "0")
new_order <- c( "0",  "2", "6", "8", "12", "13")

f <- top7[top7$cluster %in% new_order[1],]
for (i in 2:length(new_order)){
  f1 <- top7[top7$cluster %in% new_order[i],]
  f <- rbind(f, f1)
}

f <- c(f$gene)
f <- f[!duplicated(f)]
#3x9 inches 
DotPlot(temp, assay = "RNA", features = f, scale = TRUE) + RotatedAxis()
#==========================================================================
# plot NOTCH genes
pdf("sm_trajectory_NOTCH_genes_dotplot1.pdf", width = 12, height = 3.5)
#mypar(4,4)
a <- DotPlot(temp, assay = "RNA", features = c("NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4", #receptors
                                          "HES1", "HES2", "HES3", "HES5", "HES6", "HES7", "HEY1", "HEY2", "HEYL", "NRARP", #targets
                                          "DLL1", "DLL4", "JAG1", "JAG2", #ligands
                                          "POFUT1", "POFUT2", "LFNG", "MFNG", "RFNG", "MAML1", "MAML2", "MAML3", "MAMLD1", "RBPJ", #transducers
                                          "DLL3", "DLK1", "DLK2", "NUMB", "NUMBL"))+ theme(axis.text.x = element_text(angle = 90, hjust = 1)) # inhibitors
             plot(a)
             rm(a)
             dev.off()
Notch_genes <-c("NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4", #receptors
                "HES1", "HES2", "HES3", "HES5", "HES6", "HES7", "HEY1", "HEY2", "HEYL", "NRARP", #targets
                "DLL1", "DLL4", "JAG1", "JAG2", #ligands
                "POFUT1", "POFUT2", "LFNG", "MFNG", "RFNG", "MAML1", "MAML2", "MAML3", "MAMLD1", "RBPJ", #transducers
                "DLL3", "DLK1", "DLK2", "NUMB", "NUMBL")              
Notch_genes <- data.frame(Notch_genes)
write.csv(Notch_genes, file="Ext_Fig2_H_Notch_genes.csv")
