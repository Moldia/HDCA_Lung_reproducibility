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

selected_clusters <- c(0, 2,6, 8,20, 12, 13)

lung <- lung1
DefaultAssay(lung) <- "integrated"
rm(lung1)
sel <- lung$seurat_clusters %in% selected_clusters

temp <- lung[,sel]
levels(temp) <- c(0, 2,6, 8,20, 12, 13)

save.image(file="arw_sm_org_dataset.RData")

temp <- RunUMAP(temp, 
                assay = "integrated",
                n.neighbors = 10L, 
                min.dist = 0.5, 
                dims = 1:50, 
                spread = 0.5,
                metric = "cosine",
                repulsion.strength = 0.1,
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
                        start.clus = "0",
                        end.clus=c("20", "13"))


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

save.image(file="arw_sm_analysis_final_data.RData")

# because the datase is big, we randomly select 1000 cells
set.seed(1)
sel2 <- sample(colnames(temp) , size = min(1000,ncol(temp)))
counts <- lung@assays$RNA@counts[Matrix::rowSums(lung@assays$RNA@counts[,sel2] > 1) > 50 ,sel2]
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
              nknots = 6, verbose = T,parallel=T)

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
head(patternRes1,20)

# select the first 100 genes
gene_list <- rownames(patternRes1)[1:min(100,nrow(patternRes1))]
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
png(filename = paste0(cell_type, "_trajectory_branch_markers.png"),width = 800*3,height = 800*6,res = 300)
x <- t(apply(res,1,function(x){scale(x,T,T)}))
to_order <- apply(x,1,function(x){ mean( (1:length(x))[order(x,decreasing = T)][1] )})
to_order <- apply( (res[,1:20] - res[,22:41]),1,function(x){ mean( (1:length(x))[order(x,decreasing = F)][1] )})
to_order <- apply(x,1,function(x){ which.max(x) })
o <- order(to_order)
o <- hclust(as.dist((1-cor(t(res),use = "pairwise.complete.obs"))/2),method = method )$order
o <- a[["order"]]
par(mar=c(2,2,3,6))
image( t(x[o,][nrow(res):1,]) ,axes=F,frame=F , col=c("grey95",colorRampPalette(c("grey95","grey70","firebrick","firebrick4"))(90)))
text( par("usr")[2], seq(nrow(res),0,length.out = nrow(res)) / nrow(res) , labels = rownames(x)[o] , pos=4 ,xpd=T , cex=.8)
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
plot(gene_dendro, cex = 0.5, hang = -1)
k=7
rect.hclust(gene_dendro, k = k, border = "red")

cb_TF <- clusterboot(as.dist((1-cor(t(res),use = "pairwise.complete.obs"))/2), distances=TRUE, B=1000, bootmethod=
                       "boot",clustermethod=disthclustCBI,
                     k=k, cut="number",method=method, showplots=FALSE, seed=1)

write.csv(cb_TF$bootmean, file = paste0(cell_type,"_stability_values_markers.csv"))
gene_dendro1 <- cutree(gene_dendro, k=k)[gene_dendro$order]
write.csv(gene_dendro1, file = paste0(cell_type,"_marker_modules.csv"))
plot(gene_dendro1)

# reorder the leaves of the dendrogram
a <- which(gene_dendro1==3)
b <- which(gene_dendro1==c(7))
c <- which(gene_dendro1==c(6))
d <- which(gene_dendro1==c(5))
e <- which(gene_dendro1==c(4))
f <- which(gene_dendro1==c(2))
g <- which(gene_dendro1==c(1))

plot(rotate(gene_dendro, c(a,b, c, d, e, f, g)))
 a <-rotate(gene_dendro, c(a,b, c, d, e, f, g))
 
 png(filename = paste0(cell_type, "_trajectory_branch_markers_reordered.png"),width = 800*3,height = 800*6,res = 300)
 x <- t(apply(res,1,function(x){scale(x,T,T)}))
 to_order <- apply(x,1,function(x){ mean( (1:length(x))[order(x,decreasing = T)][1] )})
 to_order <- apply( (res[,1:20] - res[,22:41]),1,function(x){ mean( (1:length(x))[order(x,decreasing = F)][1] )})
 to_order <- apply(x,1,function(x){ which.max(x) })
 o <- a[["order"]]
 par(mar=c(2,2,3,6))
 image( t(x[o,][nrow(res):1,]) ,axes=F,frame=F , col=c("grey95",colorRampPalette(c("grey95","grey70","firebrick","firebrick4"))(90)))
 text( par("usr")[2], seq(nrow(res),0,length.out = nrow(res)) / nrow(res) , labels = rownames(x)[o] , pos=4 ,xpd=T , cex=.8)
 text( par("usr")[1], par("usr")[4]-diff(par("usr")[4:3])/30, labels = "lineage 2" , adj= c(0,0),xpd=T , cex=1 , pos=3)
 text( par("usr")[2], par("usr")[4]-diff(par("usr")[4:3])/30, labels = "lineage 1" , adj= c(1,0),xpd=T , cex=1, pos=3)
 text( mean(par("usr")[1:2]), par("usr")[4]-diff(par("usr")[4:3])/30, labels = "origin" , adj= c(1,0),xpd=T , cex=.8, pos=3)
 x <- curves@metadata[["lineages"]][["Lineage1"]]
 text( seq(mean(par("usr")[1:2]),par("usr")[1],length.out = length(x)), par("usr")[4], labels = x , adj= c(0,0),xpd=T , cex=.8 , pos=3)
 x <- curves@metadata[["lineages"]][["Lineage2"]]
 text( seq(mean(par("usr")[1:2]),par("usr")[2],length.out = length(x)), par("usr")[4], labels = x , adj= c(0,0),xpd=T , cex=.8 , pos=3)
 dev.off()
 

#===========================================================================
# select only the TFs
TFs <- read.csv("TFs.csv")
patternRes1 <- patternRes1[row.names(patternRes1) %in% TFs$gene,]

# filter the results
patternRes1$pvalue[is.na(patternRes1$pvalue)] <- 1
patternRes1 <- patternRes1[order(patternRes1$waldStat,decreasing = T),]
patternRes1$FDR <- p.adjust(patternRes1$pvalue,method = "BH")
patternRes1 <- patternRes1[ patternRes1$FDR < 0.0001 , ]
head(patternRes1,20)

# select no more than 100 TFs
gene_list <- rownames(patternRes1)[1:min(100,nrow(patternRes1))]
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
png(filename = paste0(cell_type, "_trajectory_branch_TFs.png"),width = 800*4,height = 800*6,res = 300)
x <- t(apply(res,1,function(x){scale(x,T,T)}))
to_order <- apply(x,1,function(x){ mean( (1:length(x))[order(x,decreasing = T)][1] )})
to_order <- apply( (res[,1:20] - res[,22:41]),1,function(x){ mean( (1:length(x))[order(x,decreasing = F)][1] )})
to_order <- apply(x,1,function(x){ which.max(x) })
o <- order(to_order)
o <- hclust(as.dist((1-cor(t(res),use = "pairwise.complete.obs"))/2),method = method )$order
par(mar=c(2,2,3,6))
image( t(x[o,][nrow(res):1,]) ,axes=F,frame=F , col=c("grey95",colorRampPalette(c("grey95","grey70","firebrick","firebrick4"))(90)))
text( par("usr")[2], seq(nrow(res),0,length.out = nrow(res)) / nrow(res) , labels = rownames(x)[o] , pos=4 ,xpd=T , cex=.8)
text( par("usr")[1], par("usr")[4]-diff(par("usr")[4:3])/30, labels = "lineage 2" , adj= c(0,0),xpd=T , cex=.8 , pos=3)
text( par("usr")[2], par("usr")[4]-diff(par("usr")[4:3])/30, labels = "lineage 1" , adj= c(1,0),xpd=T , cex=.8, pos=3)
text( mean(par("usr")[1:2]), par("usr")[4]-diff(par("usr")[4:3])/30, labels = "origin" , adj= c(1,0),xpd=T , cex=.8, pos=3)
x <- curves@metadata[["lineages"]][["Lineage1"]]
text( seq(mean(par("usr")[1:2]),par("usr")[1],length.out = length(x)), par("usr")[4], labels = x , adj= c(0,0),xpd=T , cex=.8 , pos=3)
x <- curves@metadata[["lineages"]][["Lineage2"]]
text( seq(mean(par("usr")[1:2]),par("usr")[2],length.out = length(x)), par("usr")[4], labels = x , adj= c(0,0),xpd=T , cex=.8 , pos=3)
dev.off()

#set number of gene modules
gene_dendro <- hclust(as.dist((1-cor(t(res),use = "pairwise.complete.obs"))/2),method = method )
plot(gene_dendro, cex = 0.5, hang = -1)
k=9
rect.hclust(gene_dendro, k = k, border = "red")

cb_TF <- clusterboot(as.dist((1-cor(t(res),use = "pairwise.complete.obs"))/2), distances=TRUE, B=1000, bootmethod=
                       "boot",clustermethod=disthclustCBI,
                     k=k, cut="number",method=method, showplots=FALSE, seed=1)

write.csv(cb_TF$bootmean, file = paste0(cell_type,"_stability_values_TF.csv"))
gene_dendro1 <- cutree(gene_dendro, k=k)[gene_dendro$order]
write.csv(gene_dendro1, file = paste0(cell_type,"_TF_modules.csv"))

#-------------------------------------------------------------------------------------------------------
# Find the changed coTFs
patternRes1 <- patternRes
# import the list of coTFs
coTFs <- read.csv("TF_cofactors.csv")
patternRes1 <- patternRes1[row.names(patternRes1) %in% coTFs$gene,]

# filter the results
patternRes1$pvalue[is.na(patternRes1$pvalue)] <- 1
patternRes1 <- patternRes1[order(patternRes1$waldStat,decreasing = T),]
patternRes1$FDR <- p.adjust(patternRes1$pvalue,method = "BH")
patternRes1 <- patternRes1[ patternRes1$FDR < 0.0001 , ]
head(patternRes1,20)

# select no more than 100 coTFs
gene_list <- rownames(patternRes1)[1:min(100,nrow(patternRes1))]
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
png(filename = paste0(cell_type, "_trajectory_branch_coTFs.png"),width = 800*4,height = 800*6,res = 300)
x <- t(apply(res,1,function(x){scale(x,T,T)}))
to_order <- apply(x,1,function(x){ mean( (1:length(x))[order(x,decreasing = T)][1] )})
to_order <- apply( (res[,1:20] - res[,22:41]),1,function(x){ mean( (1:length(x))[order(x,decreasing = F)][1] )})
to_order <- apply(x,1,function(x){ which.max(x) })
o <- order(to_order)
o <- hclust(as.dist((1-cor(t(res),use = "pairwise.complete.obs"))/2),method = method )$order
par(mar=c(2,2,3,6))
image( t(x[o,][nrow(res):1,]) ,axes=F,frame=F , col=c("grey95",colorRampPalette(c("grey95","grey70","firebrick","firebrick4"))(90)))
text( par("usr")[2], seq(nrow(res),0,length.out = nrow(res)) / nrow(res) , labels = rownames(x)[o] , pos=4 ,xpd=T , cex=.8)
text( par("usr")[1], par("usr")[4]-diff(par("usr")[4:3])/30, labels = "lineage 2" , adj= c(0,0),xpd=T , cex=.8 , pos=3)
text( par("usr")[2], par("usr")[4]-diff(par("usr")[4:3])/30, labels = "lineage 1" , adj= c(1,0),xpd=T , cex=.8, pos=3)
text( mean(par("usr")[1:2]), par("usr")[4]-diff(par("usr")[4:3])/30, labels = "origin" , adj= c(1,0),xpd=T , cex=.8, pos=3)
x <- curves@metadata[["lineages"]][["Lineage1"]]
text( seq(mean(par("usr")[1:2]),par("usr")[1],length.out = length(x)), par("usr")[4], labels = x , adj= c(0,0),xpd=T , cex=.8 , pos=3)
x <- curves@metadata[["lineages"]][["Lineage2"]]
text( seq(mean(par("usr")[1:2]),par("usr")[2],length.out = length(x)), par("usr")[4], labels = x , adj= c(0,0),xpd=T , cex=.8 , pos=3)
dev.off()

#set number of gene modules
gene_dendro <- hclust(as.dist((1-cor(t(res),use = "pairwise.complete.obs"))/2),method = method )
plot(gene_dendro, cex = 0.5, hang = -1)
k=7
rect.hclust(gene_dendro, k = k, border = "red")

cb_TF <- clusterboot(as.dist((1-cor(t(res),use = "pairwise.complete.obs"))/2), distances=TRUE, B=1000, bootmethod=
                       "boot",clustermethod=disthclustCBI,
                     k=k, cut="number",method=method, showplots=FALSE, seed=1)

write.csv(cb_TF$bootmean, file = paste0(cell_type,"_stability_values_coTF.csv"))
gene_dendro1 <- cutree(gene_dendro, k=k)[gene_dendro$order]
write.csv(gene_dendro1, file = paste0(cell_type,"_coTF_modules.csv"))

#------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------
# Find the changed kinases
patternRes1 <- patternRes
# import the list of kinases
kinases <- data.frame(read.delim("protein_class_kinases.tsv"))
kinases <- kinases[1:10]
patternRes1 <- patternRes1[row.names(patternRes1) %in% kinases$gene,]

# filter the results
patternRes1$pvalue[is.na(patternRes1$pvalue)] <- 1
patternRes1 <- patternRes1[order(patternRes1$waldStat,decreasing = T),]
patternRes1$FDR <- p.adjust(patternRes1$pvalue,method = "BH")
patternRes1 <- patternRes1[ patternRes1$FDR < 0.0001 , ]
head(patternRes1,20)

# select no more than 100 kinases
gene_list <- rownames(patternRes1)[1:min(100,nrow(patternRes1))]
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
png(filename = paste0(cell_type, "_trajectory_branch_kinasess.png"),width = 800*2,height = 800*3,res = 300)
x <- t(apply(res,1,function(x){scale(x,T,T)}))
to_order <- apply(x,1,function(x){ mean( (1:length(x))[order(x,decreasing = T)][1] )})
to_order <- apply( (res[,1:20] - res[,22:41]),1,function(x){ mean( (1:length(x))[order(x,decreasing = F)][1] )})
to_order <- apply(x,1,function(x){ which.max(x) })
o <- order(to_order)
o <- hclust(as.dist((1-cor(t(res),use = "pairwise.complete.obs"))/2),method = method )$order
par(mar=c(2,2,3,6))
image( t(x[o,][nrow(res):1,]) ,axes=F,frame=F , col=c("grey95",colorRampPalette(c("grey95","grey70","firebrick","firebrick4"))(90)))
text( par("usr")[2], seq(nrow(res),0,length.out = nrow(res)) / nrow(res) , labels = rownames(x)[o] , pos=4 ,xpd=T , cex=.8)
text( par("usr")[1], par("usr")[4]-diff(par("usr")[4:3])/30, labels = "lineage 2" , adj= c(0,0),xpd=T , cex=.8 , pos=3)
text( par("usr")[2], par("usr")[4]-diff(par("usr")[4:3])/30, labels = "lineage 1" , adj= c(1,0),xpd=T , cex=.8, pos=3)
text( mean(par("usr")[1:2]), par("usr")[4]-diff(par("usr")[4:3])/30, labels = "origin" , adj= c(1,0),xpd=T , cex=.8, pos=3)
x <- curves@metadata[["lineages"]][["Lineage1"]]
text( seq(mean(par("usr")[1:2]),par("usr")[1],length.out = length(x)), par("usr")[4], labels = x , adj= c(0,0),xpd=T , cex=.8 , pos=3)
x <- curves@metadata[["lineages"]][["Lineage2"]]
text( seq(mean(par("usr")[1:2]),par("usr")[2],length.out = length(x)), par("usr")[4], labels = x , adj= c(0,0),xpd=T , cex=.8 , pos=3)
dev.off()

#set number of gene modules
gene_dendro <- hclust(as.dist((1-cor(t(res),use = "pairwise.complete.obs"))/2),method = method )
plot(gene_dendro, cex = 0.5, hang = -1)
k=8
rect.hclust(gene_dendro, k = k, border = "red")

cb_TF <- clusterboot(as.dist((1-cor(t(res),use = "pairwise.complete.obs"))/2), distances=TRUE, B=1000, bootmethod=
                       "boot",clustermethod=disthclustCBI,
                     k=k, cut="number",method=method, showplots=FALSE, seed=1)

write.csv(cb_TF$bootmean, file = paste0(cell_type,"_stability_values_kinases.csv"))
gene_dendro1 <- cutree(gene_dendro, k=k)[gene_dendro$order]
write.csv(gene_dendro1, file = paste0(cell_type,"_kinases_modules.csv"))

#--------------------------------------------------------------
# Find the changed CDs
patternRes1 <- patternRes
# import the list of CD markers
cds <- data.frame(read.delim("protein_class_CD.tsv"))
cds <- cds[1:10]
patternRes1 <- patternRes1[row.names(patternRes1) %in% cds$gene,]

# filter the results
patternRes1$pvalue[is.na(patternRes1$pvalue)] <- 1
patternRes1 <- patternRes1[order(patternRes1$waldStat,decreasing = T),]
patternRes1$FDR <- p.adjust(patternRes1$pvalue,method = "BH")
patternRes1 <- patternRes1[ patternRes1$FDR < 0.0001 , ]
head(patternRes1,20)

# select no more than 100 CDs
gene_list <- rownames(patternRes1)[1:min(100,nrow(patternRes1))]
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
png(filename = paste0(cell_type, "_trajectory_branch_CDs.png"),width = 800*2,height = 800*3,res = 300)
x <- t(apply(res,1,function(x){scale(x,T,T)}))
to_order <- apply(x,1,function(x){ mean( (1:length(x))[order(x,decreasing = T)][1] )})
to_order <- apply( (res[,1:20] - res[,22:41]),1,function(x){ mean( (1:length(x))[order(x,decreasing = F)][1] )})
to_order <- apply(x,1,function(x){ which.max(x) })
o <- order(to_order)
o <- hclust(as.dist((1-cor(t(res),use = "pairwise.complete.obs"))/2),method = method )$order
par(mar=c(2,2,3,6))
image( t(x[o,][nrow(res):1,]) ,axes=F,frame=F , col=c("grey95",colorRampPalette(c("grey95","grey70","firebrick","firebrick4"))(90)))
text( par("usr")[2], seq(nrow(res),0,length.out = nrow(res)) / nrow(res) , labels = rownames(x)[o] , pos=4 ,xpd=T , cex=.8)
text( par("usr")[1], par("usr")[4]-diff(par("usr")[4:3])/30, labels = "lineage 2" , adj= c(0,0),xpd=T , cex=.8 , pos=3)
text( par("usr")[2], par("usr")[4]-diff(par("usr")[4:3])/30, labels = "lineage 1" , adj= c(1,0),xpd=T , cex=.8, pos=3)
text( mean(par("usr")[1:2]), par("usr")[4]-diff(par("usr")[4:3])/30, labels = "origin" , adj= c(1,0),xpd=T , cex=.8, pos=3)
x <- curves@metadata[["lineages"]][["Lineage1"]]
text( seq(mean(par("usr")[1:2]),par("usr")[1],length.out = length(x)), par("usr")[4], labels = x , adj= c(0,0),xpd=T , cex=.8 , pos=3)
x <- curves@metadata[["lineages"]][["Lineage2"]]
text( seq(mean(par("usr")[1:2]),par("usr")[2],length.out = length(x)), par("usr")[4], labels = x , adj= c(0,0),xpd=T , cex=.8 , pos=3)
dev.off()

#set number of gene modules
gene_dendro <- hclust(as.dist((1-cor(t(res),use = "pairwise.complete.obs"))/2),method = method )
plot(gene_dendro, cex = 0.5, hang = -1)
k=6
rect.hclust(gene_dendro, k = k, border = "red")

cb_TF <- clusterboot(as.dist((1-cor(t(res),use = "pairwise.complete.obs"))/2), distances=TRUE, B=1000, bootmethod=
                       "boot",clustermethod=disthclustCBI,
                     k=k, cut="number",method=method, showplots=FALSE, seed=1)

write.csv(cb_TF$bootmean, file = paste0(cell_type,"_stability_values_CDs.csv"))
gene_dendro1 <- cutree(gene_dendro, k=k)[gene_dendro$order]
write.csv(gene_dendro1, file = paste0(cell_type,"_CDs_modules.csv"))

#--------------------------------------------------------------
# Find the changed secreted
patternRes1 <- patternRes
# import the list of secretory markers
secr <- data.frame(read.delim("protein_class_secreted.tsv"))
secr <- secr[1:10]
patternRes1 <- patternRes1[row.names(patternRes1) %in% secr$gene,]

# filter the results
patternRes1$pvalue[is.na(patternRes1$pvalue)] <- 1
patternRes1 <- patternRes1[order(patternRes1$waldStat,decreasing = T),]
patternRes1$FDR <- p.adjust(patternRes1$pvalue,method = "BH")
patternRes1 <- patternRes1[ patternRes1$FDR < 0.0001 , ]
head(patternRes1,20)

# select no more than 100 secreted proteins
gene_list <- rownames(patternRes1)[1:min(100,nrow(patternRes1))]
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
png(filename = paste0(cell_type, "_trajectory_branch_secr.png"),width = 800*4,height = 800*6,res = 300)
x <- t(apply(res,1,function(x){scale(x,T,T)}))
to_order <- apply(x,1,function(x){ mean( (1:length(x))[order(x,decreasing = T)][1] )})
to_order <- apply( (res[,1:20] - res[,22:41]),1,function(x){ mean( (1:length(x))[order(x,decreasing = F)][1] )})
to_order <- apply(x,1,function(x){ which.max(x) })
o <- order(to_order)
o <- hclust(as.dist((1-cor(t(res),use = "pairwise.complete.obs"))/2),method = method )$order
par(mar=c(2,2,3,6))
image( t(x[o,][nrow(res):1,]) ,axes=F,frame=F , col=c("grey95",colorRampPalette(c("grey95","grey70","firebrick","firebrick4"))(90)))
text( par("usr")[2], seq(nrow(res),0,length.out = nrow(res)) / nrow(res) , labels = rownames(x)[o] , pos=4 ,xpd=T , cex=.8)
text( par("usr")[1], par("usr")[4]-diff(par("usr")[4:3])/30, labels = "lineage 2" , adj= c(0,0),xpd=T , cex=.8 , pos=3)
text( par("usr")[2], par("usr")[4]-diff(par("usr")[4:3])/30, labels = "lineage 1" , adj= c(1,0),xpd=T , cex=.8, pos=3)
text( mean(par("usr")[1:2]), par("usr")[4]-diff(par("usr")[4:3])/30, labels = "origin" , adj= c(1,0),xpd=T , cex=.8, pos=3)
x <- curves@metadata[["lineages"]][["Lineage1"]]
text( seq(mean(par("usr")[1:2]),par("usr")[1],length.out = length(x)), par("usr")[4], labels = x , adj= c(0,0),xpd=T , cex=.8 , pos=3)
x <- curves@metadata[["lineages"]][["Lineage2"]]
text( seq(mean(par("usr")[1:2]),par("usr")[2],length.out = length(x)), par("usr")[4], labels = x , adj= c(0,0),xpd=T , cex=.8 , pos=3)
dev.off()

#set number of gene modules
gene_dendro <- hclust(as.dist((1-cor(t(res),use = "pairwise.complete.obs"))/2),method = method )
plot(gene_dendro, cex = 0.5, hang = -1)
k=8
rect.hclust(gene_dendro, k = k, border = "red")

cb_secr <- clusterboot(as.dist((1-cor(t(res),use = "pairwise.complete.obs"))/2), distances=TRUE, B=1000, bootmethod=
                         "boot",clustermethod=disthclustCBI,
                       k=k, cut="number",method=method, showplots=FALSE, seed=1)

write.csv(cb_secr$bootmean, file = paste0(cell_type,"_stability_values_secr.csv"))
gene_dendro1 <- cutree(gene_dendro, k=k)[gene_dendro$order]
write.csv(gene_dendro1, file = paste0(cell_type,"_secr_modules.csv"))

#===========================================================================
# all cell markers plot 50 genes
# filter the results
patternRes1 <- patternRes
patternRes1$pvalue[is.na(patternRes1$pvalue)] <- 1
patternRes1 <- patternRes1[order(patternRes1$waldStat,decreasing = T),]
patternRes1$FDR <- p.adjust(patternRes1$pvalue,method = "BH")
patternRes1 <- patternRes1[ patternRes1$FDR < 0.0001 , ]
head(patternRes1,20)

# select the first 50 genes
gene_list <- rownames(patternRes1)[1:min(50,nrow(patternRes1))]
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
png(filename = paste0(cell_type, "_trajectory_branch_50markers.png"),width = 800*2,height = 800*3,res = 300)
x <- t(apply(res,1,function(x){scale(x,T,T)}))
to_order <- apply(x,1,function(x){ mean( (1:length(x))[order(x,decreasing = T)][1] )})
to_order <- apply( (res[,1:20] - res[,22:41]),1,function(x){ mean( (1:length(x))[order(x,decreasing = F)][1] )})
to_order <- apply(x,1,function(x){ which.max(x) })
o <- order(to_order)
o <- hclust(as.dist((1-cor(t(res),use = "pairwise.complete.obs"))/2),method = method )$order
par(mar=c(2,2,3,6))
image( t(x[o,][nrow(res):1,]) ,axes=F,frame=F , col=c("grey95",colorRampPalette(c("grey95","grey70","firebrick","firebrick4"))(90)))
text( par("usr")[2], seq(nrow(res),0,length.out = nrow(res)) / nrow(res) , labels = rownames(x)[o] , pos=4 ,xpd=T , cex=.8)
text( par("usr")[1], par("usr")[4]-diff(par("usr")[4:3])/30, labels = "lineage 2" , adj= c(0,0),xpd=T , cex=1 , pos=3)
text( par("usr")[2], par("usr")[4]-diff(par("usr")[4:3])/30, labels = "lineage 1" , adj= c(1,0),xpd=T , cex=1, pos=3)
text( mean(par("usr")[1:2]), par("usr")[4]-diff(par("usr")[4:3])/30, labels = "origin" , adj= c(1,0),xpd=T , cex=.8, pos=3)
x <- curves@metadata[["lineages"]][["Lineage1"]]
text( seq(mean(par("usr")[1:2]),par("usr")[1],length.out = length(x)), par("usr")[4], labels = x , adj= c(0,0),xpd=T , cex=.8 , pos=3)
x <- curves@metadata[["lineages"]][["Lineage2"]]
text( seq(mean(par("usr")[1:2]),par("usr")[2],length.out = length(x)), par("usr")[4], labels = x , adj= c(0,0),xpd=T , cex=.8 , pos=3)
dev.off()

#=============================================================================================================================
DefaultAssay(temp) <- "RNA"
levels(temp) <- c("0", "2", "6", "8", "20", "12", "13")
# plot ACTA2 (3x6 inches)
VlnPlot(temp, features = c("ACTA2"), cols = colors3$color, pt.size = 0)

# Run differential expression analysis with MAST
min_cluster <- min(table(temp@active.ident))
lung_cluster.markers <- FindAllMarkers(temp, test.use = "MAST", min.diff.pct = 0.1, max.cells.per.ident = min_cluster, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.1)
# calculate the percentage difference of the positive cells in the indicated cluster, in comparison to the others
lung_cluster.markers$Dpct <- (lung_cluster.markers$pct.1 - lung_cluster.markers$pct.2)

c <-subset(temp, idents =8 )
counts <- as.matrix(c[["RNA"]]@data)
m <- data.frame(rowMeans(counts))
gene <- row.names(c)
m <- cbind(gene, m)
names(m)[2] <- paste0("cluster", "0") 

for (i in c("2", "6", "8", "20", "12", "13")) {
  c1 <-subset(lung, idents =i )
  counts1 <- as.matrix(c1[["RNA"]]@data)
  m1 <- data.frame(rowMeans(counts1))
  names(m1)[1] <- paste0("cluster", (i))
  m <- cbind(m, m1)
}
rm(c, counts, counts1, c1, m1)


lung_cluster.markers <-join(lung_cluster.markers, m, by="gene")
write.csv(lung_cluster.markers, file="arw_sm_markers.csv")


lung_cluster.markers <- lung_cluster.markers[lung_cluster.markers$p_val_adj<0.001, ]
# analysis based of the top markers according to the avg_log2fc(10 genes) and then the difference in the percent of positives (5genes)
top2 <- lung_cluster.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top2 <- top2 %>% group_by(cluster) %>% top_n(n = 5, wt = Dpct)
top2 <-top2[, c(7,6)]
top2 <- top2[!duplicated(top2[ , "gene"]),]

levels(temp) <- c("13", "12", "20", "8", "6", "2", "0")
new_order <- c("0", "2", "6", "8", "20", "12", "13")

e <- top2[top2$cluster %in% new_order[1],]
for (i in 2:length(new_order)){
  e1 <- top2[top2$cluster %in% new_order[i],]
  e <- rbind(e, e1)
}

e <- c(e$gene)
e <- e[!duplicated(e)]
#3x9 inches 
DotPlot(temp, assay = "RNA", features = e, scale = TRUE) + RotatedAxis()


levels(temp) <- c("13", "12", "20", "8", "6", "2", "0")
#2.5x5
DotPlot(temp, assay = "RNA", features = c("COL9A1", "MATN2", "FBLN7", "FBN2", "FBN3", "WNT5A", "MKI67", "PCNA"), scale = TRUE) + RotatedAxis()

pdf(file = "dotplot_ecm.pdf", width = 5, height = 2.5)
dotplot_ecm <- DotPlot(temp, assay = "RNA", features = c("ACTA2", "TAGLN", "COL9A1", "MATN2", "FBLN7", "FBN2", "FBN3" , "MKI67", "PCNA"), scale = TRUE) + RotatedAxis()
plot(dotplot_ecm)
dev.off()
rm(dotplot_ecm)

levels(temp) <- c("0",  "2",  "6",  "8",  "20", "12", "13")
#2.5x5
DotPlot(temp, assay = "RNA", features = c("COL9A1", "MATN2", "FBLN7", "FBN2", "FBN3", "WNT5A", "MKI67", "PCNA"), scale = TRUE) + RotatedAxis()

pdf(file = "dotplot_ecm.pdf1", width = 5, height = 2.5)
dotplot_ecm <- DotPlot(temp, assay = "RNA", features = c("ACTA2", "TAGLN", "COL9A1", "MATN2", "FBLN7", "FBN2", "FBN3" , "MKI67", "PCNA"), scale = TRUE) + RotatedAxis()
plot(dotplot_ecm)
dev.off()
rm(dotplot_ecm)
#========================================================
# plot top5 TFs
TFs <- read.csv("TFs.csv")

# create a new file by merging the TF-list and the enriched genes in each cluster, to identify the TFs, only
lung_cluster_tfs <- join(lung_cluster.markers, TFs, by="gene", type="inner")

# analysis based of the top markers according to the avg_log2fc(10 genes) and then the difference in the percent of positives (5genes)
top3 <- lung_cluster_tfs %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top3 <- top3 %>% group_by(cluster) %>% top_n(n = 5, wt = Dpct)
top3 <-top3[, c(7,6)]
top3 <- top3[!duplicated(top3[ , "gene"]),]

levels(temp) <- c("13", "12", "20", "8", "6", "2", "0")
new_order <- c("0", "2", "6", "8", "20", "12", "13")

f <- top3[top3$cluster %in% new_order[1],]
for (i in 2:length(new_order)){
  f1 <- top3[top3$cluster %in% new_order[i],]
  f <- rbind(f, f1)
}

f <- c(f$gene)
f <- f[!duplicated(f)]
#4x8 inches 
DotPlot(temp, assay = "RNA", features = f, scale = TRUE) + RotatedAxis()

#========================================================
# plot top5 secr
secr <- data.frame(read.delim("protein_class_secreted.tsv"))
secr <- secr[1:10]

# create a new file by merging the secreted protein-list and the enriched genes in each cluster, to identify the secreted proteins
lung_cluster_secr <- join(lung_cluster.markers, secr, by="gene", type="inner")

# analysis based of the top markers according to the avg_log2fc(10 genes) and then the difference in the percent of positives (5genes)
top4 <- lung_cluster_secr %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top4 <- top4 %>% group_by(cluster) %>% top_n(n = 5, wt = Dpct)
top4 <-top4[, c(7,6)]
top4 <- top4[!duplicated(top4[ , "gene"]),]

levels(temp) <- c("13", "12", "20", "8", "6", "2", "0")
new_order <- c("0", "2", "6", "8", "20", "12", "13")

h <- top4[top4$cluster %in% new_order[1],]
for (i in 2:length(new_order)){
  h1 <- top4[top4$cluster %in% new_order[i],]
  h <- rbind(h, h1)
}

h <- c(h$gene)
h <- h[!duplicated(h)]
#5x20 inches 
DotPlot(temp, assay = "RNA", features = h, scale = TRUE) + RotatedAxis()

#========================================================
# plot top5 CDs
cds <- data.frame(read.delim("protein_class_CD.tsv"))
cds <- cds[1:10]

# create a new file by merging the CD-list and the enriched genes in each cluster, to identify onlt the CDs
lung_cluster_cds <- join(lung_cluster.markers, cds, by="gene", type="inner")

# analysis based of the top markers according to the avg_log2fc(10 genes) and then the difference in the percent of positives (5genes)
top5 <- lung_cluster_cds %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top5 <- top5 %>% group_by(cluster) %>% top_n(n = 5, wt = Dpct)
top5 <-top5[, c(7,6)]
top5 <- top5[!duplicated(top5[ , "gene"]),]

levels(temp) <- c("13", "12", "20", "8", "6", "2", "0")
new_order <- c("0", "2", "6", "8", "20", "12", "13")

g <- top5[top5$cluster %in% new_order[1],]
for (i in 2:length(new_order)){
  g1 <- top5[top5$cluster %in% new_order[i],]
  g <- rbind(g, g1)
}

g <- c(g$gene)
g <- g[!duplicated(g)]
#5x20 inches 
DotPlot(temp, assay = "RNA", features = g, scale = TRUE) + RotatedAxis()

save.image(file="arw_sm_analysis_data.RData")

save(temp, file="hdca_mes_sm_traj_object.RData")
write.csv(temp@meta.data, file="hdca_mes_traj_metadata.csv")
