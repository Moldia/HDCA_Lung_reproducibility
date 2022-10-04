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

# increase the maximum of allowed RAM
options(future.globals.maxSize = 30000 * 1024^2)

load("all_cells_subset_dataset_without_doublets.RData")
lung <-lung1
rm(lung1)
DefaultAssay(lung) <- "RNA"
lung@assays[["SCT"]]<- NULL
lung@assays[["integrated"]]<- NULL
lung@reductions[["pca"]] <-NULL
lung@reductions[["umap"]] <-NULL
lung@reductions[["umap2d"]] <-NULL
lung@reductions[["tsne"]] <-NULL

#====================================================
# plot the 
levels(lung@active.ident)
lung <- subset(lung, ident= c("epi_cl2", "epi_cl10", "epi_cl3", "epi_cl9", "epi_cl8", "epi_cl13", "epi_cl5", 
                              "mes_cl0", "mes_cl2","mes_cl6", "mes_cl8", "mes_cl20", "mes_cl12"))

levels(lung) <- c("epi_cl2", "epi_cl10", "epi_cl3", "epi_cl9", "epi_cl8", "epi_cl13", "epi_cl5", 
                  "mes_cl0", "mes_cl2","mes_cl6", "mes_cl8", "mes_cl20", "mes_cl12")
# to be able to do it with the correct values, we have to normalize and scale the RNA counts
lung <- NormalizeData(lung, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(lung)
lung <- ScaleData(lung, features = all.genes)

# script has been based on https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/CellChat-vignette.html
# extract the data from Seurat
data.input <- GetAssayData(lung, assay = "RNA", slot = "data") # normalized data matrix

#===============================================
clusters <- levels(lung@active.ident)
c <-subset(lung, idents =clusters[1] )
counts <- as.matrix(c[["RNA"]]@data)
m <- data.frame(rowMeans(counts))
gene <- row.names(c)
m <- cbind(gene, m)
names(m)[2] <- clusters[1]


for (i in 2:13) {
  name <- clusters[i]
  c1 <-subset(lung, idents =name )
  counts1 <- as.matrix(c1[["RNA"]]@data)
  m1 <- data.frame(rowMeans(counts1))
  names(m1)[1] <- name
  m <- cbind(m, m1)
}
rm(c, counts, counts1, c1, m1, name)

max_value <- data.frame(apply(m[, 2:14], 1, max))
names(max_value) <- "mean"
max_value$gene <- row.names(max_value)
max_value <- max_value[max_value$mean >= 0.3,]
#=============================================
# keep only the filtered genes
data.input <- data.input[row.names(data.input) %in% row.names(max_value),]

labels <- Idents(lung)
meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels

# create cellchat object
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")

levels(cellchat) <-c("epi_cl2", "epi_cl10", "epi_cl3", "epi_cl9", "epi_cl8", "epi_cl13", "epi_cl5", 
                     "mes_cl0", "mes_cl2","mes_cl6", "mes_cl8", "mes_cl20", "mes_cl12")

#Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
# save the ligand-receptor_interations_file
save(CellChatDB, file= "CellChatDB_human_20210520.RData")

# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)
#> Rows: 1,939
#> Columns: 11

# use a subset of CellChatDB for cell-cell communication analysis
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# set the used database in the object
cellchat@DB <- CellChatDB.use

# Preprocessing the expression data for cell-cell communication analysis

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
#future::plan("multiprocess", workers = 4) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
# project gene expression data onto PPI network (optional)
cellchat <- projectData(cellchat, PPI.human)

# Part II: Inference of cell-cell communication network
# The function computeAveExpr can help to check the average expression of signaling genes of interest, e.g, 
computeAveExpr(cellchat, features = c("FGF20","BMP2"), type =  "truncatedMean", trim = 0.1)

# Compute the communication probability and infer cellular communication network
cellchat <- computeCommunProb(cellchat)
save.image("after_communProbe1.RData")

load("after_communProbe1.RData")
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
# obtain the size of analyzed clusters
table(lung@active.ident)
min_size <- as.integer(min(table(lung@active.ident))/2)
# select only the communications that are present in a number of cells that is defined as the 50% of size of the smallest cluster 
cellchat <- filterCommunication(cellchat, min.cells = min_size)

# Extract the inferred cellular communication network as a data frame
# CellChat provides a function subsetCommunication to easily access the inferred cell-cell communications of interest. For example,
df.net <- subsetCommunication(cellchat) # returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors. 

# Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)

# Calculate the aggregated cell-cell communication network
# CellChat can calculate the aggregated cell-cell communication network by counting the number of links or summarizing the communication probability. 
# USER can also calculate the aggregated network among a subset of cell groups by setting sources.use and targets.use.
cellchat <- aggregateNet(cellchat)

# We can also visualize the aggregated cell-cell communication network. For example, showing the number of interactions or the total interaction strength 
# (weights) between any two cell groups using circle plot.
groupSize <- as.numeric(table(cellchat@idents))

# plot and export the number of interactions
par(mfrow=c(1,1))
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
tiff("number_of_interactions_circ.tiff",width = 1000, height = 1000, compression = "lzw")
number_of_interactions_circ <-netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
plot(number_of_interactions_circ)
rm(number_of_interactions_circ)
dev.off()
pdf("number_of_interactions_circ.pdf",width = 6, height = 6, paper = "special")
number_of_interactions_circ <-netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
plot(number_of_interactions_circ)
rm(number_of_interactions_circ)
dev.off()

# plot and export the Interaction weights
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
tiff("weights_of_interactions_circ.tiff",width = 1000, height = 1000, compression = "lzw")
weights_of_interactions_circ <-netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
plot(weights_of_interactions_circ)
rm(weights_of_interactions_circ)
dev.off()
pdf("weights_of_interactions_circ.pdf",width = 6, height = 6, paper = "special")
weights_of_interactions_circ <-netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
plot(weights_of_interactions_circ)
rm(weights_of_interactions_circ)
dev.off()

# Due to the complicated cell-cell communication network, we can examine the signaling sent from each cell group. Here we also control the parameter edge.weight.max so 
# that we can compare edge weights between different networks.

mat <- cellchat@net$weight

for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  name <- paste0("single_interaction_plots/", rownames(mat)[i], "_interactions_circ.pdf")
  pdf(name,width = 6, height = 6, paper = "special")
  print(netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i]))
  dev.off()
}

# Part III: Visualization of cell-cell communication network
#  Visualize each signaling pathway using Hierarchy plot, Circle plot or Chord diagram

# All the signaling pathways showing significant communications can be accessed by:
cellchat@netP$pathways

# one signaling pathway can be selected as an example for downstream analysis
pathways.show <- c("WNT") 
# Hierarchy plot
# Here the `vertex.receive` is defined so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
# obtain the order of the clusters, to select the ones of interest
rownames(cellchat@netP[["prob"]])

vertex.receiver = seq(1,7) # a numeric vector.
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)

# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

tiff(paste0(pathways.show, "_signaling_circ.tiff") ,width = 1000, height = 1000, compression = "lzw")
a<- netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
print(a)
rm(a)
dev.off()

pdf(paste0(pathways.show, "_signaling_circ.pdf"),width = 6, height = 6, paper = "special")
a<- netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
print(a)
rm(a)
dev.off()

# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
#> Note: The first link end is drawn out of sector 'Inflam. FIB'.

tiff(paste0(pathways.show, "_signaling_chord.tiff") ,width = 1000, height = 1000, compression = "lzw")
a<- netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
print(a)
rm(a)
dev.off()

png(paste0(pathways.show, "_signaling_chord.png") ,width = 1000, height = 1000)
a<- netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
print(a)
rm(a)
dev.off()


pdf(paste0(pathways.show, "_signaling_chord.pdf"),width = 6, height = 6, paper = "special")
a<- netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
print(a)
rm(a)
dev.off()

# Chord diagram
group.cellType <- c(rep("epi", 7), rep("mes", 6)) # grouping cell clusters into mesenchymal and epithelial
names(group.cellType) <- levels(cellchat@idents)
netVisual_chord_cell(cellchat, signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network"))

pdf(paste0(pathways.show, "_signaling_chord_grouped.pdf"),width = 6, height = 6, paper = "special")
a<- netVisual_chord_cell(cellchat, signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network"))
print(a)
rm(a)
dev.off()


# Heatmap
pdf(paste0(pathways.show, "_signaling_heatmap.pdf"),width = 11.69, height = 8.27, paper = "special")
a <- netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
print(a)
rm(a)
dev.off()

# Compute the contribution of each ligand-receptor pair to the overall signaling pathway and visualize cell-cell communication mediated by a single ligand-receptor pair
netAnalysis_contribution(cellchat, signaling = pathways.show)
pdf(paste0(pathways.show, "_contribution_LR_pairs_barplot.pdf"),width = 6, height = 6, paper = "special")
a<- netAnalysis_contribution(cellchat, signaling = pathways.show)
print(a)
rm(a)
dev.off()

# CellChat can also visualize the cell-cell communication mediated by a single ligand-receptor pair. It provides a function:
# extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.selected <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.selected[1,] # show one ligand-receptor pair

# obtain the order of the clusters, to select the ones of interest
rownames(cellchat@netP[["prob"]])
vertex.receiver = seq(1,7) # a numeric vector
# Hierarchy plot
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
pdf(paste0(pathways.show, "_significant_LR_pairs_heirarchy_plot.pdf"),width = 8, height = 8, paper = "special")
a <- netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
print(a)
rm(a)
dev.off()

# Circle plot of the selected L-R pair
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
pdf(paste0(pathways.show, "_significant_LR_pairs_circ_plot.pdf"),width = 8, height = 8, paper = "special")
a <- netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
print(a)
rm(a)
dev.off()

# Automatically save the plots of the all inferred network for quick exploration
# In practical use, USERS can use 'for . loop' to automatically save the all inferred network for quick exploration using netVisual. netVisual supports an output in the formats of svg, png and pdf.

# Access all the signaling pathways showing significant communications

pathways.show.all <- cellchat@netP$pathways
pathways.show.all
levels(cellchat@idents)
vertex.receiver = seq(1,7)
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.svg"), plot=gg, width = 6, height = 4, units = 'in', dpi = 300)
}


#===================================================================================================
save.image("after_communProbe1.RData")

# Visualize cell-cell communication mediated by multiple ligand-receptors or signaling pathways
# Bubble plot
# CellChat can also show all the significant interactions (L-R pairs) from some cell groups to other cell groups using netVisual_bubble.

# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
# obtain the order of the clusters, to select the ones of interest
rownames(cellchat@netP[["prob"]])
netVisual_bubble(cellchat, sources.use = c(8:13), targets.use = c(1:7), remove.isolate = FALSE, thresh = 0.001)

# show all the significant interactions (L-R pairs) associated with certain signaling pathways
pathways.show.all <- cellchat@netP$pathways
netVisual_bubble(cellchat, sources.use = c(8:13), targets.use = c(1:7), signaling = pathways.show.all, remove.isolate = FALSE,thresh = 0.001)

pdf("pathways_LR_probability_dotplot_mes_target.pdf",width = 7, height = 14, paper = "special")
a<- netVisual_bubble(cellchat, sources.use = c(1:7), targets.use = c(8:13), signaling = pathways.show.all, remove.isolate = FALSE,thresh = 0.001,font.size = 8)
print(a)
rm(a)
dev.off()

pdf("pathways_LR_probability_dotplot_epi_target.pdf",width = 7, height = 14, paper = "special")
a<- netVisual_bubble(cellchat, sources.use = c(8:13), targets.use = c(1:7), signaling = pathways.show.all, remove.isolate = FALSE,thresh = 0.001,font.size = 8)
print(a)
rm(a)
dev.off()

pdf("pathways_LR_probability_dotplot_mes.pdf",width = 7, height = 14, paper = "special")
a<- netVisual_bubble(cellchat, sources.use = c(8:13), targets.use = c(8:13), signaling = pathways.show.all, remove.isolate = FALSE,thresh = 0.001,font.size = 8)
print(a)
rm(a)
dev.off()

pdf("pathways_LR_probability_dotplot_epi.pdf",width = 7, height = 14, paper = "special")
a<- netVisual_bubble(cellchat, sources.use = c(1:7), targets.use = c(1:7), signaling = pathways.show.all, remove.isolate = FALSE,thresh = 0.001,font.size = 8)
print(a)
rm(a)
dev.off()

# Chord diagram
# Similar to Bubble plot, CellChat provides a function netVisual_chord_gene for drawing Chord diagram to
# show all the interactions (L-R pairs or signaling pathways) from some cell groups to other cell groups.  

# show all the significant signaling pathways from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
pdf("significant_pathways_from_mes_neurons_chordplot.pdf",width = 18, height = 18, paper = "special")
a<- netVisual_chord_gene(cellchat, sources.use = c(8:13), targets.use = c(1:7), slot.name = "netP", legend.pos.x = 10)
print(a)
rm(a)
dev.off()

pdf("significant_pathways_from_epi_chordplot.pdf",width = 18, height = 18, paper = "special")
a<- netVisual_chord_gene(cellchat, sources.use = c(1:7), targets.use = c(8:13), slot.name = "netP", legend.pos.x = 10)
print(a)
rm(a)
dev.off()

#=====================================================================================
# Plot the signaling gene expression distribution using violin/dot plot
# CellChat can plot the gene expression distribution of signaling genes related to L-R pairs or signaling pathway using a Seurat wrapper function plotGeneExpression.

for (i in pathways.show.all) {
pdf(paste0("violin_plots/",i,"_pathway_enriched_genes_violinplots.pdf"),width = 12, height = 12, paper = "special")
a<- plotGeneExpression(cellchat, signaling = i, enriched.only = TRUE)
print(a)
rm(a)
dev.off()
}

#=========================================================================================
#=========================================================================================
# Part IV: Systems analysis of cell-cell communication network

# Identify signaling roles (e.g., dominant senders, receivers) of cell groups as well as the major contributing signaling

# Compute and visualize the network centrality scores
# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

for (i in pathways.show.all) {
pdf(paste0("pathway_heatmaps/", i,"_pathway_signaling_role_heatmap.pdf"),width = 6, height = 3, paper = "special")
a<- netAnalysis_signalingRole_network(cellchat, signaling = i, width = 8, height = 2.5, font.size = 10, color.heatmap = "Reds")
print(a)
rm(a)
dev.off()
}

# Visualize the dominant senders (sources) and receivers (targets) in a 2D space
# We also provide another intutive way to visualize the dominant senders (sources) and receivers (targets) in a 2D space using scatter plot.

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
netAnalysis_signalingRole_scatter(cellchat)
pdf("all_pathways_signaling_role_scatterplot.pdf",width = 5, height = 5, paper = "special")
a<- netAnalysis_signalingRole_scatter(cellchat)
print(a)
rm(a)
dev.off()

#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
netAnalysis_signalingRole_scatter(cellchat, signaling = pathways.show)
pdf(paste0(pathways.show,"_pathway_signaling_role_scatterplot.pdf"),width = 5, height = 5, paper = "special")
a<- netAnalysis_signalingRole_scatter(cellchat, signaling = pathways.show)
print(a)
rm(a)
dev.off()

#-----------------------------------------------------------------------------------
# Identify signals contributing most to outgoing or incoming signaling of certain cell groups
# CellChat can also answer the question on which signals contributing most to outgoing or incoming signaling of certain cell groups.

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
pdf("outgoing_signaling_patterns_heatmap.pdf",width = 10, height = 10, paper = "special")
a<- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing",width = 5, height = 15, color.heatmap = "Reds")
print(a)
rm(a)
dev.off()

netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
pdf("incoming_signaling_patterns_heatmap.pdf",width = 10, height = 10, paper = "special")
a<- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming",width = 5, height = 15, color.heatmap = "Reds")
print(a)
rm(a)
dev.off()

# Signaling role analysis on the cell-cell communication networks of interest
netAnalysis_signalingRole_heatmap(cellchat, signaling = pathways.show.all)
pdf("selected_pathways_outgoing_signaling_patterns_heatmap.pdf",width = 10, height = 10, paper = "special")
a<- netAnalysis_signalingRole_heatmap(cellchat, signaling = pathways.show.all, pattern = "outgoing",width = 5, height = 15, color.heatmap = "Reds")
print(a)
rm(a)
dev.off()

pdf("selected_pathways_incoming_signaling_patterns_heatmap.pdf",width = 10, height = 10, paper = "special")
a<- netAnalysis_signalingRole_heatmap(cellchat, signaling = pathways.show.all, pattern = "incoming",width = 5, height = 15, color.heatmap = "Reds")
print(a)
rm(a)
dev.off()

#======================================================================================
#======================================================================================
# Identify global communication patterns to explore how multiple cell types and signaling pathways coordinate together

# Identify and visualize outgoing communication pattern of secreting cells
# Outgoing patterns reveal how the sender cells (i.e. cells as signal source) coordinate with each other as well as how
# they coordinate with certain signaling pathways to drive communication.
selectK(cellchat, pattern = "outgoing")

nPatterns = 4
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)

netAnalysis_river(cellchat, pattern = "outgoing")


pdf("communication_patterns_outgoing_heatmap.pdf",width = 12, height = 12, paper = "special")
a<- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns, width = 8, height = 16, font.size = 8)
print(a)
rm(a)
dev.off()

# river plot
pdf("communication_patterns_outgoing_riverplot.pdf",width = 12, height = 12, paper = "special")
a<- netAnalysis_river(cellchat, pattern = "outgoing")
print(a)
rm(a)
dev.off()

netAnalysis_dot(cellchat, pattern = "outgoing")
# dot plot
pdf("communication_patterns_outgoing_dotplot.pdf",width = 8, height = 2.5, paper = "special")
a<- netAnalysis_dot(cellchat, pattern = "outgoing",main.title = "Outgoing communication patterns")
print(a)
rm(a)
dev.off()

# Identify and visualize incoming communication pattern of target cells
# Incoming patterns show how the target cells (i.e. cells as signal receivers) coordinate with each other as well as how they 
# coordinate with certain signaling pathways to respond to incoming signals.

selectK(cellchat, pattern = "incoming")
# the number patterns is set to the point that the line starts to decrease
nPatterns = 4
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns,width = 8, height = 14, font.size = 6)

pdf("communication_patterns_incoming_heatmap.pdf",width = 12, height = 12, paper = "special")
a<- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns,width = 8, height = 16, font.size = 8)
print(a)
rm(a)
dev.off()

# river plot
pdf("communication_patterns_incoming_riverplot.pdf",width = 12, height = 12, paper = "special")
a<- netAnalysis_river(cellchat, pattern = "incoming")
print(a)
rm(a)
dev.off()

# dot plot
pdf("communication_patterns_incoming_dotplot.pdf",width = 12, height = 12, paper = "special")
a<- netAnalysis_dot(cellchat, pattern = "incoming")
print(a)
rm(a)
dev.off()
#================================================================================================
#================================================================================================
# Manifold and classification learning analysis of signaling networks
# Further, CellChat is able to quantify the similarity between all significant signaling
# pathways and then group them based on their cellular communication network similarity. 
# Grouping can be done either based on the functional or structural similarity.
save.image("after_communProbe1.RData")
# Identify signaling groups based on their functional similarity
cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
#> Manifold learning of the signaling networks for a single dataset
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space

pdf("signaling_pathways_functional_similarities_plot.pdf",width = 8, height = 8, paper = "special")
a<- netVisual_embedding(cellchat, type = "functional", label.size = 3.5)
print(a)
rm(a)
dev.off()

#Identify signaling groups based on structure similarity
cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
#> Manifold learning of the signaling networks for a single dataset
cellchat <- netClustering(cellchat, type = "structural")
#> Classification learning of the signaling networks for a single dataset

# Visualization in 2D-space
pdf("signaling_pathways_structural_similarities_plot.pdf",width = 8, height = 8, paper = "special")
a<- netVisual_embedding(cellchat, type = "structural", label.size = 3.5)
print(a)
rm(a)
dev.off()

save.image("final_analysis_data.RData")
load("final_analysis_data.RData")

saveRDS(cellchat, file = "distal_cellchat_object.rds")



#===============================================================================================================
DefaultAssay(lung) <- "RNA"
levels(lung) <- c("epi_cl2", "epi_cl10", "epi_cl3", "epi_cl9", "epi_cl8", "epi_cl13", "epi_cl5", 
                  "mes_cl0", "mes_cl2","mes_cl6", "mes_cl8", "mes_cl20", "mes_cl12")
DotPlot(lung, assay = "RNA", features = c("NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4", #receptors
                                          "HES1", "HES2", "HES3", "HES5", "HES6", "HES7", "HEY1", "HEY2", "HEYL", "NRARP", #targets
                                          "DLL1", "DLL4", "JAG1", "JAG2", #ligands
                                          "POFUT1", "POFUT2", "LFNG", "MFNG", "RFNG", "MAML1", "MAML2", "MAML3", "MAMLD1", "RBPJ", #transducers
                                          "DLL3", "DLK1", "DLK2", "NUMB", "NUMBL" # inhibitors
                                          
                                          ), scale = T) + RotatedAxis()

lung1 <- subset(lung1, ident= c("epi_cl2", "epi_cl10", "epi_cl3", "epi_cl9", "epi_cl8", "epi_cl13", "epi_cl5", 
                              "mes_cl0", "mes_cl2","mes_cl6", "mes_cl8", "mes_cl20", "mes_cl12", "mes_cl13", "mes_cl5", "mes_cl16"))

levels(lung1) <- c("epi_cl2", "epi_cl10", "epi_cl3", "epi_cl9", "epi_cl8", "epi_cl13", "epi_cl5", 
                  "mes_cl0", "mes_cl2","mes_cl6", "mes_cl8", "mes_cl20", "mes_cl12", "mes_cl13", "mes_cl5", "mes_cl16")
DefaultAssay(lung1) <- "RNA"
DotPlot(lung1, assay = "RNA", features = c("NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4", #receptors
                                          "HES1", "HES2", "HES3", "HES5", "HES6", "HES7", "HEY1", "HEY2", "HEYL", "NRARP", #targets
                                          "DLL1", "DLL4", "JAG1", "JAG2", #ligands
                                          "POFUT1", "POFUT2", "LFNG", "MFNG", "RFNG", "MAML1", "MAML2", "MAML3", "MAMLD1", "RBPJ", #transducers
                                          "DLL3", "DLK1", "DLK2", "NUMB", "NUMBL" # inhibitors
                                          
), scale = F) + RotatedAxis()


lung1 <- subset(lung1, ident= c("mes_cl0", "mes_cl2","mes_cl6", "mes_cl8", "mes_cl20", "mes_cl12", "mes_cl13"))

levels(lung1) <- c("mes_cl13", "mes_cl12", "mes_cl20", "mes_cl8",  "mes_cl6",  "mes_cl2",  "mes_cl0")
DefaultAssay(lung1) <- "RNA"
DotPlot(lung1, assay = "RNA", features = c("NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4", #receptors
                                           "HES1", "HES2", "HES3", "HES5", "HES6", "HES7", "HEY1", "HEY2", "HEYL", "NRARP", #targets
                                           "DLL1", "DLL4", "JAG1", "JAG2", #ligands
                                           "POFUT1", "POFUT2", "LFNG", "MFNG", "RFNG", "MAML1", "MAML2", "MAML3", "MAMLD1", "RBPJ", #transducers
                                           "DLL3", "DLK1", "DLK2", "NUMB", "NUMBL" # inhibitors
                                           
), scale = T) + RotatedAxis()

lung1 <- subset(lung1, ident= c("epi_cl1", "epi_cl10", "epi_cl2", "epi_cl3", "epi_cl9", "epi_cl5", "epi_cl13", "epi_cl8", "mes_cl0", "mes_cl2","mes_cl6", "mes_cl8", "mes_cl20", "mes_cl12", "mes_cl13"))
levels(lung1) <- c("epi_cl1", "epi_cl10", "epi_cl2", "epi_cl3", "epi_cl9", "epi_cl5", "epi_cl13", "epi_cl8", "mes_cl0", "mes_cl2","mes_cl6", "mes_cl8", "mes_cl20", "mes_cl12", "mes_cl13")

DefaultAssay(lung1) <- "RNA"
DotPlot(lung1, assay = "RNA", features = c("FGF2", "FGF7","FGF9", "FGF10", "FGF13", "FGF14", "FGF18", "FGF20",
                                           "FGFR1", "FGFR2", "FGFR3", "FGFR4",
                                           "ETV1", "ETV3", "ETV5", "ETV6", "ETS1", "ETS2", "SPRY1", "SPRY2", "SPRY4"
                                           ), scale = T) + RotatedAxis()

DotPlot(lung1, assay = "RNA", features = c("FGF20", "FGF18", "FGF9",
                                           "FGFR1", "FGFR2", "FGFR3", "FGFR4",
                                           "SPRY2","ETV5","ETV3","ETS1"
), scale = T) + RotatedAxis()
