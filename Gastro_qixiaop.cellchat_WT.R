#https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/CellChat-vignette.html#load-the-required-libraries
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
GastroSR1.downsampled.cells_WT=subset(GastroSR1.downsampled.cells,subset=orig.ident=="WT")
data.input=GastroSR1.downsampled.cells_WT@assays$RNA@data
meta=cbind(GastroSR1.downsampled.cells_WT@meta.data,labels=GastroSR1.downsampled.cells_WT@active.ident)
cellchat2 <- createCellChat(object = data.input,meta = meta,group.by = "labels")

cellchat2 <- addMeta(cellchat2, meta = meta)
cellchat2 <- setIdent(cellchat2, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat2@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat2@idents)) # number of cells in each cell group

cellchat2DB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(cellchat2DB)


# Show the structure of the database
dplyr::glimpse(cellchat2DB$interaction)

# use a subset of cellchat2DB for cell-cell communication analysis
cellchat2DB.use <- subsetDB(cellchat2DB, search = "Secreted Signaling") # use Secreted Signaling
# use all cellchat2DB for cell-cell communication analysis
# cellchat2DB.use <- cellchat2DB # simply use the default cellchat2DB

# set the used database in the object
cellchat2@DB <- cellchat2DB.use

# subset the expression data of signaling genes for saving computation cost
cellchat2 <- subsetData(cellchat2) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 4) # do parallel
cellchat2 <- identifyOverExpressedGenes(cellchat2)
cellchat2 <- identifyOverExpressedInteractions(cellchat2)

# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
# cellchat2 <- projectData(cellchat2, PPI.human)

cellchat2 <- computeCommunProb(cellchat2)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat2 <- filterCommunication(cellchat2, min.cells = 10)

cellchat2 <- computeCommunProbPathway(cellchat2)
cellchat2 <- aggregateNet(cellchat2)
groupSize <- as.numeric(table(cellchat2@idents))
par(mfrow = c(1,1), xpd=TRUE)
cell_type_cols=c("#00BFFF","#FFA500","#DDA0DD","#FF69B4","#6A5ACD","#9932CC",
                 "#32CD32","#DB7093","#808080","#8FBC8F","#A0522D","#C71585",
                 "#87CEFA","#00CED1","#8B0000","#BC8F8F","#F08080","#DAA520",
                 "#90EE90","#3CB371","#DC143C","#191970")

netVisual_circle(cellchat2@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions",color.use = cell_type_cols,vertex.label.cex = 1)+scale_color_manual(values = cell_type_cols) 
netVisual_circle(cellchat2@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength",color.use = cell_type_cols)+scale_color_manual(values = cell_type_cols) 


mat <- cellchat2@net$weight
par(mfrow = c(2,3), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

pathways.show <- c("IFN-II") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat2, signaling = pathways.show,  vertex.receiver = vertex.receiver)
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat2, signaling = pathways.show, layout = "circle")

# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat2, signaling = pathways.show, layout = "chord",color.use = cell_type_cols)+scale_color_manual(values = cell_type_cols) 

# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat2, signaling = pathways.show, color.heatmap = "Reds")
#> Do heatmap based on a single object

# Chord diagram
group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4)) # grouping cell clusters into fibroblast, DC and TC cells
names(group.cellType) <- levels(cellchat2@idents)
netVisual_chord_cell(cellchat2, signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network"))
#> Plot the aggregated cell-cell communication network at the signaling pathway level

netAnalysis_contribution(cellchat2, signaling = pathways.show)

pairLR.CXCL <- extractEnrichedLR(cellchat2, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[1,] # show one ligand-receptor pair
# Hierarchy plot
vertex.receiver = seq(1,4) # a numeric vector
netVisual_individual(cellchat2, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
#> [[1]]
# Circle plot
netVisual_individual(cellchat2, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

#> [[1]]
# Chord diagram
netVisual_individual(cellchat2, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")

# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat2@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat2@idents)
vertex.receiver = seq(1,4)
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(cellchat2, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat2, signaling = pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
}

# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_bubble(cellchat2, sources.use = 4, targets.use = c(5:11), remove.isolate = FALSE)
#> Comparing communications on a single object

# show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_bubble(cellchat2, sources.use = 4, targets.use = c(5:11), signaling = c("CCL","CXCL"), remove.isolate = FALSE)
#> Comparing communications on a single object

pairLR.use <- extractEnrichedLR(cellchat2, signaling = c("CCL","CXCL","FGF"))
netVisual_bubble(cellchat2, sources.use = c(3,4), targets.use = c(5:8), pairLR.use = pairLR.use, remove.isolate = TRUE)

# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
# show all the interactions sending from Inflam.FIB
netVisual_chord_gene(cellchat2, sources.use = 4, targets.use = c(5:11), lab.cex = 0.5,legend.pos.y = 30)

# show all the interactions received by Inflam.DC
netVisual_chord_gene(cellchat2, sources.use = c(1,2,3,4), targets.use = 8, legend.pos.x = 15)

# show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_chord_gene(cellchat2, sources.use = c(1,2,3,4), targets.use = c(5:11), signaling = c("CCL","CXCL"),legend.pos.x = 8)


plotGeneExpression(cellchat2, signaling = "CXCL")

# Compute the network centrality scores
cellchat2 <- netAnalysis_computeCentrality(cellchat2, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat2, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)


# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat2)
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cellchat2, signaling = c("CXCL", "CCL"))
#> Signaling role analysis on the cell-cell communication network from user's input
gg1 + gg2

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat2, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat2, pattern = "incoming")
ht1 + ht2

# Signaling role analysis on the cell-cell communication networks of interest
ht <- netAnalysis_signalingRole_heatmap(cellchat2, signaling = c("CXCL", "CCL"))

library(NMF)
library(ggalluvial)
selectK(cellchat2, pattern = "outgoing")
nPatterns = 3
cellchat2 <- identifyCommunicationPatterns(cellchat2, pattern = "outgoing", k = nPatterns)

# river plot
netAnalysis_river(cellchat2, pattern = "outgoing")
#> Please make sure you have load `library(ggalluvial)` when running this function


# dot plot
netAnalysis_dot(cellchat2, pattern = "outgoing")

selectK(cellchat2, pattern = "incoming")
nPatterns = 4
cellchat2 <- identifyCommunicationPatterns(cellchat2, pattern = "incoming", k = nPatterns)

# river plot
netAnalysis_river(cellchat2, pattern = "incoming")
#> Please make sure you have load `library(ggalluvial)` when running this function

netAnalysis_dot(cellchat2, pattern = "incoming")

cellchat2 <- computeNetSimilarity(cellchat2, type = "functional")

cellchat2 <- netEmbedding(cellchat2, type = "functional")
#> Manifold learning of the signaling networks for a single dataset
cellchat2 <- netClustering(cellchat2, type = "functional")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(cellchat2, type = "functional", label.size = 3.5)

#> Manifold learning of the signaling networks for a single dataset
cellchat2 <- netClustering(cellchat2, type = "functional")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(cellchat2, type = "functional", label.size = 3.5)

cellchat2 <- computeNetSimilarity(cellchat2, type = "structural")
cellchat2 <- netEmbedding(cellchat2, type = "structural")
#> Manifold learning of the signaling networks for a single dataset
cellchat2 <- netClustering(cellchat2, type = "structural")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(cellchat2, type = "structural", label.size = 3.5)

netVisual_embeddingZoomIn(cellchat2, type = "structural", nCol = 2)

#saveRDS(cellchat2, file = "cellchat2_humanSkin_LS.rds")

