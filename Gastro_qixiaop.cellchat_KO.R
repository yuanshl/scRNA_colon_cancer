#https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/CellChat-vignette.html#load-the-required-libraries
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
GastroSR1.downsampled.cells_KO=subset(GastroSR1.downsampled.cells,subset=orig.ident=="Ncf4_KO")
data.input=GastroSR1.downsampled.cells_KO@assays$RNA@data
meta=cbind(GastroSR1.downsampled.cells_KO@meta.data,labels=GastroSR1.downsampled.cells_KO@active.ident)
cellchat1 <- createCellChat(object = data.input,meta = meta,group.by = "labels")

cellchat1 <- addMeta(cellchat1, meta = meta)
cellchat1 <- setIdent(cellchat1, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat1@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat1@idents)) # number of cells in each cell group

cellchat1DB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(cellchat1DB)


# Show the structure of the database
dplyr::glimpse(cellchat1DB$interaction)

# use a subset of cellchat1DB for cell-cell communication analysis
cellchat1DB.use <- subsetDB(cellchat1DB, search = "Secreted Signaling") # use Secreted Signaling
# use all cellchat1DB for cell-cell communication analysis
# cellchat1DB.use <- cellchat1DB # simply use the default cellchat1DB

# set the used database in the object
cellchat1@DB <- cellchat1DB.use
cellchat1DB.use$interaction

# subset the expression data of signaling genes for saving computation cost
cellchat1 <- subsetData(cellchat1) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 4) # do parallel
cellchat1 <- identifyOverExpressedGenes(cellchat1)
cellchat1 <- identifyOverExpressedInteractions(cellchat1)

# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
# cellchat1 <- projectData(cellchat1, PPI.human)

cellchat1 <- computeCommunProb(cellchat1)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat1 <- filterCommunication(cellchat1, min.cells = 10)

cellchat1 <- computeCommunProbPathway(cellchat1)
cellchat1 <- aggregateNet(cellchat1)
groupSize1 <- as.numeric(table(cellchat1@idents))
par(mfrow = c(1,1), xpd=TRUE)


cell_type_cols1=c("#00BFFF","#FFA500","#DDA0DD","#FF69B4","#6A5ACD","#9932CC",
                  "#32CD32","#DB7093","#808080","#A0522D","#C71585",
                  "#87CEFA","#00CED1","#8B0000","#BC8F8F","#F08080","#DAA520",
                  "#90EE90","#3CB371","#DC143C","#191970")
netVisual_circle(cellchat1@net$count, vertex.weight = groupSize1, weight.scale = T, label.edge= F, title.name = "Number of interactions",color.use = cell_type_cols1)
netVisual_circle(cellchat1@net$weight, vertex.weight = groupSize1, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength",color.use = cell_type_cols1)


mat <- cellchat1@net$weight
par(mfrow = c(2,3), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

pathways.show <- c("TGFb") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat1, signaling = pathways.show,  vertex.receiver = vertex.receiver)
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat1, signaling = pathways.show, layout = "circle")

# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat1, signaling = pathways.show, layout = "chord",color.use = cell_type_cols1)+scale_color_manual(values = cell_type_cols1) 

# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat1, signaling = pathways.show, color.heatmap = "Reds")
#> Do heatmap based on a single object

# Chord diagram
group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4)) # grouping cell clusters into fibroblast, DC and TC cells
names(group.cellType) <- levels(cellchat1@idents)
netVisual_chord_cell(cellchat1, signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network"))
#> Plot the aggregated cell-cell communication network at the signaling pathway level

netAnalysis_contribution(cellchat1, signaling = pathways.show)

pairLR.CXCL <- extractEnrichedLR(cellchat1, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[1,] # show one ligand-receptor pair
# Hierarchy plot
vertex.receiver = seq(1,4) # a numeric vector
netVisual_individual(cellchat1, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
#> [[1]]
# Circle plot
netVisual_individual(cellchat1, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

#> [[1]]
# Chord diagram
netVisual_individual(cellchat1, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")

# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat1@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat1@idents)
vertex.receiver = seq(1,4)
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(cellchat1, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat1, signaling = pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
}

# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_bubble(cellchat1, sources.use = 4, targets.use = c(5:11), remove.isolate = FALSE)
#> Comparing communications on a single object

# show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_bubble(cellchat1, sources.use = 4, targets.use = c(5:11), signaling = c("CCL","CXCL"), remove.isolate = FALSE)
#> Comparing communications on a single object

pairLR.use <- extractEnrichedLR(cellchat1, signaling = c("CCL","CXCL","FGF"))
netVisual_bubble(cellchat1, sources.use = c(3,4), targets.use = c(5:8), pairLR.use = pairLR.use, remove.isolate = TRUE)

# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
# show all the interactions sending from Inflam.FIB
netVisual_chord_gene(cellchat1, sources.use = 4, targets.use = c(5:11), lab.cex = 0.5,legend.pos.y = 30)

# show all the interactions received by Inflam.DC
netVisual_chord_gene(cellchat1, sources.use = c(1,2,3,4), targets.use = 8, legend.pos.x = 15)

# show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_chord_gene(cellchat1, sources.use = c(1,2,3,4), targets.use = c(5:11), signaling = c("CCL","CXCL"),legend.pos.x = 8)


plotGeneExpression(cellchat1, signaling = "CXCL")

# Compute the network centrality scores
cellchat1 <- netAnalysis_computeCentrality(cellchat1, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat1, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)


# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat1)
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cellchat1, signaling = c("CXCL", "CCL"))
#> Signaling role analysis on the cell-cell communication network from user's input
gg1 + gg2

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat1, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat1, pattern = "incoming")
ht1 + ht2

# Signaling role analysis on the cell-cell communication networks of interest
ht <- netAnalysis_signalingRole_heatmap(cellchat1, signaling = c("CXCL", "CCL"))

library(NMF)
library(ggalluvial)
selectK(cellchat1, pattern = "outgoing")
nPatterns = 3
cellchat1 <- identifyCommunicationPatterns(cellchat1, pattern = "outgoing", k = nPatterns)

# river plot
netAnalysis_river(cellchat1, pattern = "outgoing")
#> Please make sure you have load `library(ggalluvial)` when running this function


# dot plot
netAnalysis_dot(cellchat1, pattern = "outgoing")

selectK(cellchat1, pattern = "incoming")
nPatterns = 4
cellchat1 <- identifyCommunicationPatterns(cellchat1, pattern = "incoming", k = nPatterns)

# river plot
netAnalysis_river(cellchat1, pattern = "incoming")
#> Please make sure you have load `library(ggalluvial)` when running this function

netAnalysis_dot(cellchat1, pattern = "incoming")

cellchat1 <- computeNetSimilarity(cellchat1, type = "functional")

cellchat1 <- netEmbedding(cellchat1, type = "functional")
#> Manifold learning of the signaling networks for a single dataset
cellchat1 <- netClustering(cellchat1, type = "functional")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(cellchat1, type = "functional", label.size = 3.5)

#> Manifold learning of the signaling networks for a single dataset
cellchat1 <- netClustering(cellchat1, type = "functional")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(cellchat1, type = "functional", label.size = 3.5)

cellchat1 <- computeNetSimilarity(cellchat1, type = "structural")
cellchat1 <- netEmbedding(cellchat1, type = "structural")
#> Manifold learning of the signaling networks for a single dataset
cellchat1 <- netClustering(cellchat1, type = "structural")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(cellchat1, type = "structural", label.size = 3.5)

netVisual_embeddingZoomIn(cellchat1, type = "structural", nCol = 2)

#saveRDS(cellchat1, file = "cellchat1_humanSkin_LS.rds")

