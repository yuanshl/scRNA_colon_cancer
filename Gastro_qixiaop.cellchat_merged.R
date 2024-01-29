#https://www.jianshu.com/p/49a0a0b50987
#https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/CellChat-vignette.html#load-the-required-libraries
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
GastroSR1.downsampled.cells_WT=subset(GastroSR1.downsampled.cells,subset=orig.ident=="WT")
GastroSR1.downsampled.cells_WT$celltype=GastroSR1.downsampled.cells_WT@active.ident
GastroSR1.downsampled.cells_WT=subset(GastroSR1.downsampled.cells_WT,subset=celltype!=c("Krt+ epithelia cells"))
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
cellchat2@DB <- cellchat2DB.use

# subset the expression data of signaling genes for saving computation cost
cellchat2 <- subsetData(cellchat2) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 4) # do parallel
cellchat2 <- identifyOverExpressedGenes(cellchat2)
cellchat2 <- identifyOverExpressedInteractions(cellchat2)
cellchat2 <- computeCommunProb(cellchat2)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat2 <- filterCommunication(cellchat2, min.cells = 10)
cellchat2 <- computeCommunProbPathway(cellchat2)
cellchat2 <- aggregateNet(cellchat2)
groupSize <- as.numeric(table(cellchat2@idents))





#############################
######
#############################
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

cellchat1 <- computeCommunProb(cellchat1)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat1 <- filterCommunication(cellchat1, min.cells = 10)
cellchat1 <- computeCommunProbPathway(cellchat1)
cellchat1 <- aggregateNet(cellchat1)
groupSize1 <- as.numeric(table(cellchat1@idents))


object.list=list(WT=cellchat2,Ncf4_KO=cellchat1)
cellchat=mergeCellChat(object.list,add.names = names(object.list))

cell_type_cols1=c("#00BFFF","#FFA500","#DDA0DD","#FF69B4","#6A5ACD","#9932CC",
                 "#32CD32","#DB7093","#808080","#A0522D","#C71585",
                 "#87CEFA","#00CED1","#8B0000","#BC8F8F","#F08080","#DAA520",
                 "#90EE90","#3CB371","#DC143C","#191970")

gg1=compareInteractions(cellchat,show.legend = F, group = c(1,2))
gg2=compareInteractions(cellchat,show.legend = F,measure = "weight", group = c(1,2))
gg1+gg2

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T,color.use = cell_type_cols1)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight",color.use = cell_type_cols1)

gg1 <- netVisual_heatmap(cellchat,color.use = cell_type_cols1[1:19])
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight",color.use = cell_type_cols1[1:19])
#> Do heatmap based on a merged object
gg1 + gg2

pathways.show <- c("IFN-II")
gg1 <- netVisual_heatmap(cellchat,signaling = pathways.show,color.use = cell_type_cols1[1:19])
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, signaling = pathways.show,measure = "weight",color.use = cell_type_cols1[1:19])
#> Do heatmap based on a merged object
gg1 + gg2

weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}




