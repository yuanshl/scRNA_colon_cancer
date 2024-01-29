library(Seurat)
library(dplyr)
library(patchwork)
library(mindr)
library(Matrix)
library(SeuratDisk)
library(patchwork)
library(dplyr)
library(ggplot2)
library(future)

Ncf4_KO="/media/liuhb/cell/liuhb/R_data/qixiaopeng/rawdata/qixiaopeng_singlecell/Result-X101SC22083646-Z01-J001-B1-1.20220920/Result-X101SC22083646-Z01-J001-B1-1/2.Summary/N_35/filtered_feature_bc_matrix"
WT="/media/liuhb/cell/liuhb/R_data/qixiaopeng/rawdata/qixiaopeng_singlecell/Result-X101SC22083646-Z01-J001-B1-1.20220920/Result-X101SC22083646-Z01-J001-B1-1/2.Summary/W_35/filtered_feature_bc_matrix"
Ncf4_KO.data=Read10X(data.dir = Ncf4_KO)
WT.data=Read10X(data.dir = WT)


Ncf4_KO <- CreateSeuratObject(counts = Ncf4_KO.data,min.cells = 3, min.features = 200,project = "Ncf4_KO")
WT<- CreateSeuratObject(counts = WT.data,min.cells = 3, min.features = 200,project = "WT")


WT[["percent.mt"]] <- PercentageFeatureSet(WT, pattern = "^mt-") ##calculate the ratio of mitochondria genes
WT[["percent.RP"]] <- PercentageFeatureSet(WT, pattern = "^Rp[sl]") ##calculate the ratio of ribosome genes

HB.genes <- c("Hba1","Hba2","Hbb","Hbd","Hbe1","Hbg1","Hbg2","Hbm","Hbq1","Hbz") ##calculate the ratio of RBC genes
library(homologene)
genelist<- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ") ##calculate the ratio of RBC genes
MsB=homologene(genelist, inTax = 9606, outTax = 10090)##homologous change
HB_m <- match(MsB$`10090`, rownames(WT@assays$RNA))  
HB.genes <- rownames(WT@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
WT[["percent.HB"]]<-PercentageFeatureSet(WT, features=HB.genes) 
VlnPlot(object=WT, features = c('nCount_RNA',"nFeature_RNA","percent.mt","percent.RP","percent.HB"),group.by = "orig.ident")
WT <- subset(WT, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 70)

###Normalized the data
WT <- NormalizeData(WT, normalization.method = "LogNormalize", scale.factor = 10000)
WT <- FindVariableFeatures(WT, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(WT)
WT <- ScaleData(WT, features = all.genes)
WT <- RunPCA(WT, features = VariableFeatures(object = WT))
WT <- RunUMAP(WT, dims = 1:30)

nExp <- round(ncol(WT) * 0.04)  ### expect 4% doublets
### remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
library(DoubletFinder)
#这个DoubletFinder包的输入是经过预处理（包括归一化、降维，但不一定要聚类）的 Seurat 对象
WT <- doubletFinder_v3(WT, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:30) ##update newly

###### 找Doublets  
### DF的名字是不固定的，因此从sce.all.filt@meta.data列名中提取比较保险
DF.name = colnames(WT@meta.data)[grepl("DF.classification", colnames(WT@meta.data))]
#(图右)
p5.dimplot=cowplot::plot_grid(ncol = 2, DimPlot(WT, group.by = "orig.ident") + NoAxes(), 
                              DimPlot(WT, group.by = DF.name) + NoAxes())
#ggsave(filename="doublet_dimplot.pdf",plot=p5.dimplot)
#ggsave(filename="doublet_dimplot.png",plot=p5.dimplot)
#(图左)红色部分上的黑点就是doublets
p5.vlnplot=VlnPlot(WT, features = "nFeature_RNA", group.by = DF.name, pt.size = 0.1)
#ggsave(filename="doublet_vlnplot.pdf",plot=p5.vlnplot)
#ggsave(filename="doublet_vlnplot.png",plot=p5.vlnplot)
###### 过滤doublet
WT=WT[, WT@meta.data[, DF.name] == "Singlet"]


###for Ncf4KO
Ncf4_KO[["percent.mt"]] <- PercentageFeatureSet(Ncf4_KO, pattern = "^mt-") ##calculate the ratio of mitochondria genes
Ncf4_KO[["percent.RP"]] <- PercentageFeatureSet(Ncf4_KO, pattern = "^Rp[sl]") ##calculate the ratio of ribosome genes
#HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ") ##calculate the ratio of RBC genes
HB.genes <- c("Hba1","Hba2","Hbb","Hbd","Hbe1","Hbg1","Hbg2","Hbm","Hbq1","Hbz") ##calculate the ratio of RBC genes
library(homologene)
genelist<- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ") ##calculate the ratio of RBC genes
MsB=homologene(genelist, inTax = 9606, outTax = 10090)##homologous change
HB_m <- match(MsB$`10090`, rownames(Ncf4_KO@assays$RNA))  
HB.genes <- rownames(Ncf4_KO@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
Ncf4_KO[["percent.HB"]]<-PercentageFeatureSet(Ncf4_KO, features=HB.genes) 
VlnPlot(object=Ncf4_KO, features = c('nCount_RNA',"nFeature_RNA","percent.mt","percent.RP","percent.HB"),group.by = "orig.ident")
Ncf4_KO <- subset(Ncf4_KO, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 70)

###Normalized the data
Ncf4_KO <- NormalizeData(Ncf4_KO, normalization.method = "LogNormalize", scale.factor = 10000)
Ncf4_KO <- FindVariableFeatures(Ncf4_KO, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Ncf4_KO)
Ncf4_KO <- ScaleData(Ncf4_KO, features = all.genes)
Ncf4_KO <- RunPCA(Ncf4_KO, features = VariableFeatures(object = Ncf4_KO))
Ncf4_KO <- RunUMAP(Ncf4_KO, dims = 1:30)

nExp <- round(ncol(Ncf4_KO) * 0.04)  ### expect 4% doublets
### remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
library(DoubletFinder)
#这个DoubletFinder包的输入是经过预处理（包括归一化、降维，但不一定要聚类）的 Seurat 对象
Ncf4_KO <- doubletFinder_v3(Ncf4_KO, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:30) ##update newly

###### 找Doublets  

DF.name = colnames(Ncf4_KO@meta.data)[grepl("DF.classification", colnames(Ncf4_KO@meta.data))]
#(图右)
p5.dimplot=cowplot::plot_grid(ncol = 2, DimPlot(Ncf4_KO, group.by = "orig.ident") + NoAxes(), 
                              DimPlot(Ncf4_KO, group.by = DF.name) + NoAxes())

p5.vlnplot=VlnPlot(Ncf4_KO, features = "nFeature_RNA", group.by = DF.name, pt.size = 0.1)

###### 过滤doublet
Ncf4_KO=Ncf4_KO[, Ncf4_KO@meta.data[, DF.name] == "Singlet"]


######merge WT and NCF4 KO
GastroSR=merge(Ncf4_KO,y=WT, add.cell.ids=c("Ncf4_KO","WT"),merge.data=TRUE)
############

GastroSR[["percent.mt"]] <- PercentageFeatureSet(GastroSR, pattern = "^mt-") ##calculate the ratio of mitochondria genes
GastroSR[["percent.RP"]] <- PercentageFeatureSet(GastroSR, pattern = "^Rp[sl]") ##calculate the ratio of ribosome genes
#HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ") ##calculate the ratio of RBC genes
HB.genes <- c("Hba1","Hba2","Hbb","Hbd","Hbe1","Hbg1","Hbg2","Hbm","Hbq1","Hbz") ##calculate the ratio of RBC genes
library(homologene)
genelist<- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ") ##calculate the ratio of RBC genes
MsB=homologene(genelist, inTax = 9606, outTax = 10090)##homologous change
HB_m <- match(MsB$`10090`, rownames(GastroSR@assays$RNA))  
HB.genes <- rownames(GastroSR@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
GastroSR[["percent.HB"]]<-PercentageFeatureSet(GastroSR, features=HB.genes) 
VlnPlot(object=GastroSR, features = c('nCount_RNA',"nFeature_RNA","percent.mt","percent.RP","percent.HB"),group.by = "orig.ident")
GastroSR <- subset(GastroSR, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 70)

###Normalized the data
GastroSR <- NormalizeData(GastroSR, normalization.method = "LogNormalize", scale.factor = 10000)
GastroSR <- FindVariableFeatures(GastroSR, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(GastroSR), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(GastroSR)
plot1
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

###
all.genes <- rownames(GastroSR)
GastroSR <- ScaleData(GastroSR, features = all.genes)
GastroSR <- RunPCA(GastroSR, features = VariableFeatures(object = GastroSR))
print(GastroSR[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(GastroSR, dims = 1:2, reduction = "pca")
DimPlot(GastroSR, reduction = "pca")
DimHeatmap(GastroSR, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(GastroSR, dims = 1:20, cells = 500, balanced = TRUE)
GastroSR <- JackStraw(GastroSR, num.replicate = 100)
GastroSR <- ScoreJackStraw(GastroSR, dims = 1:10)
JackStrawPlot(GastroSR, dims = 1:10)
ElbowPlot(GastroSR)

GastroSR <- FindNeighbors(GastroSR, dims = 1:30)
GastroSR <- FindClusters(GastroSR, resolution = 0.7)
head(Idents(GastroSR), 5)
head(GastroSR@meta.data)
table(GastroSR@meta.data$seurat_clusters)

GastroSR1=GastroSR
GastroSR1 <- RunUMAP(GastroSR1, dims = 1:30)
GastroSR1=RunTSNE(GastroSR1, dims = 1:30)
##########################

###########################
library(RColorBrewer)
cell_type_cols <- c(brewer.pal(9, "Set1"), 
                    "#FF34B3","#BC8F8F","#20B2AA","#00F5FF","#FFA500","#ADFF2F",
                    "#FF6A6A","#7FFFD4", "#AB82FF","#90EE90","#00CD00","#008B8B",
                    "#6495ED","#FFC1C1","#CD5C5C","#8B008B","#FF3030", "#7CFC00")  
DimPlot(GastroSR1, reduction = "umap",label = T,label.size = 5,cols =cell_type_cols)
FeaturePlot(GastroSR1,features = c("Ncf4"),reduction = "umap",pt.size = 0.5,split.by = "orig.ident")
VlnPlot(GastroSR1,features = c("Ncf4"),pt.size = 0.5,split.by = "orig.ident")
DimPlot(GastroSR1, reduction = "umap",label = T,label.size = 5,split.by = "orig.ident",cols =cell_type_cols)
Idents(GastroSR1)=GastroSR1$seurat_clusters

##cell type annotation

StromalCell=c('Dcn','Col1a2','Col3a1') #will not included 20230320
FeaturePlot(GastroSR1,features = StromalCell)


macrophage=c("Il1b","Cd68",'Tyrobp',"Aif1","Lyz2","Mpeg1","Pld4")
FeaturePlot(GastroSR1,features = macrophage)

B_cell=c('Cd79a','Cd79b','Cd19','Ms4a1')
FeaturePlot(GastroSR1,features = B_cell)


Plasma=c('Jchain')
FeaturePlot(GastroSR1,features =Plasma)

T_cell=c('Cd3d','Cd3e','Cd3g','Nkg7',"Klrd1","Ly6a") ##T_NK_ILCs

NK=c('Nkg7',"Klrd1")
FeaturePlot(GastroSR1,features = T_cell)


CD8_T_cell=c("Cd8a","Cd3g","Cd3d","Gzmb")
FeaturePlot(GastroSR1,features = CD8_T_cell)

CD4_T_cell=c("Cd4","Cd3g","Cd3d")
FeaturePlot(GastroSR1,features = CD4_T_cell)

NK=c('Nkg7',"Klrd1","Il2rb")
FeaturePlot(GastroSR1,features = NK)

EEC = c('Chga','Chgb','Pcsk1n')
FeaturePlot(GastroSR1,features = EEC)


Tuft= c('Trpm5','Dclk1')
FeaturePlot(GastroSR1,features = Tuft)


TA=c('Mki67','Top2a', 'Ccna2')
FeaturePlot(GastroSR1,features = TA,label.size = 5)

###****
TA_carcinoma=c("Slc12a2","Rgmb","Ephb2","Lrig1")
FeaturePlot(GastroSR1,features = Stem_cell)
###

Goblet =c('Tff3','Fcgbp','Muc2','Agr2',"Clca1")
FeaturePlot(GastroSR1,features = Goblet,max.cutoff = 10)
###

#**
seEPC=c("Reg3g","Hmgcs2","Fabp2","Sox9") #Stem-Early enterocyte precursor cell 
FeaturePlot(GastroSR1,features = seEPC)

Distal_Enterocytes_1=c("Car4","Phgr1","Slc26a3","Cdhr5")
FeaturePlot(GastroSR1,features = Distal_Enterocytes_1)

Distal_Enterocytes_2=c("Car4","Ly6g","Slc26a2")
FeaturePlot(GastroSR1,features = Distal_Enterocytes_2)

Proximal_Enterocytes=c("Fabp2","Cyp2c55","Emp1")
FeaturePlot(GastroSR1,features = Proximal_Enterocytes)

#**
Lig_Enterocytes=c("Lypd8","Prdx6","Phgr1") #Large intestine gland-Enterocytes
FeaturePlot(GastroSR1,features = Lig_Enterocytes)

###
mast_cells=c("Mcpt2","Cpa3","Cd34","Ms4a2")
FeaturePlot(GastroSR1,features = mast_cells)

####Krt+ epithelia cells
krt_epithelia=c("Krt14","Krt5","Krt6a","Krt13")
FeaturePlot(GastroSR1,features = krt_epithelia)

# ASC cell – adenoma specific cells, colonic stem and progenitor cells, 
ASC= c("Nkd1","Wif1","Notum","Axin2")
FeaturePlot(GastroSR1,features = ASC)

##newly added : to exchange TA2
Enterocyte_carcinoma=c("Slc12a2","Nop56","Fermt1")
FeaturePlot(GastroSR1,features = Enterocyte_carcinoma)

#GastroSR1_backup=GastroSR1
#GastroSR1=GastroSR1_backup
#########
C0="Enterocytes carcinoma" #"TA2"   #"Plastical Enterocytes" #"B-cells"
C1="Distal Enterocytes 2"   #"Distal Enterocytes 1" #"Stem_cell"
C2="Plastical Enterocytes"   #"Goblet" #"T-cells"
C3="Distal Enterocytes 1"   #"B-cells" #
C4="B-cells"   #"seEPC" #"Stem cell" 
C5="CD4+ T-cells"   #"CD4+ T-cells" #"T_NK_ILCs" #"Stem_cell"
C6="seEPC"   #"Distal Enterocytes 2" #"Goblet"
C7="TA carcinoma"   #"TA2" #"Goblet" #"Stem cell" #"Enterocytes"
C8="TA"  #"Lig Enterocytes" #"TA"
C9="NK cells"   #"Stem cells" #
C10="Goblet"  #"TA1" ##Goblet"    
C11="Goblet"  #"NK cells" #"T_NK_ILCs" #"Paneth"
C12="Lig Enterocytes"  #"Plasma cells"
C13="Plasma cells"  #"CD8+ T-cells" #"T_NK_ILCs" #"EEC"
C14="CD8+ T-cells"  #"Macrophage"  #
C15="Macrophage"  #"Proximal Enterocytes" #"Stem_cell"
C16="Distal Enterocytes 1"  #"Macrophage" #
C17="Macrophage"  #"Krt+ epithelia cells"
C18="B-cells"  #"B-cells" #"T-cells"
C19="Tuft"  #"Tuft" #T-cells"
C20="ASC cells"  #"EEC" #
C21="Mast cells"  #"ASC cells" #"others"
C22="Krt+ epithelia cells"  #"Mast cells" # newly update
C23="EEC"  #"Stromal cells" #"Tuft" #T-cells

GastroSR1 <- SetIdent(GastroSR1, value = "seurat_clusters")
new.cluster.ids=c(C0,C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,C13,C14,C15,C16,C17,C18,C19,C20,C21,C22,C23)
names(new.cluster.ids) <- levels(GastroSR1)
GastroSR1 <- RenameIdents(GastroSR1, new.cluster.ids)

fact=c("TA","TA carcinoma","Enterocytes carcinoma","Goblet","Plastical Enterocytes","Distal Enterocytes 1","Distal Enterocytes 2","Lig Enterocytes","seEPC","Krt+ epithelia cells","ASC cells","EEC","Tuft","Macrophage","Plasma cells","B-cells","CD8+ T-cells","CD4+ T-cells","NK cells","Mast cells")  
GastroSR1@active.ident=factor(GastroSR1@active.ident,levels=fact)
GastroSR1$orig.ident

library(RColorBrewer)
cell_type_cols=c("#00BFFF","#FFA500","#DDA0DD","#FF69B4","#6A5ACD","#9932CC",
                 "#32CD32","#DB7093","#808080","#8FBC8F","#A0522D","#C71585",
                 "#87CEFA","#00CED1","#8B0000","#BC8F8F","#F08080","#DAA520",
                 "#90EE90","#3CB371","#DC143C","#191970")

GastroSR1$orig.ident=factor(GastroSR1$orig.ident,levels = c("WT","Ncf4_KO"))
DimPlot(GastroSR1,reduction = "umap",label = F,pt.size = 0.2,label.size = 4,split.by = "orig.ident",repel=T,cols =cell_type_cols)
DimPlot(GastroSR1,reduction = "umap",label = T,pt.size = 0.2,label.size = 4,repel=T,cols =cell_type_cols)



#####
##downsample 
#####
set.seed(123)
WD_WT_down=sample(colnames(GastroSR1)[GastroSR1$orig.ident=="Ncf4_KO"],length((GastroSR1$orig.ident)[GastroSR1$orig.ident=="WT"]),replace = F)
WT_down=colnames(GastroSR1)[GastroSR1$orig.ident=="WT"]
downsampled.cells=c(WD_WT_down,WT_down)
GastroSR1.downsampled.cells=GastroSR1[,downsampled.cells]
DimPlot(GastroSR1.downsampled.cells, reduction = "umap",repel = "TRUE",label = T,label.size = 4,cols =cell_type_cols,split.by = "orig.ident")
DimPlot(GastroSR1.downsampled.cells, reduction = "umap",repel = "TRUE",label = F,label.size = 4,cols =cell_type_cols,split.by = "orig.ident")
####
GastroSR1.downsampled.cells$celltype=GastroSR1.downsampled.cells@active.ident
####

library(scRNAtoolVis)
GastroSR1.downsampled.cells$cell_type=GastroSR1.downsampled.cells@active.ident
clusterCornerAxes(object = GastroSR1.downsampled.cells,reduction = 'umap',
                  noSplit = FALSE, groupFacet = "orig.ident",clusterCol = "cell_type",pSize=0.4,cellLabel=T,cellLabelSize=3.5,relLength = 0.3,axes = "one",arrowType="closed")+scale_color_manual(values = cell_type_cols) 

clusterCornerAxes(object = GastroSR1.downsampled.cells,reduction = 'umap',
                  noSplit = FALSE, groupFacet = "orig.ident",clusterCol = "cell_type",pSize=0.05,cellLabel=F,cellLabelSize=3.5,relLength = 0.3,axes = "one",arrowType="closed")+scale_color_manual(values = cell_type_cols) 

clusterCornerAxes(object = GastroSR1.downsampled.cells,reduction = 'umap',
                  aspect.ratio = 1,keySize=7,base_size=18,themebg = 'bwCorner',noSplit = FALSE, groupFacet = "orig.ident",clusterCol = "cell_type",pSize=0.05,cellLabel=F,cellLabelSize=3,relLength = 0.3,arrowType="closed")+scale_color_manual(values = cell_type_cols) 

ggsave("p4.pdf",width = 15,height = 6)
dev.off()


#GastroSR1=GastroSR1.downsampled.cells   ##newly added,2023,1,31  ##cell number; WT vs Ncf4_KO: 9127 vs 12069 
##calculator the cell type ratio
WT_35=subset(GastroSR1.downsampled.cells,subset=orig.ident=="WT")
NCF4_KO=subset(GastroSR1.downsampled.cells,subset=orig.ident=="Ncf4_KO")

WT=table(WT_35@active.ident)/sum(table(WT_35@active.ident))
Ncf4_KO=table(NCF4_KO@active.ident)/sum(table(NCF4_KO@active.ident))
#Ncf4_KO=c(Ncf4_KO,0)
#names(Ncf4_KO)[length(Ncf4_KO)]="Krt+ epithelia cells"
seqx=match("Krt+ epithelia cells",row.names(WT))
Ncf4_KO=c(Ncf4_KO[1:seqx-1],0,Ncf4_KO[(seqx):length(Ncf4_KO)])
names(Ncf4_KO)[seqx]="Krt+ epithelia cells"

dat=cbind(WT,Ncf4_KO)
dat=dat[rev(rownames(dat)),]
library(WGCNA)
barplot(dat,beside = F,col=rev(cell_type_cols[1:length(rownames(dat))]),xlim = c(0,5),ylim=c(0,1),legend=T,border = F) 

par(mfrow=c(1,2))
barplot(-dat[,1]*100,axes=F,horiz=T,axisnames=FALSE, space = 0,
        xlim=c(-22,0),col=rev(cell_type_cols[1:length(rownames(dat))]) )  ##WT
text(round(100*dat[,1],1),x = -100*dat[,1]-2,y = c(1:length(rownames(dat)))-0.5)
#creating a new axis with desired labels
axis(side=1,at=seq(-30,0,10),labels=c(30,20,10,0))
####
barplot(100*dat[,2],horiz=T,axes=F,las=1,xlim=c(0,22),axisnames=FALSE,space = 0,col=rev(cell_type_cols[1:length(rownames(dat))])) ##NCF4 KO
text(round(100*dat[,2],1),x = 100*dat[,2]+2,y = c(1:length(rownames(dat)))-0.5)
axis(side=1,at=seq(0,30,10),labels=c(0,10,20,30))





####
marker_gene=c(TA,TA_carcinoma,Enterocyte_carcinoma,Goblet,Plastical_Enterocytes,Distal_Enterocytes_1,Distal_Enterocytes_2,Lig_Enterocytes,seEPC,krt_epithelia,ASC,StromalCell,EEC,Tuft,macrophage,Plasma,B_cell,CD8_T_cell,CD4_T_cell,NK,mast_cells)
dat=list(TA,TA_carcinoma,Enterocyte_carcinoma,Goblet,Plastical_Enterocytes,Distal_Enterocytes_1,Distal_Enterocytes_2,Lig_Enterocytes,seEPC,krt_epithelia,ASC,StromalCell,EEC,Tuft,macrophage,Plasma,B_cell,CD8_T_cell,CD4_T_cell,NK,mast_cells)
cluster_gene=NULL
for (i in 1:length(dat)){
  a=rep(i,length(dat[[i]]))
  cluster_gene=c(cluster_gene,a)
}
#cluster_gene
dat1=data.frame(marker_gene,cluster_gene)
dat1=dat1[!duplicated(dat1$marker_gene),]
dat1=dat1[order(dat1$cluster_gene,decreasing = T),]
AverageHeatmap(object = GastroSR1.downsampled.cells,column_split = 1:length(unique(GastroSR1.downsampled.cells$cell_type)),
               markerGene = dat1$marker_gene,annoCol=T,myanCol = cell_type_cols,border = T,
               row_split=dat1$cluster_gene,fontsize=9) 
ggsave("p4.pdf",width = 20,height = 10)
dev.off()
#AverageHeatmap(object = GastroSR1.downsampled.cells,column_split = 1:21,markerGene = marker_gene,annoCol=T,myanCol = cell_type_cols,border = T) 

p1=jjDotPlot(object = GastroSR1.downsampled.cells,gene = dat1$marker_gene,id = 'cell_type',rescale = T,rescale.min = 0,rescale.max = 1,split.by = 'orig.ident')
p1=jjDotPlot(object = GastroSR1.downsampled.cells,gene = dat1$marker_gene,id = 'cell_type',rescale = T,ytree = F)
p1
ggsave("p1.pdf",width = 17,height =17)
dev.off()


######
library(dplyr)
all.markers <- FindAllMarkers(GastroSR1.downsampled.cells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
head(all.markers, n = 20)

library("clusterProfiler")
library("org.Mm.eg.db")
all.markers_filter=all.markers[all.markers$avg_log2FC>0.5,]
#write.csv(all.markers,"all.markers_bycellcluster_20230510.csv",row.names = T)
geneanno=list()
genename=NULL
fact=levels(GastroSR1.downsampled.cells$celltype)
cluster=fact[c(10)] #fact[c(5,18)] #"all" #fact[1:6]   #all , #1, #2, Name of the file keep as similar to the next line;
for (i in fact[c(10)]) {
  ID_gene=all.markers_filter$gene[all.markers_filter$cluster==i]
  genename=c(genename,ID_gene)
  geneanno[[i]]=bitr(geneID=ID_gene, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb='org.Mm.eg.db', drop = TRUE)
  geneanno[[i]]=geneanno[[i]]$ENTREZID
}

geneGO <- compareCluster(geneCluster   = geneanno,
                         fun='enrichGO', OrgDb='org.Mm.eg.db',pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH",ont="BP")
write.csv(geneGO,paste(paste(cluster,collapse = "_"),'_marker_geneGO.csv'),row.names = F)

geneKEGG <- compareCluster(geneCluster   = geneanno,
                           fun           = "enrichKEGG",organism= 'mouse',
                           pvalueCutoff  = 0.05,pAdjustMethod = "BH")
write.csv(geneGO,paste(paste(cluster,collapse = "_"),'_marker_geneKEGG.csv'),row.names = F)

geneReactomePA<- compareCluster(geneCluster   = geneanno, fun="enrichPathway",
                                organism = "mouse",
                                pvalueCutoff=0.05,pAdjustMethod = "BH")
write.csv(geneGO,paste(paste(cluster,collapse = "_"),'_marker_geneReactomePA.csv'),row.names = F)

#as.data.frame(geneReacom)  ##output with the results;

#visulized 
#setwd("/Volumes/GoogleDrive/My Drive/data_BCH/R_data/Gastro_qixp/20221223/marker_gene")
# Run GO enrichment test and merge terms 
# that are close to each other to remove result redundancy
lineage1_ego <- simplify(
  geneGO, 
  cutoff=0.7, 
  by="p.adjust", 
  select_fun=min
) 
lineage1_ego@compareClusterResult$p.adjust=-log10(lineage1_ego@compareClusterResult$p.adjust)
max_GO=quantile(lineage1_ego@compareClusterResult$p.adjust,0.95)
geneKEGG@compareClusterResult$p.adjust=-log10(geneKEGG@compareClusterResult$p.adjust)
max_KEGG=quantile(geneKEGG@compareClusterResult$p.adjust,0.95)
geneReactomePA@compareClusterResult$p.adjust=-log10(geneReactomePA@compareClusterResult$p.adjust)
max_ReactomePA=quantile(geneReactomePA@compareClusterResult$p.adjust,0.95)


p2=dotplot(lineage1_ego,title="GO",showCategory=20)+scale_color_gradient(limits=c(0,max_GO),low="blue",high = "red",na.value = "red",name="-Log10 (p.adjust)")+theme_bw()+theme(axis.text.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5,size = 11), 
                                                                                                                                                                                axis.text.y = element_text(size = 11),plot.title = element_text(size=12,hjust=0.5))+ xlab("") + scale_y_discrete(labels = function(x) str_wrap(x, width = 80))  
ggsave(paste(paste(cluster,collapse = "_"),'_marker_geneGO.pdf'),width = 15,height = 20)

p3=dotplot(geneKEGG,title="KEGG",showCategory=20)+scale_color_gradient(limits=c(0,max_KEGG),low="blue",high = "red",na.value = "red",name="-Log10 (p.adjust)")+theme_bw()+theme(axis.text.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5,size = 11),
                                                                                                                                                                                axis.text.y = element_text(size = 11),plot.title = element_text(size=12,hjust=0.5))+ xlab("") +scale_y_discrete(labels = function(x) str_wrap(x, width = 80))
ggsave(paste(paste(cluster,collapse = "_"),'_marker_geneKEGG.pdf'),width = 15,height = 20)

p4=dotplot(geneReactomePA,title="ReactomePA",showCategory=20)+scale_color_gradient(limits=c(0,max_ReactomePA),low="blue",high = "red",na.value = "red",name="-Log10 (p.adjust)")+theme_bw()+theme(axis.text.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5,,size = 11),
                                                                                                                                                                                                  axis.text.y = element_text(size = 11),plot.title = element_text(size=12,hjust=0.5))+ xlab("") +scale_y_discrete(labels = function(x) str_wrap(x, width = 80))
ggsave(paste(paste(cluster,collapse = "_"),'_marker_ReactomePA.pdf'),width = 15,height = 20)
##showCategory=c("ATP metabolic process") 

AverageHeatmap(object = aor2,markerGene = genename,annoCol=T,myanCol = cell_type_cols[1:12],showRowNames=F) 





########################################## GO for each cell type
#####barplot to show DEGs between cluster I vs cluter 2,3 or 4......
library("clusterProfiler")
library("org.Mm.eg.db")
library("ReactomePA")
#cluster_markers <- FindMarkers(GastroSR1, ident.1 = "4",ident.2 = c("1","6","8","15"), min.pct = 0.25,only.pos = TRUE, logfc.threshold = 1)
fact=levels(GastroSR1.downsampled.cells$celltype)
iden1=fact[10]
iden2=fact[-10]
min.pct=0.25
logfc.threshold=0.5
cluster_markers <- FindMarkers(GastroSR1.downsampled.cells, ident.1 = iden1, min.pct = min.pct,only.pos = TRUE, logfc.threshold = logfc.threshold)
#cluster_markers <- FindMarkers(GastroSR1.downsampled.cells, ident.1 = iden1,ident.2 =iden2, min.pct = min.pct,only.pos = TRUE, logfc.threshold = logfc.threshold)
geneanno1=list()
geneanno2=list()
geneanno1[[1]]=row.names(cluster_markers)
ID_gene=geneanno1[[1]]
geneanno2[[1]]=bitr(geneID=ID_gene, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb='org.Mm.eg.db', drop = TRUE)
geneanno2[[1]]=geneanno2[[1]]$ENTREZID

geneGO=enrichGO(as.numeric(geneanno2[[1]]),OrgDb='org.Mm.eg.db', ont="BP", pvalueCutoff=0.01)
lineage1_ego <- simplify(
  geneGO, 
  cutoff=0.7, 
  by="p.adjust", 
  select_fun=min
) 
geneKEGG=enrichKEGG(as.numeric(geneanno2[[1]]),organism = "mouse", pvalueCutoff=0.01)
geneReactomePA=enrichPathway(as.numeric(geneanno2[[1]]), pvalueCutoff=0.01,organism = "mouse",)

lineage1_ego@result$p.adjust=log10(lineage1_ego@result$p.adjust)
min_GO=quantile(lineage1_ego@result$p.adjust,0.05)   ##scale_fill_gradient(limits=c(min_GO,max_GO)
max_GO=quantile(lineage1_ego@result$p.adjust,0.95)
lineage1_ego@result$p.adjust[lineage1_ego@result$p.adjust<min_GO]=min_GO
lineage1_ego@result$p.adjust[lineage1_ego@result$p.adjust>max_GO]=max_GO

showCategory=c("epidermis development","wound healing","keratinocyte differentiation","intermediate filament-based process","skin development","epithelial cell proliferation","regulation of response to wounding","regulation of epithelial cell proliferation","epithelial cell migration","cell-cell junction organization")
t1=barplot(lineage1_ego,x = "GeneRatio",color = "p.adjust",showCategory = showCategory,title = paste("GO",iden1,"maker gene",sep = "_"))+ scale_fill_gradient(limits=c(-10,0),low="red",high = "blue",name="Log10 (p.adjust)")+scale_y_discrete(labels = function(x) str_wrap(x, width = 80))+
  theme(axis.text.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5,size = 11), axis.text.y = element_text(size = 10),plot.title = element_text(size=12,hjust=0.5))  

t2=barplot(geneKEGG,x = "GeneRatio",showCategory = 100,title = paste("KEGG",iden1,"vs",paste(iden2,collapse = "_"),sep = "_"))+ scale_y_discrete(labels = function(x) str_wrap(x, width = 80))+ 
  theme(axis.text.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5,size = 11), axis.text.y = element_text(size = 10),plot.title = element_text(size=12,hjust=0.5))  

t3=barplot(geneReactomePA,x = "GeneRatio",showCategory = 50,title = paste("ReactomePA",iden1,"vs",paste(iden2,collapse = "_"),sep = "_"))+ scale_y_discrete(labels = function(x) str_wrap(x, width = 80))+ 
  theme(axis.text.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5,size = 11), axis.text.y = element_text(size = 10),plot.title = element_text(size=12,hjust=0.5))  

library(patchwork)
library(cowplot)
library("stringr")  
plot_grid(t1,t2,t3,labels=LETTERS[c(1,2,3)], ncol=3)
ggsave(paste("Functional enrichments",iden1,"vs",paste(iden2,collapse = "_"),"min.pct",min.pct,"logfc.threshold",logfc.threshold,'merge_enrichments.pdf',sep = "_"),width = 25,height = 10)




########subtype of immune cells
immu=subset(GastroSR1.downsampled.cells,idents = c("CD8+ T-cells","CD4+ T-cells","NK cells"))
immu$cell_type=immu@active.ident

immu <- NormalizeData(immu, normalization.method = "LogNormalize", scale.factor = 10000)
immu <- FindVariableFeatures(immu, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(immu)
immu <- ScaleData(immu, features = all.genes)
immu <- RunPCA(immu, features = VariableFeatures(object = immu))
#immu <- ScoreJackStraw(immu, dims = 1:5)
#JackStrawPlot(immu, dims = 1:5)
#ElbowPlot(immu)

immu <- FindNeighbors(immu, dims = 1:10)
immu <- FindClusters(immu, resolution = 0.8)  #0.5
immu <- RunUMAP(immu, dims = 1:10)
DimPlot(immu,reduction = 'umap',cols = cell_type_cols,label = T)

immu=subset(immu,idents = c(0:9))
###
C0="C0"
C1="C5"
C2="C2"
C3="C3"
C4="C4"
C5="C1"
C6="C6"
C7="C0"
C8="C5"
C9="C7"
#C10="C8"

#####
immu <- SetIdent(immu, value = "RNA_snn_res.0.8")
new.cluster.ids=""

new.cluster.ids=c(C0,C1,C2,C3,C4,C5,C6,C7,C8,C9)
names(new.cluster.ids) <- levels(immu)
immu <- RenameIdents(immu, new.cluster.ids)


fact=c("C0","C1","C2","C3","C4","C5","C6","C7")
immu@active.ident=factor(immu@active.ident,levels=fact)
immu$subtype=immu@active.ident

p1=DimPlot(immu,reduction = 'umap',cols = cell_type_cols,label = T,label.size = 5) & xlim(c(-8,8)) & ylim(c(-5,5))
p2=DimPlot(immu,reduction = 'umap',group.by  = "cell_type",cols = cell_type_cols[17:19]) & xlim(c(-8,8)) & ylim(c(-5,5))
p3=DimPlot(immu,reduction = 'umap',group.by  = "orig.ident",cols = c("red","cyan")) & xlim(c(-8,8)) & ylim(c(-5,5))
p1+p2+p3
ggsave("p5.pdf",width = 15,height =4)
dev.off()

immu$new_type=immu$orig.ident
########
immu$new_type=immu$orig.ident
WT=subset(immu,subset=new_type=="WT")
Ncf4_KO=subset(immu,subset=new_type=="Ncf4_KO")
WT=table(WT@active.ident)
Ncf4_KO=table(Ncf4_KO@active.ident)
dat=rbind(WT,Ncf4_KO)
#dat=(dat/rep(rowSums(dat),ncol(dat)))*10000 ##normalized with total cells numbers in each samples;
dat=dat/rbind((colSums(dat)),(colSums(dat)))  ##for three samples
dat=dat[rev(rownames(dat)),]
dat=dat*100
library(WGCNA)
par(las = 2)
barplot(dat,beside = F,col=rev(c("red","cyan")),xlim = c(0,13),ylim=c(0,100),legend=F,border = F,ylab = "Percentage of cells (%)",las=1,xlab = "sub_clusters") 
legend("topright",col=c("red","cyan"),c("WT","Ncf4_KO"),pch = 15,box.lty=0 )

############
all.markers <- FindAllMarkers(immu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
head(all.markers, n = 20)
##show the top two genes in each cluster according to avg_log2FC;
Top=all.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
write.csv(all.markers,"immue_all.markers_bycelltype.csv",row.names = T)
###
library("clusterProfiler")
library("org.Mm.eg.db")
all.markers_filter=all.markers[all.markers$avg_log2FC>0.5,]
#write.csv(all.markers,"all.markers_bycellcluster.csv",row.names = T)
geneanno=list()
genename=NULL
fact=levels(immu@active.ident)
cluster=fact #fact[c(5,18)] #"all" #fact[1:6]   #all , #1, #2, Name of the file keep as similar to the next line;
for (i in fact) {
  ID_gene=all.markers_filter$gene[all.markers_filter$cluster==i]
  genename=c(genename,ID_gene)
  geneanno[[i]]=bitr(geneID=ID_gene, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb='org.Mm.eg.db', drop = TRUE)
  geneanno[[i]]=geneanno[[i]]$ENTREZID
}

geneGO <- compareCluster(geneCluster   = geneanno,
                         fun='enrichGO', OrgDb='org.Mm.eg.db',pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH",ont="BP")
#write.csv(geneGO,paste(paste(cluster,collapse = "_"),'_marker_geneGO.csv'),row.names = F)

geneKEGG <- compareCluster(geneCluster   = geneanno,
                           fun           = "enrichKEGG",organism= 'mouse',
                           pvalueCutoff  = 0.05,pAdjustMethod = "BH")
#write.csv(geneGO,paste(paste(cluster,collapse = "_"),'_marker_geneKEGG.csv'),row.names = F)

geneReactomePA<- compareCluster(geneCluster   = geneanno, fun="enrichPathway",
                                organism = "mouse",
                                pvalueCutoff=0.05,pAdjustMethod = "BH")
#write.csv(geneGO,paste(paste(cluster,collapse = "_"),'_marker_geneReactomePA.csv'),row.names = F)

#visulized 
#setwd("/Volumes/GoogleDrive/My Drive/data_BCH/R_data/Gastro_qixp/20221223/marker_gene")
# Run GO enrichment test and merge terms 
# that are close to each other to remove result redundancy
lineage1_ego <- clusterProfiler::simplify(
  geneGO,
  cutoff = 0.7, 
  by = "p.adjust", 
  select_fun = min
) 

lineage1_ego@compareClusterResult$p.adjust=log10(lineage1_ego@compareClusterResult$p.adjust)
min_GO=quantile(lineage1_ego@compareClusterResult$p.adjust,0.05)
max_GO=quantile(lineage1_ego@compareClusterResult$p.adjust,0.95)
lineage1_ego@compareClusterResult$p.adjust[lineage1_ego@compareClusterResult$p.adjust<min_GO]=min_GO
lineage1_ego@compareClusterResult$p.adjust[lineage1_ego@compareClusterResult$p.adjust>max_GO]=max_GO
p2=dotplot(lineage1_ego,title="GO",showCategory=7)+scale_color_gradient(limits=c(min_GO,0),low="red",high = "blue",na.value = "red",name="Log10 (p.adjust)")+theme_bw()+theme(axis.text.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5,size = 11), 
                                                                                                                                                                              axis.text.y = element_text(size = 11),plot.title = element_text(size=12,hjust=0.5))+ xlab("") + scale_y_discrete(labels = function(x) str_wrap(x, width = 80))  
p2
ggsave(paste(paste(cluster,collapse = "_"),'_marker_geneGO.pdf'),width = 15,height = 20)






lineage1_ego@compareClusterResult$p.adjust=-log10(lineage1_ego@compareClusterResult$p.adjust)
max_GO=quantile(lineage1_ego@compareClusterResult$p.adjust,0.95)
geneKEGG@compareClusterResult$p.adjust=-log10(geneKEGG@compareClusterResult$p.adjust)
max_KEGG=quantile(geneKEGG@compareClusterResult$p.adjust,0.95)
geneReactomePA@compareClusterResult$p.adjust=-log10(geneReactomePA@compareClusterResult$p.adjust)
max_ReactomePA=quantile(geneReactomePA@compareClusterResult$p.adjust,0.95)


p2=dotplot(lineage1_ego,title="GO",showCategory=20)+scale_color_gradient(limits=c(0,max_GO),low="blue",high = "red",na.value = "red",name="-Log10 (p.adjust)")+theme_bw()+theme(axis.text.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5,size = 11), 
                                                                                                                                                                                axis.text.y = element_text(size = 11),plot.title = element_text(size=12,hjust=0.5))+ xlab("") + scale_y_discrete(labels = function(x) str_wrap(x, width = 80))  
ggsave(paste(paste(cluster,collapse = "_"),'_marker_geneGO.pdf'),width = 15,height = 20)

p3=dotplot(geneKEGG,title="KEGG",showCategory=20)+scale_color_gradient(limits=c(0,max_KEGG),low="blue",high = "red",na.value = "red",name="-Log10 (p.adjust)")+theme_bw()+theme(axis.text.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5,size = 11),
                                                                                                                                                                                axis.text.y = element_text(size = 11),plot.title = element_text(size=12,hjust=0.5))+ xlab("") +scale_y_discrete(labels = function(x) str_wrap(x, width = 80))
ggsave(paste(paste(cluster,collapse = "_"),'_marker_geneKEGG.pdf'),width = 15,height = 20)

p4=dotplot(geneReactomePA,title="ReactomePA",showCategory=20)+scale_color_gradient(limits=c(0,max_ReactomePA),low="blue",high = "red",na.value = "red",name="-Log10 (p.adjust)")+theme_bw()+theme(axis.text.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5,,size = 11),
                                                                                                                                                                                                  axis.text.y = element_text(size = 11),plot.title = element_text(size=12,hjust=0.5))+ xlab("") +scale_y_discrete(labels = function(x) str_wrap(x, width = 80))
ggsave(paste(paste(cluster,collapse = "_"),'_marker_ReactomePA.pdf'),width = 15,height = 20)
##showCategory=c("ATP metabolic process") 

labelgene=c("Fcer1g","Xcl1","Klrd1","Cd160","Malat1","Gzmk","Gzmb","Gzma","Ifng","mt-Co1","mt-Atp6","Zg16","Izumo1r","Maf","Icos","Ctla4","Cd28","Igfbp4","Rps24","Rps20","Sell","Tmem176a","Rorc","Il17a","Slc44a4","Vegfa")

genename1=genename[-match(c("mt-Co1","mt-Atp6"),genename)] ##only show the the second;
AverageHeatmap(object = immu,markerGene = genename1,markGenes = labelgene,annoCol=T,myanCol = cell_type_cols[1:8],showRowNames=F) 
par(las=2)
barplot(table(all.markers_filter$cluster),col = cell_type_cols[1:8],space = 0,ylim = c(0,180),las=1,,ylab="cell-type specific genes")
text(table(all.markers_filter$cluster),x = c(1:length(table(all.markers_filter$cluster))-0.5),y = table(all.markers_filter$cluster)+10)

#########
Idents(immu)=immu$subtype
immu$celltype1=paste(immu$orig.ident,immu@active.ident,sep = "_")
Idents(immu)=immu$celltype1

for (k in c("C3","C5","C7")){
  cluster=k   ##only set this 
  different_express_gene<- FindMarkers(immu,
                                       logfc.threshold = 0.5,
                                       test.use = "wilcox",
                                       ident.1 = paste("Ncf4_KO",cluster,sep = "_"),
                                       ident.2 = paste("WT",cluster,sep = "_"))
  #detach(package:COSG,unload=TRUE)
  #library(COSG) ##developmed by DaiMin
  #different_express_gene1=cosg(immu,groups='all',assay='RNA',slot='data',mu=100,remove_lowly_expressed=TRUE,expressed_pct=0.1)  
  #different_express_gene2=cosg(immu,groups=c(seq(0:22)),assay='RNA',slot='data',mu=100,remove_lowly_expressed=TRUE,expressed_pct=0.1)  
  #AverageHeatmap(object = immu,markerGene = as.character(as.matrix((different_express_gene2[[1]]))),annoCol=T,myanCol = cell_type_cols[1:24],showRowNames=F,use_raster = FALSE)
  diff_gona=different_express_gene
  diff_gona=diff_gona[!(is.na(diff_gona$avg_log2FC) | is.na(diff_gona$avg_log2FC)),]
  gene_down=rownames(diff_gona[diff_gona$avg_log2FC < -0.5 & diff_gona$p_val_adj<0.001,])
  gene_up=rownames(diff_gona[diff_gona$avg_log2FC > 0.5 & diff_gona$p_val_adj<0.001,])
  keyvals <- ifelse(
    diff_gona$avg_log2FC < -0.5 & diff_gona$p_val_adj<0.001, 'royalblue',
    ifelse(diff_gona$avg_log2FC > 0.5 & diff_gona$p_val_adj<0.001, 'gold',
           'black'))
  keyvals[is.na(keyvals)] <- 'black'
  names(keyvals)[keyvals == 'gold'] <- 'Upregulated'
  names(keyvals)[keyvals == 'black'] <- 'No.sig'
  names(keyvals)[keyvals == 'royalblue'] <- 'Downregulated'
  library("EnhancedVolcano")
  p1=EnhancedVolcano (diff_gona,
                      lab=rownames(diff_gona),
                      x="avg_log2FC",
                      y="p_val_adj",
                      pCutoff=0.001,
                      FCcutoff=0.5,
                      caption=NULL,
                      title=paste("Ncf4_KO vs WT for","cluster",cluster,sep = " "),
                      subtitle = paste(length(gene_down),"vs",length(gene_up),sep = " "),
                      colCustom = keyvals,
                      pointSize = 3.5,
                      labSize = 4.5,
                      colAlpha = 0.8,
                      legendLabSize = 15,
                      legendIconSize = 5.0,
                      drawConnectors = TRUE,
                      widthConnectors = 0.5,
                      colConnectors = 'black',
                      arrowheads = FALSE,
                      border = 'partial',
                      borderWidth = 0.5,
                      borderColour = 'black',legendPosition ="right")+theme(plot.title = element_text(size=12,hjust=0.5),plot.subtitle  = element_text(size=12,hjust=0.5),legend.spacing.x = unit(0.2, "lines"),legend.box.margin=margin(0.1,0.1,0.1,-0.8,'cm'))+guides(fill = guide_legend(byrow = TRUE))
  ggsave(paste(cluster,'_NcfKO_vs_WT_Volcano.pdf'),width = 10,height = 7)
}


Idents(immu)=immu$subtype

VlnPlot(immu,features = c("mt-Co1","mt-Atp6","Rps24","Rps20"),split.by="orig.ident",ncol = 2) & ylim(c(0,10))


immu35=subset(immu,idents = c("C3"))
VlnPlot(immu35,features = c("mt-Co1","mt-Atp6"),split.by="orig.ident",ncol = 2) & ylim(c(2.5,7.5))
immu35=subset(immu,idents = c("C5"))
VlnPlot(immu35,features = c("Rps24","Rps20"),split.by="orig.ident",ncol = 2) & ylim(c(3,6))

immu35=subset(immu,idents = c("C3"))
VlnPlot(immu35,features = c("mt-Co1","mt-Atp6","mt-Co2","mt-Co3","mt-Nd4","mt-Cytb","mt-Nd4l","mt-Nd1","mt-Nd5"),split.by="orig.ident",ncol = 5)  & ylim(c(2.5,7.5))

immu35=subset(immu,idents = c("C5"))
VlnPlot(immu35,features = c("Rps24","Rps20","Rps19","Rps29","Rps3a1","Rpl35a","Rpl8","Rps7","Rps27a"),split.by="orig.ident",ncol = 5) & ylim(c(3,6))


jjDotPlot(object = immu,gene = c("mt-Co1","mt-Atp6","Rps24","Rps20"),id = 'subtype',rescale = T,ytree = F)
ggsave("p1.pdf",width = 17,height =17)
dev.off()



##T-cells;
Inhibitory_receptors=c("Ctla4","Icos","Tigit","Cd28","Havcr2","Pdcd1","Lag3")
FeaturePlot(immu,features = Inhibitory_receptors[1:4],split.by = "orig.ident")
FeaturePlot(immu,features = Inhibitory_receptors[5:7],split.by = "orig.ident",)

Effector=c("Ifng","Ccl3","Ccl4","Prf1","Gzmb","Tnfsf10","Csf1","Fasl","Ccl5","Ccl2","Tgfbi","Cxcl10","Il2","Il10")
FeaturePlot(immu,features = Effector[1:4],split.by = "orig.ident")
FeaturePlot(immu,features = Effector[5:8],split.by = "orig.ident")
FeaturePlot(immu,features = Effector[9:12],split.by = "orig.ident")
FeaturePlot(immu,features = Effector[13:14],split.by = "orig.ident")
FeaturePlot(immu,features = c("Ifng","Gzmb"),split.by = "orig.ident")& 
  theme(legend.position = "right")
#####
TFs=c("Runx1","Hif1a","Nfkb1","Satb1","Id2","Bcl6","Prdm1","Runx3","Tox","Tbx21","Eomes","Id3","Bach2","Tcf7","Nfatc1")
FeaturePlot(immu,features = TFs[1:4],split.by = "orig.ident")
FeaturePlot(immu,features = TFs[5:8],split.by = "orig.ident")
FeaturePlot(immu,features = TFs[9:12],split.by = "orig.ident")
FeaturePlot(immu,features = TFs[13:16],split.by = "orig.ident")

#####
FeaturePlot(immu,features = c("Gata3","Il4"),split.by = "orig.ident",pt.size = 0.001)
FeaturePlot(immu,features = c("Il17a"),split.by = "orig.ident",pt.size = 0.001)
FeaturePlot(immu,features = c("Foxp3"),split.by = "orig.ident",pt.size = 0.001)
FeaturePlot(immu,features = c("Cxcr5","Il21"),split.by = "orig.ident",pt.size = 0.001)
FeaturePlot(immu,features = c("Il22","Il9","Cd4"),split.by = "orig.ident",pt.size = 0.001)



clusterCornerAxes(object = immu,reduction = 'umap',
                  noSplit = FALSE, groupFacet = "orig.ident",pSize=0.4,cellLabel=T,cellLabelSize=3.5,relLength = 0.3,axes = "one",arrowType="closed")

clusterCornerAxes(object = immu,reduction = 'umap',
                  noSplit = T, groupFacet = "orig.ident",pSize=0.4,cellLabel=T,cellLabelSize=3.5,relLength = 0.3,axes = "one",arrowType="closed")

clusterCornerAxes(object = immu,reduction = 'umap',
                  noSplit = FALSE, groupFacet = "cell_type",pSize=0.4,cellLabel=T,cellLabelSize=3.5,relLength = 0.3,axes = "one",arrowType="closed")

DimPlot(immu,reduction = 'umap',group.by  = "orig.ident")

#####
mac=subset(GastroSR1.downsampled.cells,idents = c("Macrophage"))
mac_gene1=c("Cxcr1","Retnlg","Tgfbr1","Cxcl2","Ccl4","Cxcl3","Nos2","Inhba","Hilpda","Arg1","Il7r","Pdcd1lg2","Il18bp","Il18","Il18rap","Il18r1","Ifngr1","Tgfbi","Ccr2","Ly6c2","Hp","Klf2","Vcam1","Ms4a7","Cd63","Apoe","Mertk","Cd72","Birc5","Top2a","Ube2c","Retnla","Ccl8","Ccl7","Pltp","Mrc1","Cx3cr1","Cst3","Cd209a","Ccl17")
mac_gene2=c("Il15","Tnf","Cd40","Cd74","Ptgs2","Cxcl9","Cxcl10","Pf4","Apoe","Arg1","Mrc1","Il10","Spp1","Ccl9","C1qa","C1qb")
DotPlot(mac,features = unique(c(mac_gene1,mac_gene2)),group.by = "orig.ident")+ggplot2:::coord_flip()


##
##NK-cells
Inhibitory_receptors=c("Ctla4","Klra7","Klra8","Klra4","Klra3","Klrb1b","Pdcd1","Klrc1","Ltb","Klra1","Klra6")
FeaturePlot(immu,features = Inhibitory_receptors[1:4],split.by = "orig.ident")
FeaturePlot(immu,features = Inhibitory_receptors[5:8],split.by = "orig.ident",)
FeaturePlot(immu,features = Inhibitory_receptors[9:11],split.by = "orig.ident",)
p1=FeaturePlot(immu,features = Inhibitory_receptors,split.by = "orig.ident")
p1
ggsave("NK_Inhibitory_receptors.pdf",width = 10,height = 40)
dev.off()

Effector=c("Tnf","Xcl1","Ccl3","Ccl5","Ccl4","Prf1","Gzmb","Ifng","Il2")
FeaturePlot(immu,features = Effector[1:4],split.by = "orig.ident")
FeaturePlot(immu,features = Effector[5:8],split.by = "orig.ident")
p2=FeaturePlot(immu,features = Effector,split.by = "orig.ident")
ggsave("NK_Effector.pdf",width = 7,height = 28)
dev.off()
#####
TFs=c("Tcf7","Klf2","Zeb2","Irf8","Tbx21","Klf7","Prdm1","Eomes","Cbfb")
FeaturePlot(immu,features = TFs[1:4],split.by = "orig.ident")
FeaturePlot(immu,features = TFs[5:8],split.by = "orig.ident")
p3=FeaturePlot(immu,features = TFs,split.by = "orig.ident")
p3
ggsave("NK_TF.pdf",width = 7,height = 28)
dev.off()

###
surface=c("Il7r","Il18r1","Ifngr1","Il18rap","Itga2")
FeaturePlot(immu,features = surface[1:4],split.by = "orig.ident")
FeaturePlot(immu,features = surface[5:8],split.by = "orig.ident")
p3=FeaturePlot(immu,features = surface,split.by = "orig.ident")
p3
ggsave("NK_surface.pdf",width = 7,height = 15)
dev.off()


#####subcelltype_enterocyte
fact1=c("TA","TA carcinoma","Enterocytes carcinoma","Plastical Enterocytes","Distal Enterocytes 1","Distal Enterocytes 2","Lig Enterocytes","seEPC","ASC cells")  
young11=subset(GastroSR1.downsampled.cells,ident=fact1)
FeaturePlot(young11,cols = c("lightgrey" ,"#DE1F1F"),features = c("Slc12a2","Nop56","Fermt1","Noxa1","Ncl"),pt.size = 0.01) & xlim(c(-10,7)) & ylim(c(-7,5))
FeaturePlot(young11,features = "Slc12a2",pt.size = 0.01) & xlim(c(-10,7)) & ylim(c(-7,5))
FeatureCornerAxes(young11,groupFacet = NULL,features  =  c("Slc12a2","Nop56","Fermt1","Noxa1","Ncl"),pSize =0.1) & xlim(c(-10,7)) & ylim(c(-7,5))
ggsave("subcelltype_enterocyte.pdf",width = 7,height = 15)
dev.off()






