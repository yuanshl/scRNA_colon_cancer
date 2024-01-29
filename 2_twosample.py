#!/bin/python
import scvelo as scv 
import scanpy as sc 
import cellrank as cr 
import numpy as np 
import pandas as pd 
import anndata as ad 
scv.settings.verbosity = 3 
scv.settings.set_figure_params('scvelo', facecolor='white', dpi=100, frameon=False) 
cr.settings.verbosity = 2 
adata = sc.read_h5ad('my_data.h5ad') 
# load loom files for spliced/unspliced matrices for each sample: 
ldata1 = scv.read('W_35.loom', cache=True) 
# rename barcodes in order to merge: 
barcodes = [bc.split(':')[1] for bc in ldata1.obs.index.tolist()] 
# 查看seurat的barcode命名方式，也可以是'Sample1_'
#barcodes = [bc[0:len(bc)-1] + '_10' for bc in barcodes] 
barcodes = ['WT_' + bc[0:len(bc)-1] + '-1'  for bc in barcodes] ##change it for EFE002, 'EFE_' not 'EFE001_'
ldata1.obs.index = barcodes 
# make variable names unique 
ldata1.var_names_make_unique() 

ldata2 = scv.read('N_35.loom', cache=True)
# rename barcodes in order to merge:
barcodes = [bc.split(':')[1] for bc in ldata2.obs.index.tolist()]
# 查看seurat的barcode命名方式，也可以是'Sample1_'
#barcodes = [bc[0:len(bc)-1] + '_10' for bc in barcodes]
barcodes = ['Ncf4_KO_' + bc[0:len(bc)-1] + '-1'  for bc in barcodes] ##change it for EFE002, 'EFE_' not 'EFE001_'
ldata2.obs.index = barcodes
# make variable names unique
ldata2.var_names_make_unique()

# concatenate the three loom 
ldata = ldata1.concatenate([ldata2],index_unique='@@') 
barcodesM = [bc.split('@@')[0] for bc in ldata.obs.index.tolist()]
ldata.obs.index = barcodesM

# merge matrices into the original adata object 
adata = scv.utils.merge(adata, ldata) 
# plot umap to check 
sc.pl.umap(adata, color='celltype', frameon=False, legend_loc='on data', title='', save='_celltypes.pdf')

fact=["EC","EndoC","Fibroblast","SMC","Cardiomyocyte","Macrophage","T-cell","Pericyte","Neuro"]
fact=["C0","C1","C2","C3","C4","C5","C6"]
fact=["TA","TA carcinoma","Enterocytes carcinoma","Goblet","Plastical Enterocytes","Distal Enterocytes 1","Distal Enterocytes 2","Lig Enterocytes","seEPC","ASC cells"]
a=adata.obs.celltype
celltypeOder=a.cat.set_categories(fact)
adata.obs['celltypeOder']=celltypeOder

color=("#00BFFF","#FFA500","#DDA0DD","#FF69B4","#6A5ACD","#9932CC","#32CD32","#DB7093","#808080","#C71585","#87CEFA","#00CED1","#8B0000","#BC8F8F","#F08080","#DAA520","#90EE90","#3CB371","#DC143C","#191970")


scv.pp.filter_and_normalize(adata)
scv.pp.moments(adata)

scv.tl.velocity(adata, mode='stochastic')
scv.tl.velocity_graph(adata)
scv.pl.velocity_embedding(adata, basis='umap', frameon=False, save='embedding.pdf')
scv.pl.velocity_embedding_grid(adata, basis='umap', color=['celltypeOder'], save='embedding_grid.pdf', title='', scale=0.25,size=10)
scv.pl.velocity_embedding_stream(adata, basis='umap', color=['celltypeOder'], save='embedding_stream.pdf', title='EFE001_EFE002',palette=color,size=10,legend_fontsize=8,legend_loc ='on data')


