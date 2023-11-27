install.packages("Seurat")
install.packages("remotes")
remotes::install_github("lazappi/clustree@develop")#install package from github
BiocManager::install("GSVA")
install.packages("msigdbr")
BiocManager::install("GSEABase")
install.packages("pheatmap")
install.packages("harmony")

library(Seurat)
library(stringr)
library(clustree)
library(dplyr)
library(GSVA)
library(msigdbr)
library(GSEABase)
library(pheatmap)
library(harmony)
rm(list = ls())
gc()#free up memory
getwd()#show the path
setwd("/lustre/work/client/users/yuannaw/singlecells")
#1.read data from each sample

filename <- paste('/lustre/work/client/users/yuannaw/singlecells/',list.files('/lustre/work/client/users/yuannaw/singlecells/'),sep = '')
filename
class(filename)
sceList <- lapply(filename, function(x){obj <- CreateSeuratObject(counts = Read10X(x),
                                                                  project = str_split(x,'/')[[1]][8])})

data_directory <- "/lustre/work/client/users/yuannaw/singlecells/T59"
counts = Read10X(data_directory)
T59_counts<-Read10X(data.dir = "/lustre/work/client/users/yuannaw/singlecells/T59")
T59<-CreateSeuratObject(counts = T59_counts,
                        min.features = 100)
head(T59@meta.data)

for (fileobj in c("T59", "T76","T77","T89","T90","X2","X3","X4")){
  seurat_data <- Read10X(data.dir = paste0("/lustre/work/client/users/yuannaw/singlecells/", fileobj))
  seurat_obj <- CreateSeuratObject(counts = seurat_data, 
                                   project = fileobj)
  assign(fileobj, seurat_obj)
}
seurat_obj@project.name
head(T59@meta.data)
head(T76@meta.data)
class(seurat_obj)
# Create a merged Seurat object
seurat_list <- list(T59,T76,T77,T89,T90,X2,X3,X4)
sce<- Reduce(function(x, y) merge(x, y), seurat_list)
dim(sce)
# 38224 genes, 67323 cells
# check the percentage of mitochondrial genes(MT)and the ribosomal genes(RPS/RPL)
grep('^MT',x=rownames(sce@assays$RNA@data),value = T)
grep('^RP[SL]',x=rownames(sce@assays$RNA@data),value = T)


sce <- PercentageFeatureSet(sce,pattern = '^MT',col.name = 'percent.MT')
sce <- PercentageFeatureSet(sce,pattern = '^RP[SL]',col.name = 'percent.RP')
VlnPlot(sce,features = 'percent.MT',pt.size = 0)
VlnPlot(sce,features = 'percent.RP',pt.size = 0)
VlnPlot(sce,features = 'nCount_RNA',pt.size = 0)
VlnPlot(sce,features = 'nFeature_RNA',pt.size = 0)
#filter cells
sce <- subset(sce,subset = nCount_RNA>3000 & nFeature_RNA>300 & percent.MT<40)
dim(sce) # 38224 genes， 29389 cells

#filter genes
sce <- sce[rowSums(sce@assays$RNA@counts>0)>3,]
dim(sce) # 28074 genes 29389 cells

#check batch effect
sce <- NormalizeData(sce) # logNormalize
sce <- FindVariableFeatures(sce) # Top2000 highly variable genes
sce <- ScaleData(sce) # use variableFeature to scale
sce <- RunPCA(sce) # use variableFeature to dedimension
DimPlot(sce,reduction = 'pca',group.by = 'orig.ident')

#use SCTransform to normalize
sce <- SCTransform(sce) # this step includes NormalizeData, ScaleData and FindVariableFeatures
DefaultAssay(sce)
sce <- RunPCA(sce)
DimPlot(sce,reduction = 'pca',group.by = 'orig.ident')
ElbowPlot(sce,ndims = 40)
#7. use umap to dedimension, then subgroup, based on PCA
sce <- RunUMAP(sce,reduction = 'pca',dims = 1:40)
sce <- FindNeighbors(sce,dims = 1:40)

sce_res <- sce
for (i in c(0.01, 0.05, 0.1, 0.15, 0.2, 0.3,0.4, 0.5,0.8,1)){
  sce_res <- FindClusters(sce_res,resolution = i)
}
clustree(sce_res,prefix = 'SCT_snn_res.')
sce <- FindClusters(sce,resolution = 0.4)
DimPlot(sce,reduction = 'umap',group.by = 'seurat_clusters',label = T)

#manual cells annotation

# immune cells  5,6,8,9,15,20 
DotPlot(sce,features = c('PTPRC','CD45'))
# T cells   0,[8],9,13
DotPlot(sce,features = c('CD3D','CD3E','CD8A')) 
# B cells   14,16,18,[20]
DotPlot(sce,features = c('CD79A', 'CD37', 'CD19', 'CD79B', 'MS4A1','CD20'))
# plasma cells [14,16,18]
DotPlot(sce,features = c('IGHG1','MZB1','SDC1','CD79A'))
# NK cells  5,[9,15]
DotPlot(sce,features = c('FGFBP2','FCG3RA','CX3CR1'))
# NK cells  8,15
DotPlot(sce,features = c('CD160','NKG7','GNLY','CD247','CCL3','GZMB','CXCR1','TYOB','PRF1'))
# Endothelial cells [11]
DotPlot(sce,features = c('PECAM1','VWF'))
# fibroblast [0,7],10,11
DotPlot(sce,features = c('FGF7','MME','DCN','LUM','GSN','PF4','PPBP'))
# epithelial cell 1,2,[3],4,[12,13,17,19]
DotPlot(sce,features = c('EPCAM','KRT19','PROM1','ALDH1A1','CD24'))
# myeloid cells 5,6,9,15
DotPlot(sce,features = c('CD68','CD163','CD14','LYZ'))
# Ovarian somatic cell  [13]
DotPlot(sce,features = c('LGR5'))

marker <- data.frame(cluster = 0:20,cell = 0:20)
marker[marker$cluster %in% c(5,6,8,9,15,20),2] <- 'Immune cells'
marker[marker$cluster %in% c(0,8,9,13),2] <- 'T cells'
marker[marker$cluster %in% c(14,16,18,20),2] <- 'B cells'
marker[marker$cluster %in% c(14,16,18),2] <- 'plasma cells'
marker[marker$cluster %in% c(5,9,15),2] <- 'NK cells'
marker[marker$cluster %in% c(11),2] <- 'Endothelial cells'
marker[marker$cluster %in% c(5,6),2] <- 'myeloid cells'
marker[marker$cluster %in% c(1,2,3,4,12,13,17,19),2] <- 'epithelial cell'
marker[marker$cluster %in% c(0,7,10),2] <- 'fibroblast'
marker[marker$cluster %in% c(13),2] <- 'Ovarian somatic cell'

sce@meta.data$cell_type <- sapply(sce@meta.data$seurat_clusters,function(x){marker[x,2]})
DimPlot(sce,reduction = 'umap',group.by = 'cell_type',label = T)

#heatmap of marker gene 
heatmap_marker_gene <- c('CD3D','CD3E','CD8A',
                         'CD79A', 'CD19', 'CD79B', 'MS4A1','CD20',
                         'MZB1','SDC1','CD79A',
                         'CD68','CD163','CD14','LYZ')
heatmap_cells <- rownames(sce@meta.data[sce@meta.data$cell_type %in% c("myeloid cells","B cells","T cells"),])
sub_sce <- subset(sce,cells = heatmap_cells)
DoHeatmap(sub_sce,features = heatmap_marker_gene,group.by = 'cell_type',size = 3)
