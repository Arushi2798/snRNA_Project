#####
#check if the data needed to be integrated and corrected for batch effect 

library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)

# perform standard workflow steps to figure out if we see any batch effects --------
seurat_hdf5 <- NormalizeData(object = seurat_hdf5)
seurat_hdf5 <- FindVariableFeatures(object = seurat_hdf5)
seurat_hdf5 <- ScaleData(object = seurat_hdf5)
seurat_hdf5 <- RunPCA(object = seurat_hdf5)
ElbowPlot(seurat_hdf5)
seurat_hdf5 <- FindNeighbors(object = seurat_hdf5, dims = 1:30)
seurat_hdf5 <- FindClusters(object = seurat_hdf5)
seurat_hdf5 <- RunUMAP(object = seurat_hdf5, dims = 1:30)

# plot
p1 <- DimPlot(seurat_hdf5, reduction = 'umap', group.by = 'Batch')
p2 <- DimPlot(seurat_hdf5, reduction = 'umap', group.by = 'Diagnosis',cols = c('red','green','blue'))
p3 <- DimPlot(seurat_hdf5, reduction = 'umap', group.by = 'SampleID',cols = c('red','green','blue',"darkblue",'black','orange','darkgreen','purple','lavender','gray','magenta','pink','lightgreen',"darkred"))
p4 <-DimPlot(seurat_hdf5, reduction = 'umap', group.by = 'seurat_clusters', label=TRUE)
p4 <-DimPlot(seurat_hdf5, reduction = 'umap', group.by = 'seurat_clusters') + NoLegend()
p5 <- DimPlot(seurat_hdf5, reduction = 'umap', group.by = 'Cell.Type',label=TRUE)
p6 <- DimPlot(seurat_hdf5, reduction = 'umap', group.by = 'cluster', label=TRUE)

p1|p2
p2|p4
p5|p6

#grid.arrange(p1, p2, ncol = 2, nrow = 2)

#####
#if separted data then perform batch correction using CCA
# perform integration to correct for batch effects ------
obj.list <- SplitObject(seurat_hdf5, split.by = 'Batch')
for(i in 1:length(obj.list)){
  obj.list[[i]] <- NormalizeData(object = obj.list[[i]])
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]])
}


# select integration features
features <- SelectIntegrationFeatures(object.list = obj.list)

# find integration anchors (CCA)
anchors <- FindIntegrationAnchors(object.list = obj.list,
                                  anchor.features = features)

# integrate data
seurat.integrated <- IntegrateData(anchorset = anchors)


# Scale data, run PCA and UMAP and visualize integrated data
seurat.integrated <- ScaleData(object = seurat.integrated)
seurat.integrated <- RunPCA(object = seurat.integrated)
seurat.integrated <- RunUMAP(object = seurat.integrated, dims = 1:50)


#p3 <- DimPlot(seurat.integrated, reduction = 'umap', group.by = 'Batch')
#p4 <- DimPlot(seurat.integrated, reduction = 'umap', group.by = 'Type',
 #             cols = c('red','green','blue'))
#p3
#grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2)
#grid.arrange(p1, p3, ncol = 2, nrow = 2)

rm(p3,seurat.integrated, anchors,features,obj.list)
