#check if the data needed to be integrated and corrected for batch effect 
#make new object sam as seuart_hdf5 for checking batch effect
output <- seurat_hdf5

# perform standard workflow steps to figure out if we see any batch effects --------
output <- NormalizeData(object = output)
output <- FindVariableFeatures(object = output)
output <- ScaleData(object = output)
output <- RunPCA(object = output)
ElbowPlot(output)
output <- FindNeighbors(object = output, dims = 1:30)
output <- FindClusters(object = output)
output <- RunUMAP(object = output, dims = 1:30)

# plot
p1 <- DimPlot(output, reduction = 'umap', group.by = 'Batch')
p2 <- DimPlot(output, reduction = 'umap', group.by = 'Diagnosis',cols = c('red','green','blue'))
p3 <- DimPlot(output, reduction = 'umap', group.by = 'SampleID',cols = c('red','green','blue',"darkblue",'black','orange','darkgreen','purple','lavender','gray','magenta','pink','lightgreen',"darkred"))
p4 <-DimPlot(output, reduction = 'umap', group.by = 'seurat_clusters', label=TRUE)
p5 <- DimPlot(output, reduction = 'umap', group.by = 'Cell.Type',label=TRUE)
p6 <- DimPlot(output, reduction = 'umap', group.by = 'cluster', label=TRUE)

p3
p1|p2
p2|p4
p5|p6
 
###########################################################################################
# Step 03A: BATCH CORRECTION USING CCA FROM R
###########################################################################################

#if separted data then perform batch correction
# perform integration to correct for batch effects ------
obj.list <- SplitObject(output, split.by = 'Batch')
for(i in 1:length(obj.list)){
  obj.list[[i]] <- NormalizeData(object = obj.list[[i]])
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]])
}

# select integration features
features <- SelectIntegrationFeatures(object.list = obj.list)

# find integration anchors
anchors <- FindIntegrationAnchors(object.list = obj.list,
                                  anchor.features = features)

# integrate data
seurat.integrated <- IntegrateData(anchorset = anchors)


# Scale data, run PCA and UMAP and visualize integrated data
seurat.integrated <- ScaleData(object = seurat.integrated)
seurat.integrated <- RunPCA(object = seurat.integrated)
seurat.integrated <- RunUMAP(object = seurat.integrated, dims = 1:50)

#visualize the changes after batch correction
p.1.1 <- DimPlot(output, reduction = 'umap', group.by = 'Batch')
p.1.2 <- DimPlot(output, reduction = 'umap', group.by = 'Diagnosis',cols = c('red','green','blue'))
p.1.3 <- DimPlot(output, reduction = 'umap', group.by = 'SampleID',cols = c('red','green','blue',"darkblue",'black','orange','darkgreen','purple','lavender','gray','magenta','pink','lightgreen',"darkred"))
p.1.4 <-DimPlot(output, reduction = 'umap', group.by = 'seurat_clusters', label=TRUE)
p.1.5 <- DimPlot(output, reduction = 'umap', group.by = 'Cell.Type',label=TRUE)
p.1.6 <- DimPlot(output, reduction = 'umap', group.by = 'cluster', label=TRUE)


rm(p3,seurat.integrated, anchors,features,obj.list)

#########################################################################################
# Step 03B: BATCH CORRECTION USING iNMF METHOD
#########################################################################################

#to perform integration between batches using LIGER, which is based on iNMF (integrative Non-negative Matrix Factorization)
#liger workflow 

#splitting the dataset based on Diagnosis
seuart_hdf5[["RNA"]] <- split(seuart_hdf5[["RNA"]], f = seuart_hdf5$Diagnosis)

#to work with seurat object we don't need to create liger object
#filtering 
seuart_hdf5 <- subset(seuart_hdf5, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 5)

#standard workflow

seuart_hdf5 <- seuart_hdf5 %>%
  normalize() %>%
  selectGenes() %>%
  scaleNotCenter()

view(seuart_hdf5)

#performing iNMF integration on the two dataset

seuart_hdf5 <- seuart_hdf5 %>%
  runINMF(k = 20) %>%
  quantileNorm()
view(seuart_hdf5)

#Visualizing the integration result when cell annotation is already provided

seuart_hdf5 <- RunUMAP(seuart_hdf5, reduction = "inmfNorm", dims = 1:20)
gg.byDataset <- DimPlot(seuart_hdf5, group.by = "Diagnosis")
gg.byCelltype <- DimPlot(seuart_hdf5, group.by = "Cell.Type")
gg.byDataset + gg.byCelltype

DimPlot(ifnb, group.by = "Cell.Type", split.by = "Diagnosis")

#in case when annotation is performed after integration

seuart_hdf5 <- seuart_hdf5 %>%
  FindNeighbors(reduction = "inmfNorm", dims = 1:20) %>%
  FindClusters()

DimPlot(seuart_hdf5, group.by = "seurat_clusters")