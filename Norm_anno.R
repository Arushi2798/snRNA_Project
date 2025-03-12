################################################################################
# Step 04: STANDARD WORKFLOW
################################################################################

#1_1. either perform SCTranform only as  this single command replaces NormalizeData(), ScaleData(), and FindVariableFeatures()
seurat_hdf5 = SCTransform(seurat_hdf5, verbose = TRUE,variable.features.n = 2000,vst.flavor = "v2s")

    #view the detection rate of genes
genes = seurat_hdf5@assays$RNA@meta.data[rownames(seurat_hdf5), ]
genes = data.frame(genes, seurat_hdf5@assays[["SCT"]]@SCTModel.list[["counts"]]@feature.attributes, row.names = rownames(genes), stringsAsFactors = FALSE)
hist(genes[, 'detection_rate'])

#1_2. or perform normalizedata instead of sctransform
seurat_hdf5 <- NormalizeData(seurat_hdf5)

# 2. Identify highly variable features
seurat_hdf5 <- FindVariableFeatures(seurat_hdf5, nfeatures = 2000)

#save the variable features for futhure downstream analysis like finding hub genes via GRN or WGCNA
hvgc<-VariableFeatures(seurat_hdf5)

    # Identify the 10 most highly variable genes
top20 <- head(VariableFeatures(seurat_hdf5), 20)

    # plot variable features with and without labels
plot1 <- VariableFeaturePlot(seurat_hdf5)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)

plot1|plot2

#5. scaling
seurat_hdf5 <- ScaleData(seurat_hdf5)

#6. PCA
seurat_hdf5 = RunPCA(seurat_hdf5, verbose = FALSE)

    # determine dimensionality of the data
ElbowPlot(seurat_hdf5,ndims = 50)

    # visualize PCA results
print(seurat_hdf5[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(seurat_hdf5, dims = 1:2, reduction = "pca")

DimPlot(seurat_hdf5, reduction= "pca") + NoLegend()

DimHeatmap(seurat_hdf5, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(seurat_hdf5, dims = 1:20, cells = 500, balanced = TRUE)

    #PCA Visualization

total_variance <- sum(matrixStats::rowVars(Seurat::GetAssayData(seurat_hdf5, assay = "SCT", layer = "scale.data")))
eigValues =  (seurat_hdf5$pca@stdev)^2

    #Plot PCA elbow
pca.pdata = data.frame(PC = 1:length(seurat_hdf5$pca@stdev),
                       eigValues = (seurat_hdf5$pca@stdev)^2,  ## EigenValues
                       varExplained = eigValues / total_variance * 100,
                       totalExplained = cumsum(eigValues / total_variance * 100))

pobj = ggplot(data = pca.pdata,aes(x = PC,y = varExplained)) + geom_point() + geom_line() + 
  geom_vline(xintercept = c(5,7,20),colour = c("red","blue","green")) + 
  theme_bw() + 
  labs(x = "PCs",y = "%. Variance Explained") + 
  theme(axis.title = element_text(size = 20),axis.text = element_text(size = 17))

print(pobj)
dev.off()

#7_1. non-linear dimensionality reduction by UMAP
seurat_hdf5 <- RunUMAP(seurat_hdf5, reduction = "inmf", dims = 1:dim(seurat_hdf5[["inmf"]])[2])
    # note that you can set `label = TRUE` or use the LabelClusters function to help label individual clusters

#7_2. non-linear dimensionality reduction by tSNE
seurat_hdf5 <- RunTSNE(seurat_hdf5, reduction = "inmf", dims = 1:dim(seurat_hdf5[["inmf"]])[2])

#visualize UMAP plot
umap1<-DimPlot(seurat_hdf5, reduction = "umap") + NoLegend()

#visualize TSNE plot
tsne1<-DimPlot(seurat_hdf5, reduction = "tsne", label = TRUE) 

#8. Clustering 
seurat_hdf5 <- FindNeighbors(seurat_hdf5, reduction = "inmf", dims = 1:dim(seurat_hdf5[["inmf"]])[2], nn.eps=0.5)

    # understanding resolution
seurat_hdf5 <- FindClusters(seurat_hdf5, resolution = c(0.1,0.2,0.3, 0.5,0.9, n.start=10))


    #used for cluster separation after reductions 
seurat_hdf5 <- FindClusters(seurat_hdf5, resolution =0.1) #(final selection)

    #visualization
cl1<-DimPlot(seurat_hdf5, group.by = "RNA_snn_res.0.1", label = TRUE)
cl2<-DimPlot(seurat_hdf5, group.by = "RNA_snn_res.0.2", label = TRUE)
cl3<-DimPlot(seurat_hdf5, group.by = "RNA_snn_res.0.3", label = TRUE)
cl4<-DimPlot(seurat_hdf5, group.by = "RNA_snn_res.0.5", label = TRUE)
cl5<-DimPlot(seurat_hdf5, group.by = "RNA_snn_res.0.9", label = TRUE)

#check the identy of the clusteers, levels of clusters formed
head(Idents(seurat_hdf5),5)

    #set the identity of the dataset for further ananlysis on the basis of cluster resolution
Idents(seurat_hdf5) <- "RNA_snn_res.0.1" #if resolution of 0.1 preferred

#9. findConserved markers
#since data layers are not joined due to batch correction.Running JoinLayers
seurat_hdf5<-JoinLayers(seurat_hdf5)
DefaultAssay(seurat_hdf5) <- "RNA"

# Create a list to store the markers for each cluster
all_cluster_markers <- list()

# Get unique cluster identifiers
clusters <- levels(Idents(seurat_hdf5)) # Assuming `Idents` contains the cluster information

# Loop through each cluster and find conserved markers
for (cluster in clusters) {
  # Find conserved markers for the current cluster
  markers <- FindConservedMarkers(
    seurat_hdf5, 
    ident.1 = cluster, 
    grouping.var = 'Diagnosis', 
    verbose = TRUE  # Optional, for detailed seurat_hdf5 during processing
  )
  
  # Store the markers in the list, using the cluster number as the key
  all_cluster_markers[[paste0("Cluster_", cluster)]] <- markers
  # Store markers for this cluster as metadata
  #seurat_hdf5@meta.data[[paste0("markers_cluster_", cluster)]] <- markers
  
}

rm(markers)

# Inspect the results for a specific cluster (e.g., Cluster 7)
head(all_cluster_markers[["Cluster_7"]])

#saveRDS(marker_results, file = "marker_results.rds")

#visualize top Features

fp1 <- FeaturePlot(seurat_hdf5,features = c("gene_name"), min.cutoff = "q10")
umap2 <- DimPlot(seurat_hdf5, reduction = 'umap', group.by = 'Cell.Type')