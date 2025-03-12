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

gc()

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
gc()

rm(markers)

# Inspect the results for a specific cluster (e.g., Cluster 7)
head(all_cluster_markers[["Cluster_7"]])
#saveRDS(marker_results, file = "marker_results.rds")

#visualize top Features

fp1 <- FeaturePlot(seurat_hdf5,features = c("RNF219-AS1"), min.cutoff = "q10")
fp+NoLegend()|fp1

umap2 <- DimPlot(seurat_hdf5, reduction = 'umap', group.by = 'Cell.Type')
fp+NoLegend()|fp1|p5
################################################################################
# Step 07: find differential expressed genes across conditions
################################################################################

theme_set(theme_cowplot())
aggregate_info <- AggregateExpression(seu1,group.by = c("cluster","Diagnosis"),return.seurat =TRUE)
genes.to.label=c(""RNF219-AS1"","","","","","","","","",'')
cs1 <- CellScatter(aggregate_info,'ASC1_AD','ASC1_Control', highlight = genes.to.label)
cs2 <- CellScatter(aggregate_info,'ASC2_AD','ASC2_Control', highlight = genes.to.label)
cs3 <- CellScatter(aggregate_info,'ASC3_AD','ASC3_Control', highlight = genes.to.label)
cs4 <- CellScatter(aggregate_info,'ASC4_AD','ASC4_Control', highlight = genes.to.label)

cs1|cs2|cs3|cs4

################################################################################
# Step 08: find differential expressed genes between conditions
################################################################################

Idents(seurat_hdf5) <- seurat_hdf5@meta.data$cluster

seurat_hdf5$cell.type.cnd <- paste0(seurat_hdf5$cluster,'_', seurat_hdf5$Diagnosis)
View(seurat_hdf5@meta.data)

Idents(seurat_hdf5) <- seurat_hdf5$cell.type.cnd

compare<-DimPlot(seurat_hdf5, reduction = 'umap', label = TRUE) +NoLegend()


# Define all cell types
cell_types <- c("INH1", "ODC2", "MG3", "ODC6", "ODC1", "ODC7", "EX2",
                "INH2", "ODC11", "ODC9", "MG1", "OPC1", "EX1", "INH3",
                "ASC1", "ODC8", "EX3", "INH4", "PER.END1", "EX4", "ASC3",
                "ASC2", "ODC4", "OPC2", "ODC5", "ODC13", "ODC3", "MG2",
                "ODC12", "EX5", "ODC10", "PER.END3", "PER.END2", "ASC4")

# Initialize an empty list to store results
marker_results <- list()

# Loop through each cell type and compute differentially expressed genes
for (cell in cell_types) {
  ident_ad <- paste0(cell, "_AD")
  ident_control <- paste0(cell, "_Control")
  
  # Find markers
  marker_results[[cell]] <- FindMarkers(seurat_hdf5, ident.1 = ident_ad, ident.2 = ident_control)
}

# Check the results
marker_results[["ASC1"]]  # Example: Accessing results for ASC1



# find markers for astrocytes
d1<- FindMarkers(seurat_hdf5, ident.1 = 'ASC1_AD', ident.2 = 'ASC1_Control')
d2<- FindMarkers(seurat_hdf5, ident.1 = 'ASC2_AD', ident.2 = 'ASC2_Control')
d3<- FindMarkers(seurat_hdf5, ident.1 = 'ASC3_AD', ident.2 = 'ASC3_Control')
d4<- FindMarkers(seurat_hdf5, ident.1 = 'ASC4_AD', ident.2 = 'ASC4_Control')

# find markers for inhibitory neurons
i1 <-FindMarkers(seurat_hdf5, ident.1 = 'INH1_AD', ident.2 = 'INH1_Control')
i2 <-FindMarkers(seurat_hdf5, ident.1 = 'INH2_AD', ident.2 = 'INH2_Control')
i3 <-FindMarkers(seurat_hdf5, ident.1 = 'INH3_AD', ident.2 = 'INH3_Control')
i4 <-FindMarkers(seurat_hdf5, ident.1 = 'INH4_AD', ident.2 = 'INH4_Control')

# find markers for excitatory neurons
# find markers for microglia
# find markers for oligo.Progenators
# find markers for oligodendrocytes
O1 <-FindMarkers(seurat_hdf5, ident.1 = 'ODC_AD', ident.2 = 'INH1_Control')

# find markers for Pericytes Endothelial
d2<- FindMarkers(seurat_hdf5, ident.1 = 'EX_AD', ident.2 = 'EX_Control')
d3 <- FindMarkers(seurat_hdf5, ident.1 = 'INH_AD', ident.2 = 'INH_Control')
d4 <- FindMarkers(seurat_hdf5, ident.1 = 'MG_AD', ident.2 = 'MG_Control')
d5 <- FindMarkers(seurat_hdf5, ident.1 = 'ODC_AD', ident.2 = 'ODC_Control')
d6 <- FindMarkers(seurat_hdf5, ident.1 = 'OPC_AD', ident.2 = 'OPC_Control')
d7 <- FindMarkers(seurat_hdf5, ident.1 = 'PER.END_AD', ident.2 = 'PER.END_Control')

head(d1)

# plotting conserved features vs DE features between conditions
head(all_cluster_markers[["Cluster_"]])

Z2<-FeaturePlot(seurat_hdf5, features = c("GINS3","MALAT1","CNN3","CYP4F3")
                                     , split.by = 'Diagnosis', min.cutoff = "q10")

Z1<-FeaturePlot(seurat_hdf5, features = c("AC063952.1","TBC1D3D","SIM2","GRTP1"), split.by = 'Diagnosis', min.cutoff = "q10")
gc()
vn3<-VlnPlot(seurat_hdf5, features= c("AC063952.1","TBC1D3D","SIM2","GRTP1"),split.by="Diagnosis", group.by= "cluster", pt.size=0, combine=FALSE)

vn4<- wrap_plots(plots=vn3, ncol=1) 

################################################################################
# Step 09: find all markers between clusters
################################################################################

seurat.markers <-FindAllMarkers(seurat_hdf5, only.pos = TRUE)

gc()
saveRDS(seurat.markers, file = "seurat.markers.rds")

seurat.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

cluster7.markers <- findMarkers(seurat_hdf5,ident.1=7,logfc.threshold=0.25,test.use ="VCAN")

##visualizations
vn1<-VlnPlot(seurat_hdf5, features= c('LINC01608', 'CD22'))
vn2<-VlnPlot(seurat_hdf5, features= c('MOG', 'LDB3'), layer= "counts", log=TRUE)
#some error in this next plot fp2
fp2<-FeaturePlot(seurat_hdf5, features=c('LINC01608', 'CD22','MOG', 'LDB3','OPALIN', 'MOBP'))

seurat.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC <1) %>%
  slice_head(n=10) %>%
  ungroup()-> top10

hp1<- DoHeatmap(seurat_hdf5, features =top10$gene) + NoLegend()