################################################################################
# Step 07: FIND DIFFERENTIALLY EXPRESSED GENES ACCROSS CONDITION
################################################################################

theme_set(theme_cowplot())
aggregate_info <- AggregateExpression(seu1,group.by = c("cluster","Diagnosis"),return.seurat =TRUE)
genes.to.label=c("RNF219-AS1","","","","","","","","",'')
cs1 <- CellScatter(aggregate_info,'ASC1_AD','ASC1_Control', highlight = genes.to.label)
cs2 <- CellScatter(aggregate_info,'ASC2_AD','ASC2_Control', highlight = genes.to.label)
cs3 <- CellScatter(aggregate_info,'ASC3_AD','ASC3_Control', highlight = genes.to.label)
cs4 <- CellScatter(aggregate_info,'ASC4_AD','ASC4_Control', highlight = genes.to.label)

cs1|cs2|cs3|cs4

#########################################################################################
# Step 08_1: FIND DIFFERENTIAL EXPRESSED GENES BETWEEN CONDITIONS (for each cell types)
#########################################################################################

# find markers for astrocytes
d1<- FindMarkers(seurat_hdf5, ident.1 = 'ASC_AD', ident.2 = 'ASC_Control')
# find markers for inhibitory neurons
i1 <-FindMarkers(seurat_hdf5, ident.1 = 'INH_AD', ident.2 = 'INH_Control')
# find markers for excitatory neurons
d2<- FindMarkers(seurat_hdf5, ident.1 = 'EX_AD', ident.2 = 'EX_Control')
# find markers for microglia
d4 <- FindMarkers(seurat_hdf5, ident.1 = 'MG_AD', ident.2 = 'MG_Control')
# find markers for oligo.Progenators
d6 <- FindMarkers(seurat_hdf5, ident.1 = 'OPC_AD', ident.2 = 'OPC_Control')
# find markers for oligodendrocytes
d1 <-FindMarkers(seurat_hdf5, ident.1 = 'ODC_AD', ident.2 = 'ODC_Control')
# find markers for Pericytes Endothelial
d7 <- FindMarkers(seurat_hdf5, ident.1 = 'PER.END_AD', ident.2 = 'PER.END_Control')

head(d1)


########################################################################################################
# Step 08_2: FIND DIFFERENTIAL EXPRESSED GENES BETWEEN CONDITIONS (for each sub cell types)
########################################################################################################

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