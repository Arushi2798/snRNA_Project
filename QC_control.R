library(ggplot2)
library(dplyr)
library(patchwork)
library(ggrepel)
library(scran)
library(glmGamPoi)
library(future)

seurat_hdf5@meta.data$orig.ident = NULL
seurat_hdf5@assays$RNA@meta.data = data.frame(Geneid = rownames(seurat_hdf5), Symbol = rownames(seurat_hdf5), row.names = rownames(seurat_hdf5), stringsAsFactors = FALSE)


View(seurat_hdf5@meta.data)

# Quality Control
#calculate mitochondrial percentage
seurat_hdf5[["percent.mt"]] <- PercentageFeatureSet(seurat_hdf5, pattern = "^MT-")
View(seurat_hdf5@meta.data)

summary(seurat_hdf5[[]]$percent.mt)

# to vizualise the distribution
# Distribution of nFeature_RNA[The number of detected genes per cell],nCount_RNA[The number of unique molecular identifiers (UMIs) per cell,
# representing RNA abundance],percent.mt (Mitochondrial Content)

VlnPlot(seurat_hdf5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <-FeatureScatter(seurat_hdf5, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')
plot2 <-FeatureScatter(seurat_hdf5, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot1|plot2

 
meta1 = seurat_hdf5@meta.data
write.table(meta1, file=paste0('meta.seurat.raw.tsv'), quote = FALSE, sep = '\t', row.names = FALSE)
meta1 = meta1[rev(order(meta1$Diagnosis)), ]
meta1$sid = unsplit(lapply(split(meta1, meta1$Diagnosis), function(x) paste0(c(AD = 'AD', Control = 'Ct')[x$Diagnosis[1]], as.integer(factor(x$SampleID, level = unique(x$SampleID))))), meta1$Diagnosis)
meta1$sid = factor(meta1$sid, levels = unique(meta1$sid))
theme_bw2 = function() theme_bw() + theme(panel.grid = element_blank())
p1 = ggplot(meta1, aes(x = sid, y = nFeature_RNA, fill = Diagnosis)) + geom_violin(width = 0.8, scale = 'width') + theme_bw2() + theme(legend.position = 'none', axis.text.x = element_text(angle = 45, size = 6, hjust = 1), axis.text.y = element_text(size = 7)) + xlab(NULL)
p2 = ggplot(meta1, aes(x = sid, y = nCount_RNA, fill = Diagnosis)) + geom_violin(width = 0.8, scale = 'width') + theme_bw2() + theme(legend.position = 'none', axis.text.x = element_text(angle = 45, size = 6, hjust = 1), axis.text.y = element_text(size = 7)) + xlab(NULL)
print(p1 / p2)

seurat_hdf5@meta.data$SampleID = factor(seurat_hdf5@meta.data$SampleID)
p3 = FeatureScatter(seurat_hdf5, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size=0.6, group.by='SampleID') + theme(title=element_text(size=10), legend.text=element_text(size=10), legend.title=element_text(size=12))
print(p3)

 
# Filter the Seurat object based on the suggested thresholds
seurat_hdf5 <- subset(seurat_hdf5, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 5)


#either perform SCTranform only as  this single command replaces NormalizeData(), ScaleData(), and FindVariableFeatures()
seurat_hdf5 = SCTransform(seurat_hdf5, verbose = TRUE,variable.features.n = 2000,vst.flavor = "v2s")

#view the detection rate of genes
genes = seurat_hdf5@assays$RNA@meta.data[rownames(seurat_hdf5), ]
genes = data.frame(genes, seurat_hdf5@assays[["SCT"]]@SCTModel.list[["counts"]]@feature.attributes, row.names = rownames(genes), stringsAsFactors = FALSE)
hist(genes[, 'detection_rate'])


#or perform normalizedata instead of sctransform
seurat_hdf5 <- NormalizeData(seurat_hdf5)

# 4. Identify highly variable features 
seurat_hdf5 <- FindVariableFeatures(seurat_hdf5, nfeatures = 2000)

# Identify the 10 most highly variable genes
top20 <- head(VariableFeatures(seurat_hdf5), 20)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(seurat_hdf5)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
plot1
plot2
plot1|plot2


#scaling
seurat_hdf5 <- ScaleData(seurat_hdf5)
 
# PCA
seurat_hdf5 = RunPCA(seurat_hdf5, verbose = FALSE)
# visualize PCA results
print(seurat_hdf5[["pca"]], dims = 1:5, nfeatures = 5)

# determine dimensionality of the data
ElbowPlot(seurat_hdf5,ndims = 50)


VizDimLoadings(seurat_hdf5, dims = 1:2, reduction = "pca")

DimPlot(seurat_hdf5, reduction= "pca") + NoLegend()

DimHeatmap(seurat_hdf5, dims = 1, cells = 500, balanced = TRUE)


DimHeatmap(seurat_hdf5, dims = 1:20, cells = 500, balanced = TRUE)
 
#PCA Visualization

total_variance <- sum(matrixStats::rowVars(Seurat::GetAssayData(seurat_hdf5, assay = "SCT", layer = "scale.data")))
eigValues =  (seurat_hdf5$pca@stdev)^2

####
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


 
# non-linear dimensionality reduction 
seurat_hdf5 <- RunUMAP(seurat_hdf5, dims = 1:30, verbose = FALSE,seed.use = 4867)
# note that you can set `label = TRUE` or use the LabelClusters function to help label individual clusters
seurat_hdf5 = RunTSNE(seurat_hdf5, dims = 1:20, verbose = FALSE)

DimPlot(seurat_hdf5, reduction = "umap") + NoLegend()


str(seurat_hdf5)


 
# Clustering 
seurat_hdf5 <- FindNeighbors(seurat_hdf5, dims = 1:30)

# understanding resolution
seurat_hdf5 <- FindClusters(seurat_hdf5, resolution = c(0.3, 0.5, 0.7,1))
View(seurat_hdf5@meta.data)

#used 
seurat_hdf5 <- FindClusters(seurat_hdf5, resolution =0.3)

#visualization
DimPlot(seurat_hdf5, group.by = "SCT_snn_res.0.3", label = TRUE)
DimPlot(seurat_hdf5, group.by = "SCT_snn_res.0.5", label = TRUE)
DimPlot(seurat_hdf5, group.by = "SCT_snn_res.0.7", label = TRUE)
DimPlot(seurat_hdf5, group.by = "SCT_snn_res.1", label = TRUE)

head(Idents(seurat_hdf5),5)


 

# findConserved markers 
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
    verbose = TRUE  # Optional, for detailed output during processing
  )
  
  # Store the markers in the list, using the cluster number as the key
  all_cluster_markers[[paste0("Cluster_", cluster)]] <- markers
  # Store markers for this cluster as metadata
  #seurat_hdf5@meta.data[[paste0("markers_cluster_", cluster)]] <- markers
  
}


rm(markers)

# Inspect the results for a specific cluster (e.g., Cluster 3)
head(all_cluster_markers[["Cluster_0"]])
saveRDS(all_cluster_markers, file = "all_cluster_markers.rds")

# Add the all_cluster_markers list to the misc slot
seurat_hdf5@misc$all_cluster_markers <- all_cluster_markers

# Verify that the data has been added
names(seurat_hdf5@misc$all_cluster_markers)

#visualize top Features

FeaturePlot(seurat_hdf5,features = c("ST18"), min.cutoff = "q10")


 
#Cell type annotation

library(SingleR)
library(celldex)
library(pheatmap)


View(seurat_hdf5@meta.data)
ref <- celldex::HumanPrimaryCellAtlasData()

View(as.data.frame(colData(ref)))
pbmc_counts <- GetAssayData(seurat_hdf5, layer = 'counts')
pred <- SingleR(test = pbmc_counts,
                ref = ref,
                labels = ref$label.fine)

view(pred)
seurat_hdf5$singleR.labels <- pred$labels[match(rownames(seurat_hdf5@meta.data), rownames(pred))]
DimPlot(seurat_hdf5, reduction = 'umap', group.by = 'singleR.labels',label = TRUE)

#another way
pred1 <- SingleR(test = pbmc_counts,
               ref = ref,
            labels = ref$label.main)

view(pred1)

seurat_hdf5$singleR <- pred1$labels[match(rownames(seurat_hdf5@meta.data), rownames(pred1))]
p7 <- DimPlot(seurat_hdf5, reduction = 'umap', group.by = 'singleR',label=TRUE)+NoLegend()
p5|p7


 Annotation diagnostics ----------

# ...Based on the scores within cells 
view(pred)
pred$scores
plotScoreHeatmap(pred)

# ...Based on deltas across cells 

plotDeltaDistribution(pred)

tab <- table(Assigned=pred$labels, Clusters=seurat_hdf5$seurat_clusters)
h1 <- pheatmap(log10(tab+10), color = colorRampPalette(c('white','blue'))(10))

tab <- table(Assigned=pred1$labels, Clusters=seurat_hdf5$seurat_clusters)
h2 <- pheatmap(log10(tab+10), color = colorRampPalette(c('white','blue'))(10))

##
#for Cell.Type already in the metadata given

tab <- table(Assigned=seurat_hdf5@meta.data$Cell.Type, Clusters=seurat_hdf5$seurat_clusters)
h3 <- pheatmap(log10(tab+10), color = colorRampPalette(c('white','blue'))(10))

tab <- table(Assigned=seurat_hdf5@meta.data$cluster, Clusters=seurat_hdf5$seurat_clusters)
h4 <- pheatmap(log10(tab+10), color = colorRampPalette(c('white','blue'))(10))


 

# setting Idents as Seurat annotations provided (also a sanity check!)
Idents(seurat_hdf5) <- seurat_hdf5@meta.data$seurat_annotations
Idents(seurat_hdf5)

DimPlot(seurat_hdf5, reduction = 'umap', label = TRUE)

#find markers between conditions---------------------
seurat_hdf5$cell.type.cnd <- paste0(seurat_hdf5$seurat_annotations,'_', seurat_hdf5$Diagnosis)
View(seurat_hdf5@meta.data)
Idents(seurat_hdf5) <- seurat_hdf5$celltype.cnd

DimPlot(seurat_hdf5, reduction = 'umap', label = TRUE)

# find markers
b.interferon.response <- FindMarkers(seurat_hdf5, ident.1 = 'CD16 Mono_STIM', ident.2 = 'CD16 Mono_CTRL')

head(b.interferon.response)

# plotting conserved features vs DE features between conditions
head(markers_cluster3)


FeaturePlot(seurat_hdf5, features = c('FCGR3A', 'AIF1', 'IFIT1'), split.by = 'stim', min.cutoff = 'q10')


 Extract top markers for each cluster ----

top_markers_list <- lapply(all_cluster_markers, function(cluster_markers) {
  cluster_markers %>%
    filter(avg_log2FC > 1 & p_val_adj < 0.05) %>% # Adjust thresholds as needed
    arrange(desc(avg_log2FC)) %>%
    head(20)  # Select the top 20 markers
})

# Combine the results for all clusters into a single data frame
top_markers_df <- bind_rows(
  lapply(names(top_markers_list), function(cluster) {
    data.frame(cluster = cluster, top_markers_list[[cluster]])
  })
)

# Save the results for future reference
write.csv(top_markers_df, "top_markers_per_cluster.csv")

# Prepare a gene list for enrichment analysis ----
top_markers_list <- lapply(all_cluster_markers, function(markers_df) {
  # Add a new column for the average log2 fold change across conditions
  markers_df <- markers_df %>%
    mutate(avg_log2FC = (Control_avg_log2FC + AD_avg_log2FC) / 2)
  
  # Filter and sort genes based on adjusted p-value and fold change
  filtered_df <- markers_df %>%
    filter(max_pval < 0.05) %>%  # Consider only significant genes
    arrange(desc(avg_log2FC)) %>%  # Sort by average fold change
    slice_head(n = 10)  # Select the top 10 genes
  
  # Extract the row names (gene names)
  rownames(filtered_df)
})

top_markers_df <- bind_rows(
  lapply(names(top_markers_list), function(cluster) {
    data.frame(cluster = cluster, gene = top_markers_list[[cluster]])
  })
)
View(top_markers_df)

# Save to a CSV file
write.csv(top_markers_df, "top_markers_per_cluster.csv", row.names = FALSE)

# Combine all cluster gene lists into a single unique list for enrichment
gene_list <- unique(unlist(top_markers_list))

# View the resulting gene list
head(gene_list)


library(clusterProfiler)
library(org.Hs.eg.db)  # Use appropriate organism database

# Convert gene names to Entrez IDs
gene_list <- unique(top_markers_df$gene)
gene_list <- toupper(gene_list)

gene_ids <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Find unmapped genes
unmapped_genes <- setdiff(gene_list, gene_ids$SYMBOL)

# Print unmapped genes for inspection
print(unmapped_genes)


# GO enrichment analysis----
go_results <- enrichGO(
  gene = gene_ids$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",  # Biological process
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)


# KEGG pathway enrichment
kegg_results <- enrichKEGG(
  gene = gene_ids$ENTREZID,
  organism = 'hsa',  # Change to the relevant organism code (e.g., "mmu" for mouse)
  keyType = "kegg",
  pvalueCutoff = 0.05
)

# Visualize GO results
dotplot(go_results, showCategory = 10)+ ggtitle("GO Enrichment")  # Adjust showCategory as needed

# Visualize KEGG results
dotplot(kegg_results, showCategory = 1)+ ggtitle("KEGG Pathway Enrichment")

# Group by annotated cell types and extract top markers
celltype_markers <- FindMarkers(seurat_hdf5, group.by = "singleR.labels")

# Save cell type markers
write.csv(celltype_markers, "celltype_markers.csv")

# Visualize cell type marker expression
FeaturePlot(seurat_hdf5, features = c("marker_gene1", "marker_gene2"))


library(enrichR)

# Available databases
dbs <- listEnrichrDbs()

# Perform enrichment analysis
results <- enrichr(gene_list, c("GO_Biological_Process_2021", "KEGG_2021_Human"))

# View results for GO Biological Process
head(results$GO_Biological_Process_2021)

# Visualize results
barplot(as.numeric(results$KEGG_2021_Human$Adjusted.P.value[1:10]), 
        names.arg = results$KEGG_2021_Human$Term[1:10],
        las = 2,
        col = "skyblue",
        main = "Top KEGG Pathways")

write.table(gene_list, file = "marker_genes.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Visualize enriched GO terms with a barplot
barplot(go_results, showCategory = 10)

# Visualize enriched KEGG pathways
emapplot(kegg_results)






