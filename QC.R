#set origin identity to null
seurat_hdf5@meta.data$orig.ident = NULL
seurat_hdf5@assays$RNA@meta.data = data.frame(Geneid = rownames(seurat_hdf5), Symbol = rownames(seurat_hdf5), row.names = rownames(seurat_hdf5), stringsAsFactors = FALSE)

################################################################################
# Step 02: Quality Control
################################################################################

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

# Filter the Seurat object based on the suggested thresholds
seurat_hdf5 <- subset(x = seurat_hdf5, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 5)

#check the new dimension
dim(seurat_hdf5)

#removing all mitochondrial genes
seurat_hdf5 <- seurat_hdf5[!grepl("^MT-", rownames(seurat_hdf5)),]

dim(seurat_hdf5) #dimension change from 35265 X 47856 to 35245 X 47856

# meta1 = seurat_hdf5@meta.data
# write.table(meta1, file=paste0('meta.seurat.raw.tsv'), quote = FALSE, sep = '\t', row.names = FALSE)
# meta1 = meta1[rev(order(meta1$Diagnosis)), ]
# meta1$sid = unsplit(lapply(split(meta1, meta1$Diagnosis), function(x) paste0(c(AD = 'AD', Control = 'Ct')[x$Diagnosis[1]], as.integer(factor(x$SampleID, level = unique(x$SampleID))))), meta1$Diagnosis)
# meta1$sid = factor(meta1$sid, levels = unique(meta1$sid))
# theme_bw2 = function() theme_bw() + theme(panel.grid = element_blank())
# p1 = ggplot(meta1, aes(x = sid, y = nFeature_RNA, fill = Diagnosis)) + geom_violin(width = 0.8, scale = 'width') + theme_bw2() + theme(legend.position = 'none', axis.text.x = element_text(angle = 45, size = 6, hjust = 1), axis.text.y = element_text(size = 7)) + xlab(NULL)
# p2 = ggplot(meta1, aes(x = sid, y = nCount_RNA, fill = Diagnosis)) + geom_violin(width = 0.8, scale = 'width') + theme_bw2() + theme(legend.position = 'none', axis.text.x = element_text(angle = 45, size = 6, hjust = 1), axis.text.y = element_text(size = 7)) + xlab(NULL)
# print(p1 / p2)

# seurat_hdf5@meta.data$SampleID = factor(seurat_hdf5@meta.data$SampleID)
# p3 = FeatureScatter(seurat_hdf5, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size=0.6, group.by='SampleID') + theme(title=element_text(size=10), legend.text=element_text(size=10), legend.title=element_text(size=12))
# print(p3)