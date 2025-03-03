#loading by directly selecting the file
metadata <- read.csv(file.choose(), row.names = 1)

# Add metadata to Seurat object
seurat_hdf5 <- AddMetaData(object=seurat_hdf5, metadata = metadata)

#to check the dimensions of the Seurat object after adding metadata
dim(seurat_hdf5)

#Check that the metadata was added correctly
head(seurat_hdf5@meta.data)

#check the dimension of the metadata
dim(seurat_hdf5@meta.data)

# Find barcodes in the Seurat object that are missing from the metadata
missing_barcodes <- setdiff(colnames(seurat_hdf5), rownames(metadata))
length(missing_barcodes)
head(missing_barcodes)  # View the first few missing barcodes