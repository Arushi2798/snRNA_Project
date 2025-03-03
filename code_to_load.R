library(Seurat)
library(Matrix)
library(hdf5r)
library(tidyverse)
library(MAST)
library(presto)

setwd("D:/3rd sem/major_project")
hd5_object <- Read10X_h5(filename = "GSE174367_snRNA-seq_filtered_feature_bc_matrix.h5",
                         use.names = TRUE, unique.features = TRUE)

# Assuming 'hd5_object' is your sparse matrix
# Get the column names (barcodes)
barcodes <- colnames(hd5_object)

# Filter for barcodes needed
barcode_CT <- barcodes[grep(".*\\b(16|8|13|17|14|12|9|10|6|2|15|11|18|3)\\b$", barcodes)]

# Create a logical vector for keeping the filtered barcodes
keep_CT_columns <- colnames(hd5_object) %in% barcode_CT

# Filter the sparse matrix to keep only the desired columns
filtered_matrix_CT <- hd5_object[, keep_CT_columns]

# Check the dimensions of the filtered matrix
dim(filtered_matrix_CT)
head(filtered_matrix_CT)

# Create a Seurat object from the filtered matrix
# Create a Seurat object using the filtered dataset
seurat_hdf5 <- CreateSeuratObject(counts = filtered_matrix_CT, project = "AD", min.cells = 3, min.features = 200)


str(seurat_hdf5)

head(seurat_hdf5)
dim(seurat_hdf5)

# Delete un-needed variables
rm(barcode_CT, barcodes, keep_CT_columns, hd5_object,filtered_matrix_CT)
gc()  # Optional: Clean up memory

