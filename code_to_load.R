#library needed in whole project
library(Seurat)
library(Matrix)
library(hdf5r)
library(MAST)
library(presto)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(dplyr)
library(patchwork)
library(ggrepel)
library(scran)
library(glmGamPoi)
library(future)

#set working directory
setwd("path_to_directory")

#read in the data
hd5_object <- Read10X_h5(filename = "GSE174367_snRNA-seq_filtered_feature_bc_matrix.h5",
                         use.names = TRUE, unique.features = TRUE)

#subset the data for specific samples
# Get the column names (barcodes)
barcodes <- colnames(hd5_object)

# Filter for barcodes needed
barcode_CT <- barcodes[grep(".*\\b(sample_no_for_specific_samples)\\b$", barcodes)]

# Create a logical vector for keeping the filtered barcodes
keep_CT_columns <- colnames(hd5_object) %in% barcode_CT

# Filter the sparse matrix to keep only the desired columns
filtered_matrix_CT <- hd5_object[, keep_CT_columns]

# Check the dimensions of the filtered matrix
dim(filtered_matrix_CT)

# Create a Seurat object using the filtered dataset
seurat_hdf5 <- CreateSeuratObject(counts = filtered_matrix_CT, project = "AD", min.cells = 3, min.features = 200)

str(seurat_hdf5)

head(seurat_hdf5)
dim(seurat_hdf5)

# Delete unnecessary variables
rm(barcode_CT, barcodes, keep_CT_columns, hd5_object,filtered_matrix_CT)
  #  Clean up memory