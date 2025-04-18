#######################################################################
#library needed in whole project
#######################################################################
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


################################################################################
# Step 0: Data Retrieval
################################################################################

#read in the data
hd5_object <- Read10X_h5(filename = "GSE174367_snRNA-seq_filtered_feature_bc_matrix.h5",
                         use.names = TRUE, unique.features = TRUE)

#if needed then subset the data for specific samples
# Get the column names (barcodes)
barcodes <- colnames(hd5_object)

# Filter for barcodes needed
barcode_CT <- barcodes[grep(".*\\b(sample_no_for_specific_samples)\\b$", barcodes)]

# Create a logical vector for keeping the filtered barcodes
keep_CT_columns <- colnames(hd5_object) %in% barcode_CT

# Filter the sparse matrix to keep only the desired columns
filtered_matrix_CT <- hd5_object[, keep_CT_columns]

# Create a Seurat object using the filtered dataset
seurat_hdf5 <- CreateSeuratObject(counts = filtered_matrix_CT, project = "AD", min.cells = 3, min.features = 200)

#to view the structure of the data
str(seurat_hdf5)

head(seurat_hdf5)

# Check the dimensions of the data
dim(seurat_hdf5)

# Delete unnecessary variables
rm(barcode_CT, barcodes, keep_CT_columns, hd5_object,filtered_matrix_CT)
gc()#  Clean up memory

#loading by directly selecting the file for metadata
metadata <- read.csv(file.choose(), row.names = 1)

#######################################################################
# Add metadata to Seurat object
#######################################################################
seurat_hdf5 <- AddMetaData(object=seurat_hdf5, metadata = metadata)

#Check that the metadata was added correctly
head(seurat_hdf5@meta.data)

#check the dimension of the metadata
dim(seurat_hdf5@meta.data)

# Find barcodes in the Seurat object that are missing from the metadata
missing_barcodes <- setdiff(colnames(seurat_hdf5), rownames(metadata))
length(missing_barcodes)
head(missing_barcodes)  # View the first few missing barcodes

#drop the cells if <5% data is missing

# Identify missing barcodes
missing_barcodes <- setdiff(colnames(seurat_hdf5), rownames(metadata))

# Filter Seurat object to keep only barcodes that exist in metadata
seurat_hdf5<- subset(seurat_hdf5, cells = setdiff(Cells(seurat_hdf5), missing_barcodes))

# Check the new dimensions
dim(seurat_hdf5)
