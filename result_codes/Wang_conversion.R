# ==============================================================================
# Script: 01_Wang_Data_Conversion.R
# Purpose: Extract raw count matrices and features from the Zebrafish Cell 
#          Landscape (ZCL) Seurat object.
# ==============================================================================

# Load required libraries
library(Seurat)
library(Matrix)

# 1. Load the original RData file provided by Wang et al.
# (This contains the Seurat object typically named 'pbmc' or similar)
print("Loading original ZCL dataset...")
load("ZCDL.rdata")

# 2. Extract the raw counts matrix
# Using 'counts' to ensure we get un-normalized, raw integer data
print("Extracting raw count matrix...")
counts_matrix <- GetAssayData(pbmc, layer = "counts")

# 3. Export the matrix in the universal Matrix Market format (.mtx)
print("Saving matrix to .mtx format...")
writeMM(counts_matrix, file = "Zebrafish_counts.mtx")

# 4. Export the gene names (features) and cell barcodes
print("Saving genes and barcodes...")
write.table(rownames(counts_matrix), file = "Zebrafish_genes.tsv", 
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(colnames(counts_matrix), file = "Zebrafish_barcodes.tsv", 
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

print("Extraction complete. Ready for Python processing.")