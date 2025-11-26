#!/usr/bin/env Rscript
library(Seurat)
library(Matrix)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    stop("Usage: Rscript src/seurat_scale_data.R INPUT OUTPUT")
}   

seurat_input <- args[1]
seurat_output <- args[2]

if (!file.exists(seurat_input)) {
    stop(paste("Seurat RDS input not found:", seurat_input))
}
seu <- readRDS(seurat_input)

seu <- FindVariableFeatures(seu, selection.method = "vst")
seu <- ScaleData(seu, features = VariableFeatures(seu), block.size = 200)
# Save the modified Seurat object
saveRDS(seu, seurat_output)