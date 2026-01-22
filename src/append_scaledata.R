#!/usr/bin/env Rscript
library(Seurat)
library(Matrix)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    stop("Usage: Rscript src/append_scaledata.R PATH_TO_SCALED_DATA_DIRECTORY [PREFIX] [SEURAT_RDS]\n\nExample: Rscript src/append_scaledata.R output/scanpy scaled_data output/Seurat5Shiny/seurat5.rds")
}

# Positional args: data directory, optional prefix (default 'scaled_data'), optional seurat rds path
data_dir <- args[1]
prefix <- if (length(args) >= 2 && nzchar(args[2])) args[2] else "scaled_data"
# Optional third arg: path to Seurat RDS to modify (default preserved for backward compatibility)
seurat_input <- if (length(args) >= 3 && nzchar(args[3])) args[3] else "output/Seurat5Shiny/seurat5.rds"
if (!dir.exists(data_dir)) {
    stop(paste("Directory not found:", data_dir))
}

# Construct filenames based on provided prefix
mtx_file <- file.path(data_dir, paste0(prefix, ".mtx"))
barcodes_file <- file.path(data_dir, paste0(prefix, "_barcodes.csv"))
genes_file <- file.path(data_dir, paste0(prefix, "_genes.csv"))

# Load the scaled data matrix
scaled_data <- readMM(mtx_file)

# Load cell and gene names
cell_barcodes <- read.csv(barcodes_file)$Barcode
gene_names <- read.csv(genes_file)$GeneName

# Set row and column names for the matrix
rownames(scaled_data) <- gene_names
colnames(scaled_data) <- cell_barcodes

# Create a Seurat object using the scaled data
# Note: If you have raw counts or normalized data, you'd typically create the Seurat object with that
# and then use SetAssayData for the scaled data.
# For directly creating with scaled data, ensure you understand the implications for downstream analyses
# that might expect raw counts or normalized data in the 'counts' or 'data' slot.

if (!file.exists(seurat_input)) {
    stop(paste("Seurat RDS file not found:", seurat_input))
}
seurat_obj <- readRDS(seurat_input)

# Determine target assay (use default assay if present)
assay_name <- tryCatch(DefaultAssay(seurat_obj), error = function(e) NULL)
if (is.null(assay_name) || assay_name == "") assay_name <- NULL

# Try to set the scaled data using the newer 'layer' argument first.
# Use `try(..., silent=TRUE)` (no `finally` arg) and only convert sparse->dense
# if it is required to succeed, to avoid large memory allocations when possible.
res <- try(SetAssayData(object = seurat_obj, assay = assay_name, layer = 'scale.data', new.data = scaled_data), silent = TRUE)
if (!inherits(res, "try-error")) {
    seurat_obj <- res
} else {
    message("layer-based SetAssayData failed, falling back to slot-based attempt.")
    # Attempt slot-based call without converting if scaled_data is sparse
    if (inherits(scaled_data, "sparseMatrix") || inherits(scaled_data, "dgTMatrix") || inherits(scaled_data, "dgCMatrix")) {
        rownames(scaled_data) <- rownames(GetAssayData(seurat_obj, layer = "data"))
        colnames(scaled_data) <- colnames(GetAssayData(seurat_obj, layer = "data"))
        res2 <- try(SetAssayData(object = seurat_obj, assay = assay_name, slot = 'scale.data', new.data = scaled_data), silent = TRUE)
        if (!inherits(res2, "try-error")) {
            seurat_obj <- res2
        } else {
            message("Setting sparse 'scale.data' failed: attempting dense conversion. This may require a lot of memory.")
            seurat_obj <- SetAssayData(object = seurat_obj, layer = "scale.data", as.matrix(scaled_data), assay = assay_name)
        }
    } else {
        # scaled_data already dense; set directly
        seurat_obj <- SetAssayData(object = seurat_obj, assay = assay_name, slot = 'scale.data', new.data = scaled_data)
    }
}

# Save the modified Seurat object (overwrite original)
# Save the modified Seurat object next to the input file with a `_mod.rds` suffix
seurat_dir <- dirname(seurat_input)
seurat_base <- basename(seurat_input)
# seurat_mod_name <- sub("\\.rds$", "_mod.rds", seurat_base, ignore.case = TRUE)
# seurat_path <- file.path(seurat_dir, seurat_mod_name)
saveRDS(seurat_obj, seurat_input)
cat("Wrote updated Seurat object to:", seurat_input, "\n")