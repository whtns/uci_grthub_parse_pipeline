#!/usr/bin/env Rscript
# Converts a Scanpy/AnnData h5ad file to a Seurat object and saves as .rds
# Preserves dimensionality reductions when possible (PCA, UMAP, tSNE, etc.)
# Usage:
# Rscript src/convert_to_seurat_w_integration.R.R --input <path/to/combined_harmony_integrated.h5ad> --output <path/to/output.rds>

suppressMessages({
  if(!requireNamespace("optparse", quietly=TRUE)) {
    install.packages("optparse", repos = "https://cloud.r-project.org")
  }
  library(optparse)
})

option_list <- list(
  make_option(c("-i","--input"), type = "character", default = "/dfs9/ucightf-lab/projects/OlabR/250922_0925Bio-10_OlabR_parse/output/scanpy/combined_harmony_integrated.h5ad",
              help = "Path to input .h5ad file", metavar = "character"),
  make_option(c("-o","--output"), type = "character", default = file.path(dirname("/dfs9/ucightf-lab/projects/OlabR/250922_0925Bio-10_OlabR_parse/output/scanpy/combined_harmony_integrated.h5ad"), "combined_harmony_integrated.seurat.rds"),
              help = "Path to output .rds file (Seurat object)", metavar = "character"),
  make_option(c("--h5seurat"), type = "character", default = NULL,
              help = "Optional path to write intermediate .h5seurat (if omitted will be <input>.h5seurat)", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

infile <- opt$input
outfile <- opt$output
h5seurat_out <- opt$h5seurat

message(sprintf("Input: %s", infile))
message(sprintf("Output RDS: %s", outfile))

if(!file.exists(infile)) stop("Input file does not exist: ", infile)

# Load required packages (do not auto-install heavy Bioconductor/CRAN deps beyond optparse above)
if(!requireNamespace("Seurat", quietly=TRUE)) stop("Package 'Seurat' is required but not installed.")
if(!requireNamespace("SeuratDisk", quietly=TRUE) && !requireNamespace("reticulate", quietly=TRUE)) {
  stop("Either 'SeuratDisk' or 'reticulate' (with 'anndata') is required to read .h5ad files. Please install them.")
}

library(Seurat)
reticulate::use_condaenv("scvi-tools")

try_convert_and_load <- function(infile, h5seurat_out = NULL){
  # Try to convert .h5ad -> .h5seurat then load
  if(!requireNamespace("SeuratDisk", quietly=TRUE)) stop("SeuratDisk required for conversion but not available")
  if(is.null(h5seurat_out)) h5seurat_out <- sub("\\.h5ad$", ".h5seurat", infile)
  message("Converting .h5ad -> .h5seurat (this may take a moment)...")
  tryCatch({
    SeuratDisk::Convert(infile, dest = "h5seurat", overwrite = TRUE)
    message("Load h5seurat: ", h5seurat_out)
    seu <- SeuratDisk::LoadH5Seurat(h5seurat_out)
    return(list(seu = seu, h5seurat = h5seurat_out))
  }, error = function(e){
    message("SeuratDisk conversion/load failed: ", conditionMessage(e))
    stop(e)
  })
}

fallback_read_h5ad_via_reticulate <- function(infile){
  # As a fallback, use reticulate + anndata to import reductions and counts and build a Seurat object
  if(!requireNamespace("reticulate", quietly=TRUE)) stop("reticulate required for fallback but not installed")
  message("Falling back to reticulate::anndata read_h5ad -> building Seurat object")
  anndata <- reticulate::import("anndata")
  ad <- anndata$read_h5ad(infile)

  # Extract X (counts or matrix)
  X_r <- NULL
  try({
    X_py <- ad$X
    X_r <- reticulate::py_to_r(X_py)
  }, silent = TRUE)

  if(is.null(X_r)) stop("Could not extract main matrix from AnnData via reticulate")

  # If sparse matrix object, convert to dense matrix with caution
  if(inherits(X_r, "dgCMatrix") || inherits(X_r, "sparseMatrix")){
    message("Converting sparse matrix to dense matrix for Seurat (may use lots of memory)")
    X_r <- as.matrix(X_r)
  }

  # Build Seurat object
  seu <- Seurat::CreateSeuratObject(counts = t(X_r)) # note: anndata X usually cells x features, Seurat expects features x cells -> transpose

  # Add obs (cell metadata)
  try({
    obs_df <- reticulate::py_to_r(ad$obs)
    if(is.data.frame(obs_df)){
      rownames(obs_df) <- rownames(obs_df)
      seu <- Seurat::AddMetaData(seu, metadata = obs_df)
    }
  }, silent = TRUE)

  # Add reductions from .obsm
  try({
    obsm <- ad$obsm
    if(!is.null(obsm)){
      keys <- reticulate::py_list_names(obsm)
      if(length(keys) > 0){
        for(k in keys){
          mat <- reticulate::py_to_r(obsm[[k]])
          # Ensure matrix has rows = cells
          if(nrow(mat) == ncol(seu) || ncol(mat) == ncol(seu)){
            if(ncol(mat) == ncol(seu)){
              emb <- mat
            } else {
              emb <- t(mat)
            }
            # choose key name
            key <- ""
            if(grepl("pca", k, ignore.case=TRUE)) key <- "PC_"
            else if(grepl("umap", k, ignore.case=TRUE)) key <- "UMAP_"
            else if(grepl("tsne", k, ignore.case=TRUE)) key <- "tSNE_"
            else key <- toupper(substr(k,1,3)); key <- paste0(key, "_")

            dr_obj <- Seurat::CreateDimReducObject(emb = emb, key = key, assay = Seurat::DefaultAssay(seu))
            seu[[k]] <- dr_obj
            message(sprintf("Added reduction: %s (dims=%d)", k, ncol(emb)))
          } else {
            message(sprintf("Skipping reduction %s: dimension mismatch (matrix %dx%d vs cells %d)", k, nrow(mat), ncol(mat), ncol(seu)))
          }
        }
      }
    }
  }, silent = TRUE)

  return(seu)
}

# Main
seu <- NULL
conversion_attempted <- FALSE
if(requireNamespace("SeuratDisk", quietly=TRUE)){
  conversion_attempted <- TRUE
  try({
    res <- try_convert_and_load(infile, h5seurat_out)
    seu <- res$seu
  }, silent = TRUE)
}

if(is.null(seu)){
  message("Primary conversion pathway failed or not available. Trying Seurat::ReadH5AD (if available) ...")
  if(exists("ReadH5AD", where = asNamespace("Seurat"))){
    try({
      seu <- Seurat::ReadH5AD(infile)
    }, silent = TRUE)
  }
}

if(is.null(seu)){
  message("Seurat::ReadH5AD unavailable or failed. Trying fallback via reticulate/anndata ...")
  seu <- fallback_read_h5ad_via_reticulate(infile)
}

if(is.null(seu)) stop("Failed to produce a Seurat object from the input .h5ad file.")

message("Seurat object created. Checking dimensional reductions...")
rels <- Seurat::Reductions(seu)
if(length(rels) == 0){
  message("No reductions detected on the Seurat object.")
} else {
  message(sprintf("Detected reductions: %s", paste(rels, collapse = ", ")))
}

# Save RDS
dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
saveRDS(seu, file = outfile)
message("Saved Seurat object to: ", outfile)

invisible(NULL)
