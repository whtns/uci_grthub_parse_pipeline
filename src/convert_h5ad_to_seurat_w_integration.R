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
  make_option(c("-o","--output"), type = "character", default = file.path(dirname("/dfs9/ucightf-lab/projects/OlabR/250922_0925Bio-10_OlabR_parse/output/scanpy/combined_harmony_integrated.h5ad"), "combined_harmony_integrated.rds"),
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

# Helper: safely extract iterable keys from a Python mapping-like object (KeysView, dict_keys, etc.)
safe_py_keys <- function(pyobj){
  if(is.null(pyobj)) return(character(0))
  keys <- NULL
  # 1) Prefer asking Python to return a list of stringified keys (this handles KeysView reliably)
  try({
    keys <- reticulate::py_to_r(reticulate::py_eval("list(map(str, obm.keys()))", envir = list(obm = pyobj)))
  }, silent = TRUE)
  # 2) Fallback: try to call the object's keys() method via py_get_attr + py_call
  if(is.null(keys)){
    try({
      if(reticulate::py_has_attr(pyobj, "keys")){
        f <- reticulate::py_get_attr(pyobj, "keys")
        kpy <- reticulate::py_call(f)
        keys <- reticulate::py_to_r(kpy)
      }
    }, silent = TRUE)
  }
  # 3) Fallback: coerce whole object to R and use names
  if(is.null(keys)){
    tmp <- NULL
    try({ tmp <- reticulate::py_to_r(pyobj) }, silent = TRUE)
    if(is.list(tmp)) keys <- names(tmp)
  }
  # 3) Fallback: coerce whole object to R and use names
  if(is.null(keys)){
    tmp <- NULL
    try({ tmp <- reticulate::py_to_r(pyobj) }, silent = TRUE)
    if(is.list(tmp)) keys <- names(tmp)
  }
  if(is.null(keys)) return(character(0))
  if(is.list(keys)) keys <- unlist(keys)
  as.character(keys)
}

try_convert_and_load <- function(infile, h5seurat_out = NULL){
  # Try to convert .h5ad -> .h5seurat then load
  if(!requireNamespace("SeuratDisk", quietly=TRUE)) stop("SeuratDisk required for conversion but not available")
  if(is.null(h5seurat_out)) h5seurat_out <- sub("\\.h5ad$", ".h5seurat", infile)
  message("Converting .h5ad -> .h5seurat (this may take a moment)...")
  tryCatch({
    SeuratDisk::Convert(infile, dest = "h5seurat", overwrite = TRUE)
    message("Load h5seurat: ", h5seurat_out)
    seu <- SeuratDisk::LoadH5Seurat(h5seurat_out, assays = "RNA")
    return(list(seu = seu, h5seurat = h5seurat_out))
  }, error = function(e){
    # Improve diagnostics for common 'Ambiguous assays' failure during LoadH5Seurat
    msg <- conditionMessage(e)
    message("SeuratDisk conversion/load failed: ", msg)

    # If the converted file exists, try to inspect assays inside it to provide guidance
    if(file.exists(h5seurat_out)){
      message("Inspecting generated .h5seurat file to list assays...")
      assays_found <- NULL
      # Try rhdf5 first, then hdf5r
      if(requireNamespace("rhdf5", quietly=TRUE)){
        try({
          ls <- rhdf5::h5ls(h5seurat_out)
          # look for groups under /assays
          if(any(ls$group == "/assays")){
            assays_found <- unique(sub("/assays/", "", ls$name[ls$group == "/assays"]))
          }
        }, silent = TRUE)
      }
      if(is.null(assays_found) && requireNamespace("hdf5r", quietly=TRUE)){
        try({
          f <- hdf5r::H5File$new(h5seurat_out, mode = "r")
          if("assays" %in% names(f)) assays_found <- names(f[["assays"]])
          f$close_all()
        }, silent = TRUE)
      }

      if(!is.null(assays_found) && length(assays_found) > 0){
        message(sprintf("Assays found in %s: %s", h5seurat_out, paste(assays_found, collapse=", ")))
        message("SeuratDisk may be unable to choose a default assay automatically. You can:")
        message("  - specify which assay to load by using SeuratDisk::LoadH5Seurat(..., assays = \"<ASSAY_NAME>\") if available,")
        message("  - or set the desired default assay in the resulting Seurat object: Seurat::DefaultAssay(seu) <- \"<ASSAY_NAME>\"")
      } else {
        message("Could not enumerate assays inside the .h5seurat file (rhdf5/hdf5r may not be installed or file structure unexpected).")
      }
    }

    # Re-raise the original error so outer try() can fall back to the reticulate path
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
  # If sparse matrix object, prefer keeping it sparse to avoid large dense allocations.
  # Seurat's CreateSeuratObject accepts sparse matrices (Matrix package) so we avoid coercion here.
  if(inherits(X_r, "dgCMatrix") || inherits(X_r, "sparseMatrix")){
    message("AnnData X is a sparse matrix; keeping as sparse to avoid large memory usage")
    # leave X_r as-is (dgCMatrix) and do not coerce to dense
  }

  # Build Seurat object
  seu <- Seurat::CreateSeuratObject(counts = t(X_r)) # note: anndata X usually cells x features, Seurat expects features x cells -> transpose

  # Try to set cell names from AnnData (ad.obs_names)
  try({
    cell_names <- reticulate::py_to_r(ad$obs_names)
    if(!is.null(cell_names) && length(cell_names) == ncol(seu)){
      colnames(seu) <- cell_names
      message("Set Seurat cell names from AnnData obs_names")
    } else {
      message("AnnData obs_names missing or length mismatch; keeping default Seurat cell names")
    }
  }, silent = TRUE)

  # Add obs (cell metadata)
  try({
    obs_df <- reticulate::py_to_r(ad$obs)
    if(is.data.frame(obs_df)){
      # Ensure metadata rownames match Seurat cell names
      try({
        if(exists("cell_names") && length(cell_names) == nrow(obs_df)){
          rownames(obs_df) <- cell_names
        }
      }, silent = TRUE)
      seu <- Seurat::AddMetaData(seu, metadata = obs_df)
      message("Added cell metadata from AnnData.obs")
    }
  })

  # Add reductions from .obsm
  # Add reductions from .obsm with explicit checks and informative messages
  obsm <- ad$obsm
    if(!is.null(obsm)){
      # Use helper to obtain a safe character vector of keys (handles KeysView/dict_keys)
      keys <- safe_py_keys(obsm)
      message(sprintf("Found obsm keys in AnnData: %s", paste(keys, collapse=", ")))
      if(length(keys) > 0){
        for(k in keys){
        message(sprintf("Processing obsm key: %s", k))
        mat <- NULL
        try({ mat <- reticulate::py_to_r(obsm[[k]]) })
        if(is.null(mat)){
          message(sprintf("  - Could not convert obsm[%s] to R matrix/data.frame", k))
          next
        }
        # orient so rows = cells
        if(nrow(mat) == ncol(seu)) emb <- mat
        else if(ncol(mat) == ncol(seu)) emb <- t(mat)
        else {
          message(sprintf("  - Skipping %s: dimension mismatch (mat %dx%d vs cells=%d)", k, nrow(mat), ncol(mat), ncol(seu)))
          next
        }
        # Ensure rownames of embedding match cell names
        if(exists("cell_names") && length(cell_names) == nrow(emb)){
          rownames(emb) <- cell_names
        }
        # choose key name for the internal Key parameter
        key <- ""
        if(grepl("pca", k, ignore.case=TRUE)) key <- "PC_"
        else if(grepl("umap", k, ignore.case=TRUE)) key <- "UMAP_"
        else if(grepl("tsne", k, ignore.case=TRUE)) key <- "tSNE_"
        else key <- toupper(substr(k,1,3)); key <- paste0(key, "_")

        # Create DimReduc object and attach under a cleaned name (remove leading X_ if present)
        red_name <- sub("^X_", "", k)
        dr_obj <- Seurat::CreateDimReducObject(emb = emb, key = key, assay = Seurat::DefaultAssay(seu))
        seu[[red_name]] <- dr_obj
        message(sprintf("  - Attached reduction as '%s' (dims=%d)", red_name, ncol(emb)))
      }
    }
  }

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

# Ensure reductions from the input .h5ad are attached to the Seurat object when possible
add_reductions_from_h5ad <- function(seu, infile){
  if(!requireNamespace("reticulate", quietly=TRUE)){
    message("reticulate not available; cannot import additional reductions from .h5ad")
    return(seu)
  }
  anndata <- NULL
  try({ anndata <- reticulate::import("anndata") }, silent = TRUE)
  if(is.null(anndata)){
    message("Could not import anndata via reticulate; skipping additional reductions")
    return(seu)
  }
  ad <- anndata$read_h5ad(infile)
  obsm <- ad$obsm
  if(is.null(obsm)) return(seu)
  # Retrieve keys reliably using safe_py_keys() which handles KeysView/dict_keys
  keys <- safe_py_keys(obsm)
  if(length(keys) == 0) return(seu)

  message(sprintf("add_reductions_from_h5ad: found obsm keys: %s", paste(keys, collapse=", ")))
  existing <- Seurat::Reductions(seu)
  # try to extract cell names from ad$obs_names
  cell_names <- NULL
  try({ cell_names <- reticulate::py_to_r(ad$obs_names) }, silent = TRUE)
  for(k in keys){
    # normalize name (remove leading X_ which is common in AnnData)
    red_name <- sub("^X_", "", k)
    if(red_name %in% existing) next
    mat <- NULL
    try({ mat <- reticulate::py_to_r(obsm[[k]]) })
    if(is.null(mat)){
      message(sprintf("  - could not convert obsm[%s], skipping", k)); next
    }
    if(nrow(mat) == ncol(seu)) emb <- mat
    else if(ncol(mat) == ncol(seu)) emb <- t(mat)
    else { message(sprintf("  - skipping %s: dim mismatch (mat %dx%d vs cells=%d)", k, nrow(mat), ncol(mat), ncol(seu))); next }
    if(!is.null(cell_names) && length(cell_names) == nrow(emb)) rownames(emb) <- cell_names
    key <- ""
    if(grepl("pca", k, ignore.case=TRUE)) key <- "PC_"
    else if(grepl("umap", k, ignore.case=TRUE)) key <- "UMAP_"
    else if(grepl("tsne", k, ignore.case=TRUE)) key <- "tSNE_"
    else key <- toupper(substr(k,1,3)); key <- paste0(key, "_")
    dr_obj <- Seurat::CreateDimReducObject(emb = emb, key = key, assay = Seurat::DefaultAssay(seu))
    seu[[red_name]] <- dr_obj
    message(sprintf("  - attached reduction '%s' (dims=%d)", red_name, ncol(emb)))
  }
  return(seu)
}

# Try to attach any missing reductions from the source .h5ad
seu <- add_reductions_from_h5ad(seu, infile)

# Save RDS
dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
saveRDS(seu, file = outfile)
message("Saved Seurat object to: ", outfile)

invisible(NULL)
