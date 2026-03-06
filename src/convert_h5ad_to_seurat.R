#!/usr/bin/env Rscript
library(Seurat)
library(anndataR)

args <- commandArgs(trailingOnly = TRUE)

print_usage <- function() {
  cat(paste0(
    "convert_h5ad_to_seurat.R - convert an AnnData (.h5ad) file to Seurat\n",
    "\nUsage:\n",
    "  Rscript src/convert_h5ad_to_seurat.R",
    " --data <data.h5ad> [--output <file.rds>]\n\n",
    "Options:\n",
    "  --data    Path to input .h5ad file\n",
    "  --output  Path to output .rds file",
    " (default: replaces .h5ad extension)\n"
  ))
  invisible(NULL)
}

data_file   <- NULL
output_file <- NULL

if (length(args) == 0) {
  data_file   <- "output/scanpy/combined_harmony_integrated.h5ad"
  output_file <- "output/scanpy/combined_harmony_integrated.rds"
} else {
  i <- 1
  while (i <= length(args)) {
    a <- args[i]
    if (a == "-h" || a == "--help") {
      print_usage()
      quit(status = 0)
    } else if (startsWith(a, "--data=")) {
      data_file <- sub("^--data=", "", a)
    } else if (a == "--data") {
      if ((i + 1) <= length(args)) {
        data_file <- args[i + 1]
        i <- i + 1
      } else {
        stop("--data requires a value")
      }
    } else if (startsWith(a, "--output=")) {
      output_file <- sub("^--output=", "", a)
    } else if (a == "--output") {
      if ((i + 1) <= length(args)) {
        output_file <- args[i + 1]
        i <- i + 1
      } else {
        stop("--output requires a value")
      }
    } else {
      stop(paste("Unknown argument:", a))
    }
    i <- i + 1
  }

  if (is.null(data_file)) {
    data_file <- "output/scanpy/combined_harmony_integrated.h5ad"
  }
  if (is.null(output_file)) {
    if (grepl("\\.h5ad$", data_file, ignore.case = TRUE)) {
      output_file <- sub("\\.h5ad$", ".rds", data_file, ignore.case = TRUE)
    } else {
      output_file <- paste0(data_file, ".rds")
    }
  }
}

message("Converting (all genes): ", data_file, " -> ", output_file)

seu <- anndataR::read_h5ad(data_file, as = "Seurat")

if ("batch" %in% colnames(seu[[]])) {
  Idents(seu) <- seu$batch
}

saveRDS(seu, file = output_file)

message("Done. Saved Seurat object to: ", output_file)
