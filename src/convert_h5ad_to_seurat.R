#!/usr/bin/env Rscript 
library(Seurat)
library(reticulate)
library(sceasy)

# Simple command-line argument parsing to accept --input and --output
args <- commandArgs(trailingOnly = TRUE)

print_usage <- function() {
	cat("convert_h5ad_to_seurat.R - convert an AnnData (.h5ad) file to a Seurat .rds file\n")
	cat("\nUsage:\n")
	cat("  Rscript src/convert_h5ad_to_seurat.R --input <file.h5ad> [--output <file.rds>]\n\n")
	cat("Options:\n")
	cat("  --input   Path to input .h5ad file (default: output/scanpy/bp_harmony_integrated.h5ad)\n")
	cat("  --output  Path to output .rds file (default: same as input but with .rds extension)\n")
	invisible(NULL)
}

if (length(args) == 0) {
	# No args provided: use existing defaults
	input_file <- "output/scanpy/combined_harmony_integrated.h5ad"
	output_file <- "output/scanpy/combined_harmony_integrated.rds"
} else {
	# parse args
	input_file <- NULL
	output_file <- NULL
	i <- 1
	while (i <= length(args)) {
		a <- args[i]
		if (a == "-h" || a == "--help") {
			print_usage()
			quit(status = 0)
		} else if (startsWith(a, "--input=")) {
			input_file <- sub("^--input=", "", a)
		} else if (a == "--input") {
			if ((i+1) <= length(args)) {
				input_file <- args[i+1]
				i <- i + 1
			} else {
				stop("--input requires a value")
			}
		} else if (startsWith(a, "--output=")) {
			output_file <- sub("^--output=", "", a)
		} else if (a == "--output") {
			if ((i+1) <= length(args)) {
				output_file <- args[i+1]
				i <- i + 1
			} else {
				stop("--output requires a value")
			}
		} else {
			stop(paste("Unknown argument:", a))
		}
		i <- i + 1
	}

	# set defaults if missing
	if (is.null(input_file)) {
		input_file <- "output/scanpy/bp_harmony_integrated.h5ad"
	}
	if (is.null(output_file)) {
		# replace .h5ad extension with .rds; if no .h5ad, just append .rds
		if (grepl("\\.h5ad$", input_file, ignore.case = TRUE)) {
			output_file <- sub("\\.h5ad$", ".rds", input_file, ignore.case = TRUE)
		} else {
			output_file <- paste0(input_file, ".rds")
		}
	}
}

# Ensure conda env is used for reticulate (if needed in the environment)
try({
	use_condaenv("scvi-tools", conda = "/opt/apps/mamba/24.3.0/bin/mamba")
}, silent = TRUE)

message("Converting: ", input_file, " -> ", output_file)

# perform conversion
sceasy::convertFormat(input_file, from = "anndata", to = "seurat", outFile = output_file)

seu <- readRDS(output_file)
Idents(seu) <- seu$batch

seu <- NormalizeData(seu)

seu <- FindVariableFeatures(seu, selection.method = "vst")

seu <- ScaleData(seu, features = VariableFeatures(seu), block.size = 200)

saveRDS(seu, file = output_file)

message("Done. Saved Seurat object to: ", output_file)
