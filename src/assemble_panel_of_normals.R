#!/usr/bin Rscript

args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
  eval(parse(text = args[[i]]))
}

library(fs)
library(tidyverse)
library(numbat)
library(Seurat)
library(readr)
library(magrittr)
library(glue)

# Rprof(rprof_out)

study = "yang"
seu_paths <- dir_ls(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/seurat/"), glob = "*infercnv*_seu.rds") %>% 
	str_subset(".*dropped.*", negate = TRUE)

bad_cell_types = c("RPCs", "Late RPCs", c("Red Blood Cells", "Microglia", "Muller Glia", "RPE", "Horizontal Cells", "Rod Bipolar Cells", "Pericytes", "Bipolar Cells", "Astrocytes", "Endothelial", "Schwann", "Fibroblasts"))


filter_to_bad_cells <- function(seu_path, bad_cell_types){
	# browser()
	seu <- readRDS(seu_path)
	
	seu <- seu[,seu$type %in% bad_cell_types]
	
	return(seu)
}

safe_filter_to_bad_cells <- safely(filter_to_bad_cells)

seus <- seu_paths %>% 
	map(safe_filter_to_bad_cells, bad_cell_types) %>% 
	map("result") %>% 
	unlist()

myseu <- reduce(seus, merge)

num_out_cells <- sum(myseu$type %in% bad_cell_types)

normal_reference_mat <- 
	myseu %>% 
	GetAssayData(slot = "counts")

normal_cell_annot <-
	myseu@meta.data[c("type")] %>% 
	tibble::rownames_to_column("cell") %>% 
	dplyr::rename(group = type) %>% 
	identity()

ref_internal = numbat::aggregate_counts(normal_reference_mat, normal_cell_annot)

saveRDS(ref_internal, glue("~/Homo_sapiens/numbat/{study}_ref.rds"))
