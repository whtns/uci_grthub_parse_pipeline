#!/usr/bin/env Rscript

args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
  eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages({
  library('tidyverse')
  library('fs')
  library(clustifyr)
  library(Seurat)
  library(seuratTools)
})

seu <- readRDS(seu_path)

mynb <- readRDS(numbat_rds)

seu <- Seurat::RenameCells(seu, new.names = str_replace(colnames(seu), "\\.", "-"))

nb_meta <- mynb[["clone_post"]][,c("cell", "clone_opt", "GT_opt", "compartment_opt")] %>%
  dplyr::mutate(cell = str_replace(cell, "\\.", "-")) %>%
  tibble::column_to_rownames("cell")

seu <- Seurat::AddMetaData(seu, nb_meta)

seu <- seu[,!is.na(seu$clone_opt)]

filtered_seu <- seu[,(!seu$GT_opt != "")]


filter_non_tumor_cells_from_seu <- function(seu, mynb){

  seu <- seu




}

filtered_seu <- filter_non_tumor_cells_from_seu(seu, mynb)


saveRDS(filtered_seu, seu_path_out)


