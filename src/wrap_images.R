#!/usr/bin/env Rscript

library('tidyverse')
library('fs')
library('readxl')
library(numbat)
library(Seurat)
library(seuratTools)
library(infercnv)


plot_numbat_genotype_distribution <- function(numbat_seu_path) {
    # browser()
    seu <- readRDS(numbat_seu_path)

    seu@meta.data$clone_opt <- factor(seu@meta.data$clone_opt)

    DimPlot(seu, group.by = "clone_opt")

    genotype_per_cluster <- janitor::tabyl(seu@meta.data, clone_opt, gene_snn_res.0.2) %>%
        tidyr::pivot_longer(-c("clone_opt"), names_to = "cluster", values_to = "num_cell") %>%
        dplyr::group_by(cluster) %>%
        dplyr::mutate(percent_cell = num_cell/sum(num_cell))

    # percent cell  ------------------------------
    genotype_distribution = ggplot(genotype_per_cluster, aes(percent_cell, cluster, fill = clone_opt)) +
        geom_col() +
        theme_bw() +
        labs(title = numbat_seu_path, y = "Percent cells")

    clone_dimplot <- DimPlot(seu, group.by = "clone_opt")

    cluster_dimplot <- DimPlot(seu, group.by = "gene_snn_res.0.2")

    patchwork <- genotype_distribution / (clone_dimplot + cluster_dimplot)

}

normal_reference_path <- "~/Homo_sapiens/infercnv/reference_mat.rds"

normal_reference_mat <- readRDS(normal_reference_path)
normal_seu <- Seurat::CreateSeuratObject(normal_reference_mat) %>%
  RenameCells(new.names = paste0("normal_", str_replace(colnames(.), "-", ".")))

cellranger_paths <-
  fs::dir_ls("output/cellranger/", glob = "*SRR*") %>%
  purrr::set_names(str_extract(path_file(.), "SRR[0-9]*"))

infercnv_seu_paths <-
  fs::dir_ls("output/seurat/", glob = "*SRR*cnv_seu.rds") %>%
  purrr::set_names(str_extract(path_file(.), "SRR[0-9]*"))

numbat_seu_paths <-
    fs::dir_ls("output/seurat/", glob = "*SRR*numbat_seu.rds") %>%
    purrr::set_names(str_extract(path_file(.), "SRR[0-9]*"))

seus <- purrr::map(numbat_seu_paths, readRDS)

cnv_cols <- c("proportion_dupli_1", "proportion_dupli_2", "proportion_dupli_6", "proportion_loss_16")

read_image_as_plot <- function(image_path){
  infercnv_image <- magick::image_read(image_path)

  infercnv_image <- ggplot() +
    ggpubr::background_image(infercnv_image) +
    # coord_fixed() +
    NULL

  return(infercnv_image)
}

infercnv_image_paths <- dir_ls("output/infercnv/", glob = "*infercnv.png", recurse = TRUE)

# save infercnv images ------------------------------
dir_create("results/infercnv")

new_infercnv_image_paths <- path("results/infercnv", paste0(path_file(path_dir(infercnv_image_paths)), "_infercnv.png"))

map2(infercnv_image_paths, new_infercnv_image_paths, fs::file_copy)

infercnv_images <- dir_ls("results/infercnv", glob = "*.png") %>%
  purrr::map(read_image_as_plot) %>%
  identity()

# save numbat images ------------------------------
numbat_dir <- "results/numbat"
dir_create(numbat_dir)

numbat_image_paths <- dir_ls("output/numbat", regexp = "\\/SRR[0-9]*\\/panel_2.png", recurse = TRUE)

new_numbat_image_paths <- path(numbat_dir, paste0(path_file(path_dir(numbat_image_paths)), "_numbat.png"))

numbat_image_paths <- dir_ls("output/numbat", regexp = "\\/SRR[0-9]*\\/bulk_clones_final.png", recurse = TRUE)

new_numbat_image_paths <- path(numbat_dir, paste0(path_file(path_dir(numbat_image_paths)), "_bulk_clones_numbat.png"))



map2(numbat_image_paths, new_numbat_image_paths, fs::file_copy)

numbat_images <- dir_ls("results/infercnv", glob = "*.png") %>%
  purrr::map(read_image_as_plot) %>%
  identity()

make_phylo_heatmap <- function(numbat_dir){
    nb = Numbat$new(out_dir = numbat_dir)

    nb$plot_consensus()

    mypal = scales::hue_pal()(n_distinct(nb$clone_post$clone_opt))

    names(mypal) <- seq(1:length(mypal))

    numbat_phylo_heatmap <- nb$plot_phylo_heatmap(
        clone_bar = TRUE,
        p_min = 0.5,
        pal_clone = mypal
    )

    return(numbat_phylo_heatmap)
}

numbat_dirs = path("output/numbat", str_extract(path_file(infercnv_seu_paths), "SRR[A-Z]*[0-9]*"))

possible_phylo <- purrr::possibly(make_phylo_heatmap, NA)

phylo_heatmaps <- map(numbat_dirs, possible_phylo)

numbat_phylo_images <- dir_ls("output/numbat", glob = "*panel_2.png", recurse = TRUE) %>%
    purrr::map(read_image_as_plot) %>%
    identity()

infercnv_plots <- purrr::map(seus, FeaturePlot, features = cnv_cols)

make_numbat_cnv_plot <-
    function(seu_path){
        seu <- readRDS(seu_path)
        DimPlot(seu, group.by = "clone_opt") +
            labs(title = seu_path)
    }
≥
numbat_cnv_plots <- purrr::map(numbat_seu_paths, make_numbat_cnv_plot)

numbat_clone_distributions <- purrr::map(numbat_seu_paths, plot_numbat_genotype_distribution)


pdf("results/numbat_clone_patchworks2.pdf")
numbat_clone_distributions
dev.off()

test0 <- plot_numbat_genotype_distribution(numbat_seu_paths[[1]])


pdf("results/numbat_cnv_plots.pdf")
numbat_cnv_plots0
dev.off()


marker_plots <- purrr::map(seus, ~plot_markers(.x, metavar = "gene_snn_res.0.2", marker_method = "presto", return_plotly = FALSE))

umap_plots <- purrr::imap(seus, ~(DimPlot(.x, group.by = "gene_snn_res.0.2") + labs(title = .y)))

library(patchwork)

assemble_cnv_patchwork <- function(umap, markers, cnv_image, cnv_plot){
  mypatchwork = (umap + markers + cnv_image + cnv_plot) + plot_layout(ncol = 2, widths = c(2,2))

  return(mypatchwork)
}

patchworks <- purrr::pmap(list(umap_plots, marker_plots, infercnv_images, infercnv_plots), assemble_cnv_patchwork)


pdf("results/patchworks2.pdf", width = 14, height = 14)
patchworks
dev.off()
