#!/usr/bin/env python
# coding: utf-8
# # Integrate (batch correct) single cell datasets with scanpy and scvi-tools

# %% [markdown]
# ## load required packages

# %%
import scanpy as sc
import scanpy.external as sce
from pathlib import Path
import anndata as ad
import re
import pandas as pd
import matplotlib.pyplot as plt
import os
import argparse
import grthub_tools as gt
import scvi


# %% [markdown]
# ## set scanpy default display settings 

# %%
sc.settings.verbosity = 1 # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(dpi=100, fontsize=10, dpi_save=300, figsize=(5,4), format='png')

# %% [markdown]
# ## load parse output directories for all samples
# https://support.parsebiosciences.com/hc/en-us/articles/360052794312-Scanpy-Tutorial-65k-PBMCs
# 

# %%

# Parse command line arguments

parser = argparse.ArgumentParser(description="Integrate single cell datasets with scanpy and scvi-tools", add_help=False)
parser.add_argument('--help', '-h', action='help', default=argparse.SUPPRESS,
                    help='Show this help message and exit')
parser.add_argument('--input_dir', type=str, required=True, help='Directory containing sample parse output folders', default='/dfs9/ucightf-lab/projects/LaSpA/250716_0725Bio-04_LaSpA_Parse-WT-100K/Trailmaker_results')
parser.add_argument('--output_prefix', type=str, default=None, help='Prefix for output files (default: input_dir name)')
parser.add_argument('--min_genes', type=int, default=200, help='Minimum genes per cell for filtering')
parser.add_argument('--min_cells', type=int, default=5, help='Minimum cells per gene for filtering')
parser.add_argument('--n_top_genes', type=int, default=2000, help='Number of highly variable genes to select')
parser.add_argument('--batch_key', type=str, default='batch', help='Batch key for integration')
parser.add_argument('--combine-only', action='store_true', default=False, help='Only combine and preprocess input datasets; skip Harmony integration')
args = parser.parse_args()

directory_path = Path(args.input_dir)
if args.output_prefix is None:
    output_prefix = str(directory_path.name)
else:
    output_prefix = args.output_prefix

# Look for directories containing DGE_filtered/count_matrix.mtx
data_paths = []
for root, dirs, files in os.walk(directory_path):
    if 'all-sample' not in root:
        # Check if this directory has DGE_filtered subdirectory with count_matrix.mtx
        dge_filtered_path = os.path.join(root, 'DGE_filtered')
        count_matrix_path = os.path.join(dge_filtered_path, 'count_matrix.mtx')
        if os.path.exists(count_matrix_path):
            data_paths.append(root)

print(data_paths)

# ## Reading in data
# 
# After reading in the data we'll perform basic filtering on our expression matrix to remove low-quality cells and uninformative genes. The parameter "min_genes" will keep cells that have at least 300 genes, and similarly, "min_cells" will keep genes that are expressed in at least 5 cells.

combined_adata_path = f"{output_prefix}.h5ad"
if os.path.exists(combined_adata_path):
    print(f"{combined_adata_path} already exists. Skipping data loading and concatenation.")
    combined_adata = ad.read_h5ad(combined_adata_path)
else:
    adatas = []
    for i, data_path in enumerate(data_paths):
        sample_id = os.path.basename(data_path)
        adata = gt.read_parse(data_path, min_genes=args.min_genes, min_cells=args.min_cells, n_top_genes=args.n_top_genes, batch_key=args.batch_key)
        adata.obs[args.batch_key] = sample_id
        adatas.append(adata)
    combined_adata = ad.concat(adatas, join='outer', uns_merge="first")
    combined_adata.var_names_make_unique()
    combined_adata.write(combined_adata_path)
    del adatas

print(output_prefix)
combined_adata = gt.preprocess_adata(combined_adata, n_top_genes=args.n_top_genes, batch = True)
if args.combine_only:
    # save the preprocessed combined AnnData and skip integration
    preprocessed_path = f"{output_prefix}.preprocessed.h5ad"
    combined_adata.write(preprocessed_path)
    print(f"Combine-only flag set; wrote preprocessed AnnData to {preprocessed_path} and skipped Harmony integration.")
else:
    gt.integrate_w_harmony(combined_adata, output_prefix)