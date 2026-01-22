#!/usr/bin/env python
"""pull_scaled_data.py

Usage:
  python src/pull_scaled_data.py /path/to/adata.h5ad [--out-dir OUTDIR] [--prefix PREFIX]

Reads an AnnData `.h5ad` file and writes a matrix market file of the scaled
data (genes x cells) plus CSVs for cell barcodes and gene names.
"""

import argparse
from pathlib import Path
import pandas as pd
import scanpy as sc
import sys

def main():
	parser = argparse.ArgumentParser(description="Export scaled data from AnnData to MTX and CSV files")
	parser.add_argument("adata_path", help="Path to input AnnData .h5ad file")
	parser.add_argument("--out-dir", default=None, help="Directory to write outputs (defaults to input file parent)")
	parser.add_argument("--prefix", default="scaled_data", help="Prefix for output files")
	args = parser.parse_args()

	adata_path = Path(args.adata_path)
	if not adata_path.exists():
		print(f"ERROR: input file not found: {adata_path}", file=sys.stderr)
		sys.exit(2)

	out_dir = Path(args.out_dir) if args.out_dir else adata_path.parent
	out_dir.mkdir(parents=True, exist_ok=True)
	prefix = args.prefix

	print(f"Reading AnnData from: {adata_path}")
	adata = sc.read_h5ad(str(adata_path))

	# Extract scaled data matrix. Ensure orientation genes x cells for Seurat.
	scaled = adata.X
	# AnnData X is typically cells x genes; we want genes x cells for Seurat
	# Handle sparse, HDF-backed, and dense arrays robustly.
	import numpy as np
	from scipy import sparse
	if sparse.issparse(scaled):
		scaled_t = scaled.transpose()
	else:
		# If object provides transpose/transpose method, use it; otherwise convert to numpy
		if hasattr(scaled, "transpose"):
			try:
				scaled_t = scaled.transpose()
			except Exception:
				scaled_t = np.asarray(scaled).T
		else:
			scaled_t = np.asarray(scaled).T

	cell_barcodes = adata.obs_names.to_numpy()
	gene_names = adata.var_names.to_numpy()

	# Write outputs
	from scipy import io as sio
	mtx_path = out_dir / f"{prefix}.mtx"
	print(f"Writing matrix market to: {mtx_path}")
	# If dense, convert to sparse for mmwrite
	try:
		sio.mmwrite(str(mtx_path), scaled_t)
	except Exception:
		from scipy import sparse
		sparse_mat = sparse.csr_matrix(scaled_t)
		sio.mmwrite(str(mtx_path), sparse_mat)

	barcodes_path = out_dir / f"{prefix}_barcodes.csv"
	genes_path = out_dir / f"{prefix}_genes.csv"
	pd.DataFrame(cell_barcodes, columns=["Barcode"]).to_csv(barcodes_path, index=False)
	pd.DataFrame(gene_names, columns=["GeneName"]).to_csv(genes_path, index=False)

	print("Done.")

if __name__ == '__main__':
	main()
