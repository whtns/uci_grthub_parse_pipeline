#!/usr/bin/env python3
import argparse
import os
from pathlib import Path

import scanpy as sc
import pandas as pd


def main():
	parser = argparse.ArgumentParser(
		description="Append metadata to an AnnData object's obs and write a new h5ad file."
	)
	parser.add_argument(
		"--adatapath",
		help="Path to input AnnData h5ad file",
		default="../output/scanpy/combined.h5ad",
	)
	parser.add_argument(
		"--metadatapath",
		help="Path to metadata CSV file",
		default="data/metadata.csv",
	)

	args = parser.parse_args()

	adata_path = Path(args.adatapath)
	if not adata_path.exists():
		raise FileNotFoundError(f"AnnData file not found: {adata_path}")

	meta_path = Path(args.metadatapath)
	if not meta_path.exists():
		raise FileNotFoundError(f"Metadata CSV file not found: {meta_path}")

	adata = sc.read_h5ad(str(adata_path))
	metadata = pd.read_csv(str(meta_path), index_col=0)

	# Merge metadata into adata.obs. Keep original index unless merge changes it.
	merged = adata.obs.merge(metadata, left_on="sample", right_on="sample", how="left")
	# Assign merged dataframe back to adata.obs
	adata.obs = merged

	adata.write(adata_path)


if __name__ == "__main__":
	main()
