#!/usr/bin/env python3

import argparse
from pathlib import Path
import scanpy as sc
import sys


def main():
	parser = argparse.ArgumentParser(description="Read an AnnData file and save it (works with Snakemake fallback)")
	parser.add_argument("input", nargs="?", help="Path to input .h5ad file")
	parser.add_argument("data_output", nargs="?", help="Path to output .h5ad file")
	parser.add_argument("scaled_data_output", nargs="?", help="Path to output .h5ad file")
	args = parser.parse_args()

	in_path = args.input
	out_path = args.data_output

	# If running under Snakemake, allow using snakemake.input/output when args omitted
	if (in_path is None) or (out_path is None):
		try:
			# 'snakemake' is injected when executing via Snakemake
			snakemake  # type: ignore
			in_path = in_path or snakemake.input[0]
			out_path = out_path or snakemake.output[0]
		except NameError:
			pass

	if in_path is None or out_path is None:
		parser.error("input and output paths required (or run via Snakemake with input/output)")

	in_path = Path(in_path)
	out_path = Path(out_path)
	if not in_path.exists():
		print(f"ERROR: input file not found: {in_path}", file=sys.stderr)
		sys.exit(2)
	adata = sc.read_h5ad(str(in_path))
	print(f"Wrote AnnData to: {args.data_output}")
	sc.pp.scale(adata, max_value=10)
	adata.write_h5ad(str(args.scaled_data_output))


if __name__ == "__main__":
	main()
