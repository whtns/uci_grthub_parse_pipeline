#!/usr/bin/env python3

import argparse
import scrublet as scr
import scipy.io
import h5py
import numpy as np
import pandas as pd
import os

def load_matrix_h5(matrix_h5_path):
    with h5py.File(matrix_h5_path, 'r') as f:
        # 10x h5 format: /matrix/ data, indices, indptr, shape
        data = f['matrix']['data'][:]
        indices = f['matrix']['indices'][:]
        indptr = f['matrix']['indptr'][:]
        shape = f['matrix']['shape'][:]
        from scipy.sparse import csc_matrix
        matrix = csc_matrix((data, indices, indptr), shape=shape)
        return matrix.transpose()  # cells x genes

def main():
    parser = argparse.ArgumentParser(description="Run Scrublet doublet detection on 10x filtered_feature_bc_matrix.h5")
    parser.add_argument('--input', required=True, help='Path to filtered_feature_bc_matrix.h5')
    parser.add_argument('--output', required=True, help='Output CSV file for doublet scores and calls')
    parser.add_argument('--expected_doublet_rate', type=float, default=0.06, help='Expected doublet rate (default: 0.06)')
    args = parser.parse_args()

    print(f"Loading matrix from {args.input}")
    counts_matrix = load_matrix_h5(args.input)

    print("Running Scrublet...")
    scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=args.expected_doublet_rate)
    doublet_scores, predicted_doublets = scrub.scrub_doublets()

    # Get barcodes
    with h5py.File(args.input, 'r') as f:
        barcodes = [b.decode() for b in f['matrix']['barcodes'][:]]

    df = pd.DataFrame({
        'barcode': barcodes,
        'doublet_score': doublet_scores,
        'predicted_doublet': predicted_doublets
    })
    df.to_csv(args.output, index=False)
    print(f"Scrublet results saved to {args.output}")

if __name__ == "__main__":
    main()