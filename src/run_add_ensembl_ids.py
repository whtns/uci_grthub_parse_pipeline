#!/usr/bin/env python
"""Add Ensembl gene IDs to an h5ad file using multiple methods.

Also supports Cellarium CAS for cell type annotation.

Usage:
    python run_add_ensembl_ids.py [input_file] [--method METHOD] [--species SPECIES]
    python run_add_ensembl_ids.py [input_file] --use-cellarium [--cellarium-token TOKEN]

Example:
    # Fast gene ID mapping with MyGene (recommended)
    python run_add_ensembl_ids.py output/scanpy/combined_harmony_integrated.h5ad --method mygene
    
    # Slow gene ID mapping with Ensembl REST API
    python run_add_ensembl_ids.py output/scanpy/combined_harmony_integrated.h5ad --method ensembl-rest
    
    # Fastest gene ID mapping from local GTF file
    python src/run_add_ensembl_ids.py output/scanpy/combined_harmony_integrated.h5ad --method gtf --gtf-path /dfs8/commondata/cellranger/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf
    
    # Cellarium cell type annotation
    python run_add_ensembl_ids.py output/scanpy/combined_harmony_integrated.h5ad --use-cellarium
"""
import argparse
import os
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

# Disable Arrow-backed strings (pandas 2.x+) for HDF5 compatibility with anndata
pd.options.future.infer_string = False

import grthub_tools as gt

# MyGene imports (fast batch queries)
try:
    import mygene
    MYGENE_AVAILABLE = True
except ImportError:
    MYGENE_AVAILABLE = False

# Cellarium imports
try:
    from cellarium.cas.client import CASClient
    from cellarium.cas.postprocessing import insert_cas_ontology_aware_response_into_adata
    import cellarium.cas.postprocessing.ontology_aware as pp
    from cellarium.cas.postprocessing.cell_ontology import CellOntologyCache
    from cellarium.cas._io import suppress_stderr
    CELLARIUM_AVAILABLE = True
except ImportError:
    CELLARIUM_AVAILABLE = False
    CASClient = None


def convert_arrow_to_string(adata):
    """Convert ArrowStringArray to regular strings for HDF5 compatibility.
    
    This fixes the IORegistryError when writing h5ad files with newer pandas versions
    that use ArrowStringArray by default.
    """
    print("Converting string arrays for HDF5 compatibility...")
    
    # Force index to object dtype by rebuilding from list
    adata.obs.index = pd.Index(list(adata.obs.index), dtype='object', name=adata.obs.index.name)
    adata.var.index = pd.Index(list(adata.var.index), dtype='object', name=adata.var.index.name)
    
    # Convert all string columns to object dtype
    for col in adata.obs.columns:
        try:
            if pd.api.types.is_string_dtype(adata.obs[col]):
                # Convert via list to force object dtype
                adata.obs[col] = pd.Series(list(adata.obs[col]), index=adata.obs.index, dtype='object')
        except:
            pass
    
    for col in adata.var.columns:
        try:
            if pd.api.types.is_string_dtype(adata.var[col]):
                # Convert via list to force object dtype
                adata.var[col] = pd.Series(list(adata.var[col]), index=adata.var.index, dtype='object')
        except:
            pass
    
    # Remove raw to avoid any remaining issues
    if adata.raw is not None:
        print("Note: Removing .raw attribute for HDF5 compatibility")
        adata.raw = None


def main():
    parser = argparse.ArgumentParser(
        description="Add Ensembl gene IDs to an AnnData h5ad file or annotate with Cellarium CAS"
    )
    parser.add_argument(
        "input_file",
        type=str,
        help="Path to input h5ad file"
    )
    parser.add_argument(
        "--species",
        type=str,
        default="human",
        help="Species for Ensembl lookup (default: human)"
    )
    parser.add_argument(
        "--var-name",
        type=str,
        default="gene_ids",
        help="Name of the column to add to adata.var (default: gene_ids)"
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing gene_ids column if present"
    )
    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="Output file path (default: overwrites input file)"
    )
    parser.add_argument(
        "--method",
        type=str,
        default="ensembl-rest",
        choices=["ensembl-rest", "mygene", "gtf"],
        help="Method for gene ID mapping: ensembl-rest (slow), mygene (fast, recommended), gtf (fastest, needs file)"
    )
    parser.add_argument(
        "--gtf-path",
        type=str,
        default=None,
        help="Path to GTF annotation file (required if --method gtf)"
    )
    parser.add_argument(
        "--sleep",
        type=float,
        default=0.05,
        help="Sleep duration between API calls in seconds for ensembl-rest (default: 0.05)"
    )
    parser.add_argument(
        "--timeout",
        type=float,
        default=5.0,
        help="API request timeout in seconds for ensembl-rest (default: 5.0)"
    )
    parser.add_argument(
        "--use-cellarium",
        action="store_true",
        help="Use Cellarium CAS for cell type annotation"
    )
    parser.add_argument(
        "--cellarium-token",
        type=str,
        default=os.environ.get("CELLARIUM_API_KEY"),
        help="Cellarium API token (defaults to CELLARIUM_API_KEY environment variable)"
    )
    parser.add_argument(
        "--cluster-key",
        type=str,
        default="cluster_label",
        help="Column in adata.obs containing cluster labels (for Cellarium)"
    )
    parser.add_argument(
        "--min-score",
        type=float,
        default=0.1,
        help="Minimum acceptable relevance score for Cellarium cell type calls (default: 0.1)"
    )
    parser.add_argument(
        "--top-k",
        type=int,
        default=5,
        help="Number of top cell type calls to return from Cellarium (default: 5)"
    )
    parser.add_argument(
        "--pseudobulk",
        action="store_true",
        help="Use pseudobulk approach for Cellarium annotation (aggregate by cluster)"
    )
    
    args = parser.parse_args()
    
    # Validate input file exists
    input_path = Path(args.input_file)
    if not input_path.exists():
        print(f"Error: Input file not found: {args.input_file}", file=sys.stderr)
        sys.exit(1)
    
    print(f"Loading AnnData from: {args.input_file}")
    adata = ad.read_h5ad(args.input_file)
    
    print(f"AnnData shape: {adata.shape}")
    print(f"Number of genes: {adata.n_vars}")
    print(f"Gene index sample: {adata.var.index[:5].tolist()}")
    
    # Fix HDF5 compatibility issues with Arrow string arrays BEFORE any processing
    # This is needed because Cellarium serializes adata internally
    convert_arrow_to_string(adata)
    
    # Decide on workflow
    if args.use_cellarium:
        # Cellarium workflow
        if not CELLARIUM_AVAILABLE:
            print("ERROR: Cellarium is not available. Install with: pip install cellarium-cas")
            sys.exit(1)

        # SCT/Seurat objects often have X=None; restore from counts layer
        if adata.X is None:
            if "counts" in adata.layers:
                print("adata.X is None — restoring from 'counts' layer for Cellarium...")
                adata.X = adata.layers["counts"]
            else:
                print("ERROR: adata.X is None and no 'counts' layer found.", file=sys.stderr)
                sys.exit(1)

        gt.annotate_with_cellarium(
            adata,
            api_token=args.cellarium_token,
            cluster_key=args.cluster_key,
            min_score=args.min_score,
            top_k=args.top_k,
            pseudobulk=args.pseudobulk
        )
        
        # Show sample results
        cas_cols = [c for c in adata.obs.columns if c.startswith('cas_') or c.startswith('pseudobulk_cas_')]
        if cas_cols:
            print(f"\nSample Cellarium annotations:")
            print(adata.obs[cas_cols].head(10))
    
    else:
        # Gene ID mapping workflow (choose method)
        # Check if gene_ids already exists
        if args.var_name in adata.var.columns and not args.overwrite:
            print(f"Column '{args.var_name}' already exists in adata.var")
            print("Use --overwrite to replace existing values")
            sys.exit(0)
        
        if args.method == "mygene":
            # Fast batch method using MyGene (RECOMMENDED)
            if not MYGENE_AVAILABLE:
                print("ERROR: mygene not available. Install with: pip install mygene")
                print("Falling back to ensembl-rest method...")
                args.method = "ensembl-rest"
            else:
                gt.add_ensembl_ids_mygene(
                    adata,
                    species=args.species,
                    var_name=args.var_name,
                    overwrite=args.overwrite
                )
        
        if args.method == "gtf":
            # Fastest method using local GTF file
            if not args.gtf_path:
                print("ERROR: --gtf-path required when using --method gtf")
                sys.exit(1)
            gt.add_ensembl_ids_from_gtf(
                adata,
                gtf_path=args.gtf_path,
                var_name=args.var_name,
                overwrite=args.overwrite
            )
        
        if args.method == "ensembl-rest":
            # Slowest method: one-at-a-time REST API queries
            print(f"\nMapping gene symbols to Ensembl IDs (species: {args.species})...")
            print(f"This may take a while for {adata.n_vars} genes...")
            print(f"Using {args.sleep}s sleep between API calls")
            print(f"TIP: Use --method mygene for 50-100x faster mapping!")
            
            gt.add_ensembl_ids(
                adata,
                species=args.species,
                var_name=args.var_name,
                overwrite=args.overwrite,
                sleep=args.sleep,
                timeout=args.timeout
            )
        
        # Show results
        n_mapped = adata.var[args.var_name].notna().sum()
        print(f"\nMapping complete!")
        print(f"Successfully mapped: {n_mapped}/{adata.n_vars} genes ({100*n_mapped/adata.n_vars:.1f}%)")
        print(f"Sample mappings:")
        sample_df = adata.var[[args.var_name]].head(10)
        print(sample_df)
    
    # Save the modified AnnData
    if args.output:
        output_path = args.output
    else:
        input_suffix = input_path.suffix if input_path.suffix else ".h5ad"
        if args.use_cellarium:
            output_path = str(input_path.with_name(f"{input_path.stem}_cellarium{input_suffix}"))
        else:
            output_path = str(input_path.with_name(f"{input_path.stem}_ensembl_ids{input_suffix}"))
    print(f"\nSaving modified AnnData to: {output_path}")
    
    adata.write_h5ad(output_path)
    print("Done!")


if __name__ == "__main__":
    main()
