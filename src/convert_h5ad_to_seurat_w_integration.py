#!/usr/bin/env python3
"""
Convert an AnnData .h5ad file to a Seurat object (.rds), preserving dimensional reductions when possible.

Primary method: call an Rscript that uses SeuratDisk::Convert + SeuratDisk::LoadH5Seurat then saveRDS.
Fallback: use rpy2 to call Seurat functions directly from Python (may require converting matrices to dense).

Usage:
  python src/convert_h5ad_to_seurat_w_integration.py \
    --input /path/to/combined_harmony_integrated.h5ad \
    --output /path/to/combined_harmony_integrated.seurat.rds

Notes:
- Requires R + Seurat + SeuratDisk for primary pathway.
- For fallback, requires Python packages: anndata, rpy2, numpy, pandas, scipy (optional).
"""

import argparse
import os
import sys
import tempfile
import subprocess
import textwrap
import shutil


def run_rscript_conversion(infile, outfile, h5seurat=None):
    """Write a temporary R script that uses SeuratDisk to convert and load, then run it with Rscript."""
    rscript = textwrap.dedent("""
    args <- commandArgs(trailingOnly=TRUE)
    infile <- args[1]
    outfile <- args[2]
    h5seurat <- if (length(args) >= 3 && nchar(args[3])>0) args[3] else sub('\\.h5ad$', '.h5seurat', infile)

    if(!requireNamespace('SeuratDisk', quietly=TRUE)){
      stop('SeuratDisk is required in R for conversion. Install with: remotes::install_github("mojaveazure/seurat-disk") or Bioconductor/CRAN if available')
    }
    library(SeuratDisk)

    message(sprintf('Converting %s -> %s', infile, h5seurat))
    Convert(infile, dest = 'h5seurat', overwrite = TRUE)
    message('Loading h5seurat...')
    seu <- LoadH5Seurat(h5seurat)
    message(sprintf('Saving Seurat RDS to %s', outfile))
    saveRDS(seu, file = outfile)
    message('done')
    """)

    with tempfile.TemporaryDirectory() as td:
        rpath = os.path.join(td, 'convert_and_save.R')
        with open(rpath, 'w') as fh:
            fh.write(rscript)
        cmd = ['Rscript', rpath, infile, outfile]
        if h5seurat:
            cmd.append(h5seurat)
        try:
            print('Running Rscript conversion:', ' '.join(cmd))
            subprocess.run(cmd, check=True)
            return True
        except subprocess.CalledProcessError as e:
            print('Rscript conversion failed:', e, file=sys.stderr)
            return False


def fallback_rpy2_build(infile, outfile):
    """Fallback: use rpy2 to build a Seurat object from Python/AnnData objects.

    This approach may be memory-heavy because it often converts sparse matrices to dense.
    """
    try:
        import anndata
        import numpy as np
        import pandas as pd
    except Exception as e:
        print('Fallback requires anndata, numpy, pandas: ', e, file=sys.stderr)
        return False

    try:
        from rpy2 import robjects
        from rpy2.robjects import r, pandas2ri, numpy2ri
        from rpy2.robjects.packages import importr
        pandas2ri.activate()
        numpy2ri.activate()
    except Exception as e:
        print('rpy2 is required for fallback but is not available:', e, file=sys.stderr)
        return False

    print('Reading AnnData...')
    ad = anndata.read_h5ad(infile)

    # Extract main data matrix
    X = ad.X
    is_sparse = False
    try:
        import scipy.sparse as sp
        is_sparse = sp.issparse(X)
    except Exception:
        is_sparse = False

    if is_sparse:
        print('Converting sparse matrix to dense for rpy2 transfer (may require large memory).')
        X = X.toarray()

    # AnnData X is usually cells x features; Seurat expects features x cells counts matrix.
    # We'll transpose before sending to R.
    X_t = X.T

    # Convert obs to pandas DataFrame
    obs = ad.obs.copy()

    # Launch R and create Seurat object
    try:
        base = importr('base')
        seurat = importr('Seurat')
    except Exception as e:
        print('R package Seurat must be installed for fallback pathway:', e, file=sys.stderr)
        return False

    print('Converting data to R and building Seurat object...')
    # Convert counts to R matrix
    r_matrix = robjects.r['matrix'](numpy2ri.py2rpy(X_t), nrow=X_t.shape[0], ncol=X_t.shape[1])
    # Set dimnames (features x cells)
    try:
        feat_names = list(ad.var_names.astype(str))
        cell_names = list(ad.obs_names.astype(str))
        dimnames = robjects.ListVector({'rownames': robjects.StrVector(feat_names), 'colnames': robjects.StrVector(cell_names)})
        r('dimnames')(r_matrix)[0]
        # set dimnames via R-level assignment
        robjects.r['dimnames'](r_matrix)
    except Exception:
        pass

    # Create Seurat object in R: CreateSeuratObject(counts = r_matrix)
    try:
        create_call = robjects.r['CreateSeuratObject']
        seu = create_call(robjects.r['as.matrix'](r_matrix))
    except Exception as e:
        print('Failed to create Seurat object via rpy2:', e, file=sys.stderr)
        return False

    # Add metadata
    try:
        meta_r = pandas2ri.py2rpy(obs)
        robjects.r['rownames'](meta_r)
        seu = robjects.r['AddMetaData'](seu, metadata=meta_r)
    except Exception:
        pass

    # Add reductions from ad.obsm where possible
    try:
        if hasattr(ad, 'obsm') and len(ad.obsm.keys()) > 0:
            for key in ad.obsm.keys():
                mat = ad.obsm[key]
                if is_sparse:
                    # if obsm sparse (rare), convert
                    try:
                        import scipy.sparse as sp
                        if sp.issparse(mat):
                            mat = mat.toarray()
                    except Exception:
                        pass
                # Ensure shape is cells x dims
                if mat.shape[0] == len(ad.obs_names):
                    emb = mat
                elif mat.shape[1] == len(ad.obs_names):
                    emb = mat.T
                else:
                    print(f'Skipping reduction {key}: incompatible shape {mat.shape}', file=sys.stderr)
                    continue
                # Convert to R matrix (cells x dims)
                r_emb = numpy2ri.py2rpy(emb)
                # Create DimReduc object
                try:
                    # CreateDimReducObject(emb = t(r_emb), key = 'UMAP_', assay = DefaultAssay(seu))
                    robjects.r('library(Seurat)')
                    CreateDimReduc = robjects.r['CreateDimReducObject']
                    DefaultAssay = robjects.r['DefaultAssay']
                    assay_name = DefaultAssay(seu)
                    # emb needs to be a matrix with cells in rows
                    dr = CreateDimReduc(emb=r_emb, key=key.upper()[:4]+'_', assay=assay_name)
                    # assign to reductions slot
                    robjects.r['reductions<-'](seu, robjects.ListVector({key: dr}))
                    print(f'Added reduction {key}')
                except Exception as e:
                    print('Failed to add reduction', key, e, file=sys.stderr)
    except Exception:
        pass

    # Save to RDS
    try:
        robjects.r['saveRDS'](seu, file=outfile)
        print('Saved Seurat RDS to', outfile)
        return True
    except Exception as e:
        print('Failed to save RDS via rpy2:', e, file=sys.stderr)
        return False


def main():
    p = argparse.ArgumentParser(description='Convert h5ad -> Seurat RDS (preserve reductions when possible)')
    p.add_argument('--input', '-i', required=False, default='/dfs9/ucightf-lab/projects/OlabR/250922_0925Bio-10_OlabR_parse/output/scanpy/combined_harmony_integrated.h5ad')
    p.add_argument('--output', '-o', required=False, default=None)
    p.add_argument('--h5seurat', default=None, help='Optional path for intermediate .h5seurat')
    args = p.parse_args()

    infile = args.input
    if args.output:
        outfile = args.output
    else:
        base = os.path.splitext(os.path.basename(infile))[0]
        outfile = os.path.join(os.path.dirname(infile), base + '.seurat.rds')

    if not os.path.exists(infile):
        print('Input file not found:', infile, file=sys.stderr)
        sys.exit(2)

    # Prefer Rscript + SeuratDisk conversion
    rscript_path = shutil.which('Rscript')
    if rscript_path:
        ok = run_rscript_conversion(infile, outfile, args.h5seurat)
        if ok:
            print('Conversion succeeded via Rscript/SeuratDisk. Output:', outfile)
            sys.exit(0)
        else:
            print('Rscript/SeuratDisk conversion failed, attempting fallback...', file=sys.stderr)
    else:
        print('Rscript not found on PATH; attempting fallback via rpy2', file=sys.stderr)

    ok2 = fallback_rpy2_build(infile, outfile)
    if ok2:
        print('Conversion succeeded via rpy2 fallback. Output:', outfile)
        sys.exit(0)
    else:
        print('All conversion attempts failed. Please ensure R + Seurat + SeuratDisk are installed, or install rpy2 + anndata in Python.', file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
#!/usr/bin/env python3
"""
Convert an AnnData .h5ad file to a Seurat object (.rds), preserving dimensional reductions when possible.

Primary method: call an Rscript that uses SeuratDisk::Convert + SeuratDisk::LoadH5Seurat then saveRDS.
Fallback: use rpy2 to call Seurat functions directly from Python (may require converting matrices to dense).

Usage:
  python src/convert_h5ad_to_seurat_w_integration.py \
    --input /path/to/combined_harmony_integrated.h5ad \
    --output /path/to/combined_harmony_integrated.seurat.rds

Notes:
"""

import argparse
import os
import sys
import tempfile
import subprocess
import textwrap
import shutil


def run_rscript_conversion(infile, outfile, h5seurat=None):
    """Write a temporary R script that uses SeuratDisk to convert and load, then run it with Rscript."""
    rscript = textwrap.dedent("""
    args <- commandArgs(trailingOnly=TRUE)
    infile <- args[1]
    outfile <- args[2]
    h5seurat <- if (length(args) >= 3 && nchar(args[3])>0) args[3] else sub('\\.h5ad$', '.h5seurat', infile)

    if(!requireNamespace('SeuratDisk', quietly=TRUE)){
      stop('SeuratDisk is required in R for conversion. Install with: remotes::install_github("mojaveazure/seurat-disk") or Bioconductor/CRAN if available')
    }
    library(SeuratDisk)

    message(sprintf('Converting %s -> %s', infile, h5seurat))
    Convert(infile, dest = 'h5seurat', overwrite = TRUE)
    message('Loading h5seurat...')
    seu <- LoadH5Seurat(h5seurat)
    message(sprintf('Saving Seurat RDS to %s', outfile))
    saveRDS(seu, file = outfile)
    message('done')
    """)

    with tempfile.TemporaryDirectory() as td:
        rpath = os.path.join(td, 'convert_and_save.R')
        with open(rpath, 'w') as fh:
            fh.write(rscript)
        cmd = ['Rscript', rpath, infile, outfile]
        if h5seurat:
            cmd.append(h5seurat)
        try:
            print('Running Rscript conversion:', ' '.join(cmd))
            subprocess.run(cmd, check=True)
            return True
        except subprocess.CalledProcessError as e:
            print('Rscript conversion failed:', e, file=sys.stderr)
            return False


def fallback_rpy2_build(infile, outfile):
    """Fallback: use rpy2 to build a Seurat object from Python/AnnData objects.

    This approach may be memory-heavy because it often converts sparse matrices to dense.
    """
    try:
        import anndata
        import numpy as np
        import pandas as pd
    except Exception as e:
        print('Fallback requires anndata, numpy, pandas: ', e, file=sys.stderr)
        return False

    try:
        from rpy2 import robjects
        from rpy2.robjects import r, pandas2ri, numpy2ri
        from rpy2.robjects.packages import importr
        pandas2ri.activate()
        numpy2ri.activate()
    except Exception as e:
        print('rpy2 is required for fallback but is not available:', e, file=sys.stderr)
        return False

    print('Reading AnnData...')
    ad = anndata.read_h5ad(infile)

    # Extract main data matrix
    X = ad.X
    is_sparse = False
    try:
        import scipy.sparse as sp
        is_sparse = sp.issparse(X)
    except Exception:
        is_sparse = False

    if is_sparse:
        print('Converting sparse matrix to dense for rpy2 transfer (may require large memory).')
        X = X.toarray()

    # AnnData X is usually cells x features; Seurat expects features x cells counts matrix.
    # We'll transpose before sending to R.
    X_t = X.T

    # Convert obs to pandas DataFrame
    obs = ad.obs.copy()

    # Launch R and create Seurat object
    try:
        base = importr('base')
        seurat = importr('Seurat')
    except Exception as e:
        print('R package Seurat must be installed for fallback pathway:', e, file=sys.stderr)
        return False

    print('Converting data to R and building Seurat object...')
    # Convert counts to R matrix
    r_matrix = robjects.r['matrix'](numpy2ri.py2rpy(X_t), nrow=X_t.shape[0], ncol=X_t.shape[1])
    # Set dimnames (features x cells)
    try:
        feat_names = list(ad.var_names.astype(str))
        cell_names = list(ad.obs_names.astype(str))
        dimnames = robjects.ListVector({'rownames': robjects.StrVector(feat_names), 'colnames': robjects.StrVector(cell_names)})
        r('dimnames')(r_matrix)[0]
        # set dimnames via R-level assignment
        robjects.r['dimnames'](r_matrix)
    except Exception:
        pass

    # Create Seurat object in R: CreateSeuratObject(counts = r_matrix)
    try:
        create_call = robjects.r['CreateSeuratObject']
        seu = create_call(robjects.r['as.matrix'](r_matrix))
    except Exception as e:
        print('Failed to create Seurat object via rpy2:', e, file=sys.stderr)
        return False

    # Add metadata
    try:
        meta_r = pandas2ri.py2rpy(obs)
        robjects.r['rownames'](meta_r)
        seu = robjects.r['AddMetaData'](seu, metadata=meta_r)
    except Exception:
        pass

    # Add reductions from ad.obsm where possible
    try:
        if hasattr(ad, 'obsm') and len(ad.obsm.keys()) > 0:
            for key in ad.obsm.keys():
                mat = ad.obsm[key]
                if is_sparse:
                    # if obsm sparse (rare), convert
                    try:
                        import scipy.sparse as sp
                        if sp.issparse(mat):
                            mat = mat.toarray()
                    except Exception:
                        pass
                # Ensure shape is cells x dims
                if mat.shape[0] == len(ad.obs_names):
                    emb = mat
                elif mat.shape[1] == len(ad.obs_names):
                    emb = mat.T
                else:
                    print(f'Skipping reduction {key}: incompatible shape {mat.shape}', file=sys.stderr)
                    continue
                # Convert to R matrix (cells x dims)
                r_emb = numpy2ri.py2rpy(emb)
                # Create DimReduc object
                try:
                    # CreateDimReducObject(emb = t(r_emb), key = 'UMAP_', assay = DefaultAssay(seu))
                    robjects.r('library(Seurat)')
                    CreateDimReduc = robjects.r['CreateDimReducObject']
                    DefaultAssay = robjects.r['DefaultAssay']
                    assay_name = DefaultAssay(seu)
                    # emb needs to be a matrix with cells in rows
                    dr = CreateDimReduc(emb=r_emb, key=key.upper()[:4]+'_', assay=assay_name)
                    # assign to reductions slot
                    robjects.r['reductions<-'](seu, robjects.ListVector({key: dr}))
                    print(f'Added reduction {key}')
                except Exception as e:
                    print('Failed to add reduction', key, e, file=sys.stderr)
    except Exception:
        pass

    # Save to RDS
    try:
        robjects.r['saveRDS'](seu, file=outfile)
        print('Saved Seurat RDS to', outfile)
        return True
    except Exception as e:
        print('Failed to save RDS via rpy2:', e, file=sys.stderr)
        return False


def main():
    p = argparse.ArgumentParser(description='Convert h5ad -> Seurat RDS (preserve reductions when possible)')
    p.add_argument('--input', '-i', required=False, default='/dfs9/ucightf-lab/projects/OlabR/250922_0925Bio-10_OlabR_parse/output/scanpy/combined_harmony_integrated.h5ad')
    p.add_argument('--output', '-o', required=False, default=None)
    p.add_argument('--h5seurat', default=None, help='Optional path for intermediate .h5seurat')
    args = p.parse_args()

    infile = args.input
    if args.output:
        outfile = args.output
    else:
        base = os.path.splitext(os.path.basename(infile))[0]
        outfile = os.path.join(os.path.dirname(infile), base + '.seurat.rds')

    if not os.path.exists(infile):
        print('Input file not found:', infile, file=sys.stderr)
        sys.exit(2)

    # Prefer Rscript + SeuratDisk conversion
    rscript_path = shutil.which('Rscript')
    if rscript_path:
        ok = run_rscript_conversion(infile, outfile, args.h5seurat)
        if ok:
            print('Conversion succeeded via Rscript/SeuratDisk. Output:', outfile)
            sys.exit(0)
        else:
            print('Rscript/SeuratDisk conversion failed, attempting fallback...', file=sys.stderr)
    else:
        print('Rscript not found on PATH; attempting fallback via rpy2', file=sys.stderr)

    ok2 = fallback_rpy2_build(infile, outfile)
    if ok2:
        print('Conversion succeeded via rpy2 fallback. Output:', outfile)
        sys.exit(0)
    else:
        print('All conversion attempts failed. Please ensure R + Seurat + SeuratDisk are installed, or install rpy2 + anndata in Python.', file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()