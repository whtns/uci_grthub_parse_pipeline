import scanpy as sc
import pandas as pd

adata = sc.read_h5ad("../output/scanpy/combined.h5ad")
meta = pd.read_csv("data/metadata.csv", index_col=0)

adata.obs.merge(metadata, left_on = "sample", right_on = "sample", how = "left")