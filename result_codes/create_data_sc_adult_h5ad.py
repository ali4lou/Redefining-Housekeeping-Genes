
import scanpy as sc
import pandas as pd
import scipy.io
import gc
import os
import scanpy as sc
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt

os.chdir(r"PATH_WHERE_YOU_KEEP_YOUR_DATASETS")

mat = scipy.io.mmread("Zebrafish_counts.mtx")

mat = mat.T

mat = mat.tocsr()

adata = sc.AnnData(X=mat)

del mat
gc.collect()

adata.var_names = pd.read_csv("Zebrafish_genes.tsv", header=None)[0].values
adata.obs_names = pd.read_csv("Zebrafish_barcodes.tsv", header=None)[0].values

metadata = pd.read_csv("Zebrafish_metadata.csv", index_col=0)
adata.obs = metadata

adata.write("Zebrafish_Adult_FINAL.h5ad", compression="gzip")


