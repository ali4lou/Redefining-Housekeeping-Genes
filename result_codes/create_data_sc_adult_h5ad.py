#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 11:08:23 2026

@author: alicia
"""

import scanpy as sc
import pandas as pd
import scipy.io
import gc
import os
import scanpy as sc
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt

os.chdir(r"/media/alicia/TOSHIBA EXT/zebrafish/revision")

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


