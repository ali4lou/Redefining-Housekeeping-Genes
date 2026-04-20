# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 18:39:16 2026

@author: logslab
"""

import scanpy as sc
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt

# ==========================================
# PATHS
# ==========================================
path_save_data = 'PATH_WHERE_YOU_KEEP_YOUR_DATASETS'
path_genes = 'PATH_WHERE_YOU_KEEP_YOUR_DATASETS'
path_dictionary = 'PATH_WHERE_YOU_KEEP_YOUR_DATASETS' 
h5ad_path = 'PATH_WHERE_YOU_KEEP_YOUR_DATASETS'
path_metadata = 'PATH_WHERE_YOU_KEEP_YOUR_DATASETS'


#TISSUE

# ==========================================
# 1.) TRANSLATE GENES
# ==========================================
translator = {}
with open(path_dictionary, 'r') as f:
    for line in f:
        columns = line.strip().split() 
        if len(columns) >= 4:
            symbol_lower = columns[2].lower() 
            ensembl_code = columns[3]            
            translator[ensembl_code] = symbol_lower 
            translator[columns[2].lower()] = symbol_lower    

with open(path_genes, 'r') as f:
    raw_genes = [line.strip() for line in f if line.strip()]

U_tissue_genes = []
for gene in raw_genes:
    if gene in translator:
        U_tissue_genes.append(translator[gene])
    else:
        clean_gene = gene.split('.')[0].lower()
        U_tissue_genes.append(clean_gene)
U_tissue_genes = np.unique(U_tissue_genes)

# ==========================================
# 2.) READ DATA (LIGHTWEIGHT MODE)
# ==========================================
adata = sc.read_h5ad(h5ad_path, backed='r') 

metadata = pd.read_csv(path_metadata, index_col='barcodes')
obs_names = adata.obs_names.to_series()
aligned_meta = metadata.reindex(obs_names)

var_names_clean = adata.var_names.to_series().str.split('.').str[0].str.lower()
mask_genes = var_names_clean.isin(U_tissue_genes)
pos_genes = np.where(mask_genes)[0]

# ==========================================
# 3.) STRICT PHYSICAL FILTERING
# ==========================================
print("--- STARTING STRICT FILTERING ---")
target_stages = ['22Month']

# 3.1 Identify exact cells belonging to these 3 stages
mask_stages = aligned_meta['stage'].isin(target_stages)
pos_cells = np.where(mask_stages)[0]

# 3.2 Extract ONLY those cells and ONLY the U genes to RAM
# This creates an adata_sub object completely free of embryonic cells
print("Extracting the subset into memory...")
adata_sub = adata[pos_cells, pos_genes].to_memory()

# Attach the corresponding metadata and clean gene names
adata_sub.obs = aligned_meta.iloc[pos_cells].copy()
adata_sub.var_names = var_names_clean.iloc[pos_genes].values

print(f"Total cells after filtering: {adata_sub.n_obs}")

# ==========================================
# 4.) COUNT REMAINING ACTIVE GENES/TYPES
# ==========================================
# How many cell types are TRULY present in this subset?
real_cell_types = adata_sub.obs['cell_type'].dropna().unique()
n_types = len(real_cell_types)

# How many genes are actually expressed? (Total sum in this subset > 0)
total_expression_per_gene = np.array(adata_sub.X.sum(axis=0)).flatten()
alive_genes_mask = total_expression_per_gene > 0

# Permanently discard non-expressed (dead) genes
adata_sub = adata_sub[:, alive_genes_mask].copy()
n_alive_genes = adata_sub.n_vars

print(f"\n--- FILTER RESULTS ---")
print(f"Real cell types in (21d, 3m, 22m): {n_types}")
print(f"Originally mapped U genes: {len(pos_genes)}")
print(f"TRULY ALIVE U genes in these stages: {n_alive_genes}")
print(f"Dead genes (turned off after embryonic stages): {len(pos_genes) - n_alive_genes}")
print(f"-----------------------------\n")

# ==========================================
# 5.) 3D EXPRESSION CALCULATION
# ==========================================
n_stages = len(target_stages)

matrix_3d = np.zeros((n_alive_genes, n_stages, n_types))
stage_type_presence = np.zeros((n_stages, n_types))

print("Calculating 3D matrix with purified data...")
for e_idx, stage in enumerate(target_stages):
    for t_idx, cell_type in enumerate(real_cell_types):
        
        mask = (adata_sub.obs['stage'] == stage) & (adata_sub.obs['cell_type'] == cell_type)
        
        if mask.sum() > 0:
            stage_type_presence[e_idx, t_idx] = 1
            sub_X = adata_sub.X[mask]
            
            if hasattr(sub_X, "sum"):
                sum_per_gene = np.array(sub_X.sum(axis=0)).flatten()
            else:
                sum_per_gene = np.sum(sub_X, axis=0)
                
            matrix_3d[:, e_idx, t_idx] = (sum_per_gene > 0).astype(int)

# ==========================================
# 6.) FAIR MATH AND FRACTIONS
# ==========================================
stage_count = matrix_3d.sum(axis=1) 
valid_stages_per_type = stage_type_presence.sum(axis=0)
valid_stages_per_type = np.where(valid_stages_per_type == 0, 1, valid_stages_per_type)

real_stage_fraction = stage_count / valid_stages_per_type
mean_exp_U = real_stage_fraction.mean(axis=1)

# ==========================================
# 7.) RESULTS AND HISTOGRAM
# ==========================================
threshold_10_percent = np.percentile(mean_exp_U, 10)
print(f"\nThe 10% threshold is: {threshold_10_percent:.4f}")

ind_low_U_genes = np.where(mean_exp_U <= threshold_10_percent)[0]
print(f"Highly specific genes (<= 10% threshold): {len(ind_low_U_genes)}")

plt.figure(figsize=(3, 2), dpi=600)
plt.hist(mean_exp_U, bins=50, log=True, color='violet')
plt.xlabel('Fraction of cell types\nexpressing the gene', fontweight='bold')
plt.ylabel('# Ut genes', fontsize=12, fontweight='bold')
plt.axvline(x=threshold_10_percent, color='indigo', linestyle='--', linewidth=0.8, label='10% threshold')
plt.legend()
plt.savefig(path_save_data + 'frac_Ut_genes_expressed_in_cell_types_adult.png', dpi=600, bbox_inches='tight')
plt.show()



#DEVELOPMENT

path_save_data='PATH_WHERE_YOU_KEEP_YOUR_DATASETS'
path_bulk='PATH_WHERE_YOU_KEEP_YOUR_DATASETS'

timepoints = ['10hpf', '12hpf', '14hpf', '16hpf', '19hpf', '24hpf', '2dpf', '3dpf', '5dpf', '10dpf']


#1.) Read the data
f=open(path_save_data+'genes_pseudo_more_cell_types.txt', 'r')
txt = f.read()
genes = txt.split('\n')
del txt, f
genes=np.delete(genes, len(genes)-1)
genes=np.array(genes)

f=open(path_save_data+'cell_type_unique.txt', 'r')
txt = f.read()
cell_type = txt.split('\n')
del txt, f
cell_type=np.delete(cell_type, len(cell_type)-1)
cell_type=np.array(cell_type)

f=open(path_save_data+'time_unique.txt', 'r')
txt = f.read()
time = txt.split('\n')
del txt, f
time=np.delete(time, len(time)-1)
time=np.array(time)

pseudo_matrix= np.genfromtxt(path_save_data+'pseudo_bulk_matrix.txt', dtype=None, encoding=None)
m_types= np.genfromtxt(path_save_data+'m_types.txt', dtype=None, encoding=None)


f=open(path_bulk+'hS_id_bulk.txt', 'r')
txt = f.read()
hS_id_bulk = txt.split('\n')
del txt, f
hS_id_bulk=np.delete(hS_id_bulk, len(hS_id_bulk)-1)
hS_id_bulk=np.array(hS_id_bulk)

f=open(path_bulk+'S_id_bulk.txt', 'r')
txt = f.read()
S_id_bulk = txt.split('\n')
del txt, f
S_id_bulk=np.delete(S_id_bulk, len(S_id_bulk)-1)
S_id_bulk=np.array(S_id_bulk)

f=open(path_bulk+'U_id_bulk.txt', 'r')
txt = f.read()
U_id_bulk = txt.split('\n')
del txt, f
U_id_bulk=np.delete(U_id_bulk, len(U_id_bulk)-1)
U_id_bulk=np.array(U_id_bulk)



#2. Fraction of cell types expressing the U genes in developement 
#We want to see if U genes in development are expressed in a high number of cell types
genes_U_single_cell, ind_genes_matrix_U_sc, j=np.intersect1d(genes, U_id_bulk, return_indices=True)
new_U_sc_matrix=pseudo_matrix[ind_genes_matrix_U_sc, :]

print(len(genes_U_single_cell))
print(len(genes_U_single_cell)/len(U_id_bulk))

n_genes = len(genes_U_single_cell) 
n_tiempos = len(time)
n_tipos_celulares = len(cell_type)

data_reshaped = new_U_sc_matrix.reshape((n_genes, n_tiempos, n_tipos_celulares))
print(n_tipos_celulares)

conteo_expresion = (data_reshaped > 0).sum(axis=1)

se_expresa_en_tipo = (conteo_expresion > 0).astype(int)

tiempos_validos_por_tipo = np.zeros(len(cell_type))
for i in range(len(m_types[:, 0])):
    tiempos_validos_por_tipo[i] = len(np.where(m_types[i, :]>0)[0])

total_tipos_validos = (tiempos_validos_por_tipo > 0).sum()

mean_exp_U = se_expresa_en_tipo.sum(axis=1) / total_tipos_validos


uno_por_ciento=np.percentile(mean_exp_U, 1)
print(uno_por_ciento)
ind_genes_U_bajo=np.where(mean_exp_U<uno_por_ciento)[0]
print(len(ind_genes_U_bajo))

plt.figure(figsize=(3, 2), dpi=600)
plt.hist(mean_exp_U, bins=100, log=True, color='darkorange')
plt.xlabel('Fraction of cell types\nexpressing the gene', fontweight='bold')
plt.ylabel('# U genes', fontsize=12, fontweight='bold')
plt.axvline(x=uno_por_ciento, color='sienna', linestyle='--', linewidth=0.8, label='1% threshold')
plt.legend()
plt.savefig(path_save_data + 'frac_Ut_genes_expressed_in_cell_types_development.png', dpi=600, bbox_inches='tight')
plt.show()





