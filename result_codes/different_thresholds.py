#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 16 14:42:00 2026

@author: alicia
"""

import scanpy as sc
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
import os
import warnings; warnings.simplefilter('ignore')
import anndata
from scipy.stats import pearsonr
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from sklearn.metrics.pairwise import pairwise_distances
from itertools import product
from scipy import stats
from scipy.stats import fisher_exact, false_discovery_control
import obonet
import networkx as nx

# ===================================================================================
# 1.0) PATHS AND DATA LOADING
# ===================================================================================

path_save_data='PATH_WHERE_YOU_KEEP_YOUR_DATASETS'
path='PATH_WHERE_YOU_KEEP_YOUR_DATASETS'
path_ontology='PATH_WHERE_YOU_KEEP_YOUR_DATASETS'   
path_go='PATH_WHERE_YOU_KEEP_YOUR_DATASETS'

file_path = path_save_data+"elife-30860-supp3-v1.tsv"

df = pd.read_csv(file_path, sep='\t')

genes=df['e85.Ensembl.Gene.ID'].to_numpy(dtype=str)
genes_name=df['Gene.name'].to_numpy(dtype=str)
gene_type=df['Gene.type'].to_numpy(dtype=str)
gene_type_unique, n_times_type=np.unique(gene_type, return_counts=True)

np.savetxt(path_save_data+'genes_bulk.txt', genes, fmt='%s', delimiter=',')
np.savetxt(path_save_data+'genes_bulk_name.txt', genes_name, fmt='%s', delimiter=',')

columns=df.columns.tolist()
select=columns[8:98]

time_unique=np.array(['0hpf', '0.75hpf', '2.25hpf', '3hpf', '4.3hpf', '5.25hpf', '6hpf', '8hpf', '10.3hpf', '16hpf', '19hpf', '24hpf', '30hpf', '36hpf', '2dpf', '3dpf', '4dpf', '5dpf'])
np.savetxt(path_save_data+'time_bulk.txt', time_unique, fmt='%s', delimiter=',')

embryo=np.array(['1', '2', '3', '4', '5'])
time_embryo = np.array([f"{time}-{embryo}" for time, embryo in product(time_unique, embryo)])

bulk_matrix = df[select].to_numpy()
total_tpm_per_stage1=np.sum(bulk_matrix, axis=0)

# ===================================================================================
# 2.0) CALCULATE MEAN BULK MATRIX AND TPM
# ===================================================================================

# We create a mean_bulk_matrix averaging the embryos in the same times 
mean_bulk_matrix=np.zeros((len(genes), len(time_unique)))
for k in range(len(genes)):
    count=0
    for i in range(len(time_unique)):
        sum_embryo=0
        for j in range(len(embryo)):
            sum_embryo=sum_embryo+bulk_matrix[k][int(count)]
            count=count+1
        mean_bulk_matrix[k][i]=sum_embryo/len(embryo)

# TPM mean expression
mean_bulk_matrix_tpm=np.zeros((len(genes), len(time_unique)))
for i in range(len(time_unique)):
    for j in range(len(genes)):
        mean_bulk_matrix_tpm[j][i]=mean_bulk_matrix[j][i]*1000000/(np.sum(mean_bulk_matrix[:, i]))

mean_bulk_matrix=mean_bulk_matrix_tpm
del mean_bulk_matrix_tpm

# Threshold: 1 TPM
mean_bulk_matrix_clean = np.where(mean_bulk_matrix < 1, 0, mean_bulk_matrix)
np.savetxt(path_save_data+'mean_bulk_matrix.txt', mean_bulk_matrix, fmt='%s', delimiter=',')
np.savetxt(path_save_data+'mean_bulk_matrix_clean.txt', mean_bulk_matrix_clean, fmt='%s', delimiter=',')

# We plot the distribution of stages in which genes are active
act_stage_per_gene=np.zeros(len(genes))
for i in range(len(genes)):
    non_null=len(np.where(mean_bulk_matrix_clean[i, :]>0)[0])
    act_stage_per_gene[i]=non_null

# Original Classification
null_genes_id=[]
hS_id=[]
U_id=[]
S_id=[]
for i in range(len(genes)):
    if act_stage_per_gene[i]==0:
        null_genes_id.append(genes[i])
    if (act_stage_per_gene[i]>0) & (act_stage_per_gene[i]<=8):
        hS_id.append(genes[i])
    if (act_stage_per_gene[i]>8) & (act_stage_per_gene[i]<=17):
        S_id.append(genes[i])
    if act_stage_per_gene[i]==18:
        U_id.append(genes[i])

np.savetxt(path_save_data+'hS_id_bulk.txt', hS_id, fmt='%s', delimiter=',')
np.savetxt(path_save_data+'S_id_bulk.txt', S_id, fmt='%s', delimiter=',')
np.savetxt(path_save_data+'U_id_bulk.txt', U_id, fmt='%s', delimiter=',')
np.savetxt(path_save_data+'null_genes.txt', null_genes_id, fmt='%s', delimiter=',')



# ===================================================================================
# 3.0) CONTINUOUS METRICS ANALYSIS
# ===================================================================================
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

print("\n=== Computing Continuous Metrics (Reviewer 1) ===")

# 3.1. METRIC: Fraction of stages detected (0 to 1)
fraction_stages = act_stage_per_gene / len(time_unique)

# 3.2. METRIC: Temporal Specificity Index (Tau / τ)

# Transformación logarítmica previa recomendada para Tau
log_matrix = np.log2(mean_bulk_matrix + 1)

x_max = np.max(log_matrix, axis=1)
with np.errstate(divide='ignore', invalid='ignore'):
    x_norm = log_matrix / x_max[:, None]
    tau = np.sum(1 - x_norm, axis=1) / (log_matrix.shape[1] - 1)

# 3.3. METRIC: Stability (CV of log(TPM+1))
log_tpm = np.log1p(mean_bulk_matrix)
std_log = np.std(log_tpm, axis=1)
mean_log = np.mean(log_tpm, axis=1)
with np.errstate(divide='ignore', invalid='ignore'):
    cv_log = std_log / mean_log
cv_log[np.isnan(cv_log)] = np.nan # Genes with 0 expression will be NaN

# Calculate percentiles BEFORE plotting
x1 = np.percentile(fraction_stages, 33.33) # 33rd percentile (Tertile 1)
x2 = np.percentile(fraction_stages, 66.66) # 66th percentile (Tertile 2)
print(f"33rd Percentile: {x1:.2f} (approx {x1*18:.1f} stages)")
print(f"66th Percentile: {x2:.2f} (approx {x2*18:.1f} stages)")
plt.figure(figsize=(4.5, 3), dpi=600)
# Plo the histogram of the continuous variable
plt.hist(fraction_stages, bins=18, color='lightgrey', alpha=0.8)
# Add vertical lines for the percentiles
plt.axvline(x=x1, color='tomato', linestyle='--', linewidth=2, label=f'33rd Percentile')
plt.axvline(x=x2, color='royalblue', linestyle='--', linewidth=2, label=f'66th Percentile')
plt.xlabel('Fraction of expressed stages', fontsize=14, fontweight='bold')
plt.ylabel('# genes', fontsize=14, fontweight='bold')
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(fontsize=10, loc='upper left')
plt.tight_layout()
plt.savefig(path_save_data + 'empirical_cutoff_justification_histogram.png', dpi=600, bbox_inches='tight')
plt.show()


pearsonr(fraction_stages, tau)
pearsonr(fraction_stages, cv_log)

# Create a DataFrame for easy Seaborn plotting
# Map genes to their original classes
gene_class_array = np.full(len(genes), 'Filtered', dtype=object)
gene_class_array[np.isin(genes, U_id)] = 'U'
gene_class_array[np.isin(genes, S_id)] = 'S'
gene_class_array[np.isin(genes, hS_id)] = 'hS'

df_continuous = pd.DataFrame({
    'Gene': genes,
    'Class': gene_class_array,
    'Fraction_Stages': fraction_stages,
    'Tau': tau,
    'CV_logTPM': cv_log
})


df_plot = df_continuous[df_continuous['Class'].isin(['U', 'S', 'hS'])]
df_plot['Class'] = pd.Categorical(df_plot['Class'], categories=['U', 'S', 'hS'], ordered=True)

#figure
fig, axes = plt.subplots(1, 2, figsize=(10, 4.5), dpi=600)
classes = ['U', 'S', 'hS']
colors = ['royalblue', 'mediumturquoise', 'tomato']
data_tau = [df_plot[df_plot['Class'] == c]['Tau'].dropna().values for c in classes]
data_cv = [df_plot[df_plot['Class'] == c]['CV_logTPM'].dropna().values for c in classes]

def style_violin(parts, colors):
    for pc, color in zip(parts['bodies'], colors):
        pc.set_facecolor(color)
        pc.set_edgecolor(None)
        pc.set_alpha(1.0) 
        
    for partname in ('cbars', 'cmins', 'cmaxes', 'cquantiles'):
        if partname in parts:
            vp = parts[partname]
            vp.set_edgecolor('black')
            vp.set_linewidth(0.9)
            
parts_tau = axes[0].violinplot(
    data_tau, 
    positions=[1, 2, 3], 
    showmeans=False, 
    showmedians=False,
    showextrema=True, 
    quantiles=[[0.25, 0.5, 0.75]] * len(classes)
)
style_violin(parts_tau, colors)
axes[0].set_ylabel(r'Temporal specificity ($\tau$)', fontsize=20, fontweight='bold')
axes[0].set_xticks([1, 2, 3])
axes[0].set_xticklabels(classes)

parts_cv = axes[1].violinplot(
    data_cv, 
    positions=[1, 2, 3], 
    showmeans=False, 
    showmedians=False,
    showextrema=True, 
    quantiles=[[0.25, 0.5, 0.75]] * len(classes)
)
style_violin(parts_cv, colors)
axes[1].set_ylabel('Temporal variability (CV)', fontsize=20, fontweight='bold')
axes[1].set_xticks([1, 2, 3])
axes[1].set_xticklabels(classes)

for ax in axes:
    ax.tick_params(axis='x', labelsize=24)
    for label in ax.get_xticklabels():
        label.set_fontweight('bold')
    
    ax.tick_params(axis='y', labelsize=15)
plt.tight_layout()
plt.savefig(path_save_data + 'continuous_metrics_vs_discrete_bins.png', dpi=600, bbox_inches='tight')
plt.show()


# ===================================================================================
# 4) SENSITIVITY ANALYSIS: BOUNDARIES 
# ===================================================================================
parameters_to_test = [
    (0.5, 8, 18),  # 1. Lower TPM threshold
    (1.0, 8, 18),  # 2. ORIGINAL (Reference)
    (2.0, 8, 18),  # 3. Higher TPM threshold
]

labels = ["TPM > 0.5", "Original\nTPM > 1", "TPM > 2.0"]

stability_U = []
stability_S = []
stability_hS = []

for tpm_thresh, hs_max, u_min in parameters_to_test:
    temp_matrix_clean = np.where(mean_bulk_matrix < tpm_thresh, 0, mean_bulk_matrix)
    temp_act_stage_per_gene = np.sum(temp_matrix_clean > 0, axis=1)
    
    temp_U_id = genes[temp_act_stage_per_gene >= u_min]
    temp_hS_id = genes[(temp_act_stage_per_gene > 0) & (temp_act_stage_per_gene <= hs_max)]
    temp_S_id = genes[(temp_act_stage_per_gene > hs_max) & (temp_act_stage_per_gene < u_min)]
    
    stayed_U = len(np.intersect1d(U_id, temp_U_id)) / len(U_id) * 100 if len(U_id) > 0 else 0
    stayed_S = len(np.intersect1d(S_id, temp_S_id)) / len(S_id) * 100 if len(S_id) > 0 else 0
    stayed_hS = len(np.intersect1d(hS_id, temp_hS_id)) / len(hS_id) * 100 if len(hS_id) > 0 else 0
    
    stability_U.append(stayed_U)
    stability_S.append(stayed_S)
    stability_hS.append(stayed_hS)

x = np.arange(len(labels))
width = 0.25  

fig, ax = plt.subplots(figsize=(7.5, 4.5), dpi=600)
rects1 = ax.bar(x - width, stability_U, width, label='U', color='royalblue')
rects2 = ax.bar(x, stability_S, width, label='S', color='mediumturquoise')
rects3 = ax.bar(x + width, stability_hS, width, label='hS', color='tomato')
ax.set_ylabel('% genes retained\nfrom original threshold', fontsize=20, fontweight='bold')
ax.set_xticks(x)
ax.set_xticklabels(labels, rotation=0, fontsize=18, fontweight='bold')
ax.set_ylim(0, 110) 
plt.yticks(fontsize=14)
ax.legend(fontsize=14, loc='upper left', bbox_to_anchor=(1.05, 1))
ax.axhline(90, color='gray', linestyle='--', linewidth=1)
ax.axhline(70, color='lightgray', linestyle='--', linewidth=1)
plt.tight_layout()
plt.savefig(path_save_data + 'membership_stability_plot.png', dpi=600, bbox_inches='tight')
plt.show()

# ===================================================================================
# 5) LOAD GENE AGE AND ONTOLOGY DATA
# ===================================================================================
df_gene_origin = pd.read_csv(path + "gene_origin.csv")
genes_df = np.array(list(df_gene_origin['ensembl_gene_id']))
age_df = np.array(list(df_gene_origin['gene_age']))

# GO functions
def flatten_and_unique(nested_list):
    result = []
    for sublist in nested_list:
        unique_items = set()
        for item in sublist:
            unique_items.update(item)  
        result.append(list(unique_items)) 
    return result

def create_association_matrix_gene_go_term(gene_subset, GO_specific_terms, pathlist_specific_terms, GO_specific_terms_descrip):
    common_genes = np.intersect1d(gene_subset, genes_id_type)
    GO_specific_terms = np.array(GO_specific_terms)
    matrix = np.zeros((len(common_genes), len(GO_specific_terms)))
    for i in range(len(common_genes)):
        ind_gene = np.where(genes_id_type == common_genes[i])[0]
        for j in range(len(go_type[int(ind_gene)])):
            go_ind = np.where(GO_specific_terms == go_type[int(ind_gene)][j])[0]
            if len(go_ind) > 0:
                ind_matrix = np.isin(GO_specific_terms, pathlist_specific_terms[int(go_ind)])
                matrix[i, ind_matrix] = 1
    return matrix

def enrichement_go(big_matrix_all_genes, submatrix_gene_go, go_list, go_list_term, label):
    go_n_times_all_genes = np.sum(big_matrix_all_genes, axis=0)
    odd_ratio_enrich = np.zeros(len(go_list))
    p_value_enrich = np.zeros(len(go_list))
    n_genes_subset = len(submatrix_gene_go[:, 0])
    
    go_enrich_fisher_genes_subset = []
    go_term_enrich = []
    p_value_enriched_go_term = []
    n_genes = []
    n_genes_subset_associated_go = []
    subset_analyzed = []
    
    for fen in range(len(go_list)):
        go_n_times_subset = np.sum(submatrix_gene_go[:, fen])
        tabla = [[go_n_times_subset, n_genes_subset - go_n_times_subset],
                 [go_n_times_all_genes[fen], len(big_matrix_all_genes[:, 0]) - go_n_times_all_genes[fen]]]
        odd_ratio_enrich[fen], p_value_enrich[fen] = fisher_exact(tabla, alternative="greater") 
        if p_value_enrich[fen] < 0.001:
            go_enrich_fisher_genes_subset.append(go_list[fen])
            go_term_enrich.append(go_list_term[fen])
            p_value_enriched_go_term.append(p_value_enrich[fen])
            n_genes.append(go_n_times_all_genes[fen])
            n_genes_subset_associated_go.append(go_n_times_subset)
            subset_analyzed.append(label)
            
    return np.array(subset_analyzed), np.array(go_enrich_fisher_genes_subset), np.array(go_term_enrich), np.array(p_value_enriched_go_term), np.array(n_genes), np.array(n_genes_subset_associated_go)

def convert_dictionary_to_array(dictionary):
    return np.array(list(dictionary.items()))

# Load GO Network
url = path_go + 'go.obo'
graph = obonet.read_obo(url)
id_to_name = {id_: data.get('name') for id_, data in graph.nodes(data=True)}
name_to_id = {data['name']: id_ for id_, data in graph.nodes(data=True) if 'name' in data}
gene_description = {data['def']: id_ for id_, data in graph.nodes(data=True) if 'def' in data}

GO_term = convert_dictionary_to_array(id_to_name)
go_description = convert_dictionary_to_array(gene_description)

# Generate biological process paths
pathlist_bio_process = []
GO_bio_process = []
GO_bio_process_descrip = []
for i in range(len(GO_term)):
    start = id_to_name[GO_term[i][0]]
    if 'biological_process' in name_to_id and start in name_to_id:
        try:
            paths = nx.all_simple_paths(graph, source=name_to_id[start], target=name_to_id['biological_process'])
            innerlist = list(paths)
            if len(innerlist) > 0:
                pathlist_bio_process.append(innerlist)
                GO_bio_process.append(GO_term[i][0])
                GO_bio_process_descrip.append(GO_term[i][1])
        except nx.NetworkXNoPath:
            pass
        except nx.NodeNotFound:
            pass


pathlist_cell_comp=[]
GO_cell_comp=[]
GO_cell_comp_descrip=[]
for i in range(len(GO_term)):
    start=id_to_name[GO_term[i][0]]
    
    paths = nx.all_simple_paths(
        graph,
        source=name_to_id[start],
        target=name_to_id['cellular_component']
    )
    innerlist = []
    for path in paths:
        innerlist.append(path)
        
    if len(innerlist)>0:
        pathlist_cell_comp.append(innerlist)
        GO_cell_comp.append(GO_term[i][0])
        GO_cell_comp_descrip.append(GO_term[i][1])
    
pathlist_molecular=[]
GO_molecular=[]
GO_molecular_descrip=[]
for i in range(len(GO_term)):
    start=id_to_name[GO_term[i][0]]
    
    paths = nx.all_simple_paths(
        graph,
        source=name_to_id[start],
        target=name_to_id['molecular_function']
    )
    innerlist = []
    for path in paths:
        innerlist.append(path)
    
    if len(innerlist)>0:
        pathlist_molecular.append(innerlist)
        GO_molecular.append(GO_term[i][0])
        GO_molecular_descrip.append(GO_term[i][1])


# Load ZFIN annotations
file_path = path_go + "zfin.gaf"
df_zfin = pd.read_csv(file_path, sep="\t", comment="!", header=None, on_bad_lines="skip", engine="python")
gene_id = np.array(df_zfin.iloc[:, 1], dtype=str)
GOterm_association = np.array(df_zfin.iloc[:, 4], dtype=str)
genes_id_type = np.unique(gene_id)

go_type = []
for i in range(len(genes_id_type)):
    ind_genes = np.where(gene_id == genes_id_type[i])[0]
    inner_list = [GOterm_association[idx] for idx in ind_genes]
    go_type.append(inner_list)

# Match ZFIN to ENSEMBL
f = open(path_ontology + 'zfin_to_ENS.txt', 'r')
zfin_to_ENS = f.read().split('\n')[:-1]
f.close()

zfin_all = []
ENS_all = []
for elemento in zfin_to_ENS:
    partes = elemento.split("\t")
    for p in partes:
        if p.startswith("ZDB"): zfin_all.append(p)
        if p.startswith("ENS"): ENS_all.append(p)

ENS_all = np.array(ENS_all)
zfin_all = np.array(zfin_all)

common_genes_zfin, ind_zfin_go, ind_ens = np.intersect1d(genes_id_type, zfin_all, return_indices=True)
genes_ens = ENS_all[ind_ens]
genes_unique = np.unique(genes_ens)

go_type_list_all_genes = [go_type[idx] for idx in ind_zfin_go]
go_type_list = []
for i in range(len(genes_unique)):
    ind_gene = np.where(genes_ens == genes_unique[i])[0]
    if len(ind_gene) > 1:
        go_type_list.append(np.unique(np.concatenate((go_type_list_all_genes[ind_gene[0]], go_type_list_all_genes[ind_gene[1]]))))
    else:
        go_type_list.append(go_type_list_all_genes[ind_gene[0]])

go_type = go_type_list
genes_id_type = genes_unique

pathlist_bio_process_unique = flatten_and_unique(pathlist_bio_process)
pathlist_molecular_unique=flatten_and_unique(pathlist_molecular)
pathlist_cell_comp_unique=flatten_and_unique(pathlist_cell_comp)

# ===================================================================================
# 6) FAST ROBUSTNESS CHECK: AGE VS MAX STAGE (Sliding Window) ACROSS THRESHOLDS
# ===================================================================================
print("\n=== Running Fast Sliding Window Age Analysis Across Thresholds ===")

tpm_thresholds = [0.5, 1.0, 2.0]
sw = 400 # Sliding window size
fig, axes = plt.subplots(1, 3, figsize=(16, 4), dpi=600, sharey=True)

common, ind_matrix_all, ind_age_all = np.intersect1d(genes, genes_df, return_indices=True)

for idx, tpm_thresh in enumerate(tpm_thresholds):
    ax = axes[idx]
    
    # 1. Dynamically re-classify genes
    temp_matrix_clean = np.where(mean_bulk_matrix < tpm_thresh, 0, mean_bulk_matrix)
    temp_act_stage_per_gene = np.sum(temp_matrix_clean > 0, axis=1)
    
    temp_U_id = genes[temp_act_stage_per_gene == 18]
    temp_S_id = genes[(temp_act_stage_per_gene > 8) & (temp_act_stage_per_gene <= 17)]
    temp_hS_id = genes[(temp_act_stage_per_gene > 0) & (temp_act_stage_per_gene <= 8)]
    
    # Helper function to process each class quickly
    def process_class_sw(class_ids):
        # Filter to only genes in this class
        is_in_class = np.isin(genes[ind_matrix_all], class_ids)
        ind_matrix_class = ind_matrix_all[is_in_class]
        
        # Get ages and clean '>4290'
        raw_ages = age_df[ind_age_all[is_in_class]]
        ages_class = np.array(['4290' if x == '>4290' else x for x in raw_ages], dtype=int)
        
        # Get max expression stage
        bulk_class = mean_bulk_matrix[ind_matrix_class, :]
        stage_max_class = np.argmax(bulk_class, axis=1)
        
        # Sort by stage of max expression
        sort_idx = np.argsort(stage_max_class)
        stage_max_sorted = stage_max_class[sort_idx]
        ages_sorted = ages_class[sort_idx]
        
        # Apply rolling window
        sw_stage = pd.Series(stage_max_sorted).rolling(window=sw, center=False).mean()
        sw_age = pd.Series(ages_sorted).rolling(window=sw, center=False).mean()
        
        return sw_stage, sw_age

    # Process U, S, and hS
    sw_stage_U, sw_age_U = process_class_sw(temp_U_id)
    sw_stage_S, sw_age_S = process_class_sw(temp_S_id)
    sw_stage_hS, sw_age_hS = process_class_sw(temp_hS_id)

    # Plot on the corresponding subplot
    ax.scatter(sw_stage_U, sw_age_U, s=1, color='royalblue', label='U')
    ax.scatter(sw_stage_S, sw_age_S, s=1, color='mediumturquoise', label='S')
    ax.scatter(sw_stage_hS, sw_age_hS, s=1, color='tomato', label='hS')
    
    ax.set_title(f'TPM > {tpm_thresh}', fontsize=16, fontweight='bold')
    ax.set_xlabel('Stage max. expression', fontsize=14, fontweight='bold')
    ax.set_xticks(np.arange(len(time_unique)))
    ax.set_xticklabels(time_unique, fontsize=10, rotation=90)
    
    if idx == 0:
        ax.set_ylabel('Age (myr)', fontsize=14, fontweight='bold')
    
    # Only add legend to the last plot to save space
    if idx == 2:
        ax.legend(markerscale=5, loc='upper right', fontsize=12)

plt.tight_layout()
plt.savefig(path_save_data + f'age_vs_max_stage_sw_{sw}_ALL_THRESHOLDS.png', dpi=600, bbox_inches='tight')
plt.show()
print("Saved sliding window plot across all thresholds.")
# ===================================================================================

# ===================================================================================
# 7) GROUPED BOXPLOT WITH CUSTOM COLORS: GENE AGE ACROSS THRESHOLDS 
# ===================================================================================
import matplotlib.patches as mpatches

print("\n=== Generating Custom Colored Grouped Boxplot ===")

tpm_thresholds = [0.5, 1.0, 2.0]

# Dictionary to store the arrays before plotting
age_dict = {'U': {}, 'S': {}, 'hS': {}}

# Pre-calculate intersection (if not done already)
common, ind_matrix_all, ind_age_all = np.intersect1d(genes, genes_df, return_indices=True)

for tpm_thresh in tpm_thresholds:
    # Dynamically re-classify genes
    temp_matrix_clean = np.where(mean_bulk_matrix < tpm_thresh, 0, mean_bulk_matrix)
    temp_act_stage_per_gene = np.sum(temp_matrix_clean > 0, axis=1)
    
    temp_U_id = genes[temp_act_stage_per_gene == 18]
    temp_S_id = genes[(temp_act_stage_per_gene > 8) & (temp_act_stage_per_gene <= 17)]
    temp_hS_id = genes[(temp_act_stage_per_gene > 0) & (temp_act_stage_per_gene <= 8)]
    
    # Helper function to extract ages
    def get_ages_for_class(class_ids):
        is_in_class = np.isin(genes[ind_matrix_all], class_ids)
        raw_ages = age_df[ind_age_all[is_in_class]]
        return np.array(['4290' if x == '>4290' else x for x in raw_ages], dtype=int)

    # Store ages for each threshold
    age_dict['U'][tpm_thresh] = get_ages_for_class(temp_U_id)
    age_dict['S'][tpm_thresh] = get_ages_for_class(temp_S_id)
    age_dict['hS'][tpm_thresh] = get_ages_for_class(temp_hS_id)

# Organize data into a flat list for Matplotlib (Order: U(0.5,1,2), S(0.5,1,2), hS(0.5,1,2))
data = [
    age_dict['U'][0.5], age_dict['U'][1.0], age_dict['U'][2.0],
    age_dict['S'][0.5], age_dict['S'][1.0], age_dict['S'][2.0],
    age_dict['hS'][0.5], age_dict['hS'][1.0], age_dict['hS'][2.0]
]

# Define physical positions on the X-axis to create gaps between the U, S, and hS groups
positions = [1, 1.8, 2.6,   4.6, 5.4, 6.2,   8.2, 9.0, 9.8]

fig, ax = plt.subplots(figsize=(6, 4), dpi=600)

# Create the boxplot
box = ax.boxplot(data, positions=positions, widths=0.6, patch_artist=True,
                 flierprops=dict(marker='o', color='gray', markersize=1.5, alpha=0.5),
                 boxprops=dict(edgecolor="black", linewidth=1.2))

# ---------------------------------------------------------
# CUSTOM COLOR PALETTE (Dark -> Medium -> Light per class)
# ---------------------------------------------------------
custom_colors = [
    '#002b5e', '#0050b3', '#4da6ff',  # U  : Dark Navy -> Royal Blue -> Light Blue
    '#006d77', '#2a9d8f', '#88d4cb',  # S  : Dark Teal -> Turquoise  -> Light Cyan
    '#9b2226', '#e63946', '#ff8fa3'   # hS : Dark Red  -> Bright Red -> Light Red
]

# Apply the colors to each box
for patch, color in zip(box["boxes"], custom_colors):
    patch.set_facecolor(color)

# Style medians, whiskers, and caps
plt.setp(box["medians"], color="white", linewidth=1.5) # White medians contrast nicely with dark boxes
plt.setp(box["whiskers"], color="black", linestyle="--", linewidth=1)
plt.setp(box["caps"], color="black", linewidth=1)

# Formatting the plot
plt.ylabel("Age (myr)", fontsize=18, fontweight='bold')
# plt.title("Evolutionary Age Stability Across TPM Thresholds", fontsize=20, fontweight='bold')

# Set X-ticks directly in the center of each group
ax.set_xticks([1.8, 5.4, 9.0])
ax.set_xticklabels(['U', 'S', 'hS'], fontsize=18, fontweight='bold')
plt.yticks(fontsize=14)
plt.ylim(0, 2800) 

# Create a custom legend representing the intensities
leg1 = mpatches.Patch(color='#666666', label='TPM > 0.5')
leg2 = mpatches.Patch(color='#999999', label='TPM > 1.0')
leg3 = mpatches.Patch(color='#cccccc', label='TPM > 2.0')
ax.legend(handles=[leg1, leg2, leg3], 
          fontsize=12, loc='upper right')

plt.tight_layout()
plt.savefig(path_save_data + 'age_boxplot_grouped_custom_colors.png', dpi=600, bbox_inches='tight')
plt.show()
# ===================================================================================


# ===================================================================================
# 8) DOWNSTREAM ROBUSTNESS ANALYSIS (GO Enrichment across thresholds)
# ===================================================================================
# PRE-CALCULATE BACKGROUND MATRICES (Fix to speed up the loop)
print("Calculating background GO matrices for BP, MF, and CC (this may take a moment)...")

big_matrix_bio_process = create_association_matrix_gene_go_term(genes_id_type, GO_bio_process, pathlist_bio_process_unique, GO_bio_process_descrip)
big_matrix_molecular = create_association_matrix_gene_go_term(genes_id_type, GO_molecular, pathlist_molecular_unique, GO_molecular_descrip)
big_matrix_cell_comp = create_association_matrix_gene_go_term(genes_id_type, GO_cell_comp, pathlist_cell_comp_unique, GO_cell_comp_descrip)

tpm_thresholds = [0.5, 1.0, 2.0]

# Dictionary to store the significant GO terms for comparison later across all 3 ontologies
significant_go_terms = {
    'bio_process': {'U': {}, 'S': {}, 'hS': {}},
    'molecular':   {'U': {}, 'S': {}, 'hS': {}},
    'cell_comp':   {'U': {}, 'S': {}, 'hS': {}}
}

for tpm_thresh in tpm_thresholds:
    print(f"\n--- Running downstream GO analysis for TPM > {tpm_thresh} ---")
    
    # 1. Dynamically re-classify genes
    temp_matrix_clean = np.where(mean_bulk_matrix < tpm_thresh, 0, mean_bulk_matrix)
    temp_act_stage_per_gene = np.sum(temp_matrix_clean > 0, axis=1)
    
    temp_U_id = genes[temp_act_stage_per_gene == 18]
    temp_S_id = genes[(temp_act_stage_per_gene > 8) & (temp_act_stage_per_gene <= 17)]
    temp_hS_id = genes[(temp_act_stage_per_gene > 0) & (temp_act_stage_per_gene <= 8)]
    
    # --- GO ENRICHMENT (BP, MF, CC) ---
    genes_analyze_list = [temp_U_id, temp_S_id, temp_hS_id]
    labels_list = ['U', 'S', 'hS']

    for i in range(3):
        genes_analyze = genes_analyze_list[i]
        current_class = labels_list[i]
        
        # Biological Process
        matrix_subset_bio = create_association_matrix_gene_go_term(genes_analyze, GO_bio_process, pathlist_bio_process_unique, GO_bio_process_descrip)
        subset_b, enrich_descrip_bio, enrich_term_b, p_val_b, n_g_b, n_g_sub_b = enrichement_go(
            big_matrix_bio_process, matrix_subset_bio, GO_bio_process_descrip, GO_bio_process, current_class)
        significant_go_terms['bio_process'][current_class][tpm_thresh] = set(enrich_descrip_bio)

        # Molecular Function
        matrix_subset_mol = create_association_matrix_gene_go_term(genes_analyze, GO_molecular, pathlist_molecular_unique, GO_molecular_descrip)
        subset_m, enrich_descrip_mol, enrich_term_m, p_val_m, n_g_m, n_g_sub_m = enrichement_go(
            big_matrix_molecular, matrix_subset_mol, GO_molecular_descrip, GO_molecular, current_class)
        significant_go_terms['molecular'][current_class][tpm_thresh] = set(enrich_descrip_mol)

        # Cellular Component
        matrix_subset_cell = create_association_matrix_gene_go_term(genes_analyze, GO_cell_comp, pathlist_cell_comp_unique, GO_cell_comp_descrip)
        subset_c, enrich_descrip_cell, enrich_term_c, p_val_c, n_g_c, n_g_sub_c = enrichement_go(
            big_matrix_cell_comp, matrix_subset_cell, GO_cell_comp_descrip, GO_cell_comp, current_class)
        significant_go_terms['cell_comp'][current_class][tpm_thresh] = set(enrich_descrip_cell)

# ===================================================================================
# 9) PRINT THE OVERLAP OF GO TERMS ACROSS THRESHOLDS (ALL ONTOLOGIES)
# ===================================================================================
print("\n" + "="*50)
print("             GO TERM STABILITY RESULTS")
print("="*50)

ontologies_to_print = [
    ('Biological Process', 'bio_process'),
    ('Molecular Function', 'molecular'),
    ('Cellular Component', 'cell_comp')
]

for ont_name, ont_key in ontologies_to_print:
    print(f"\n---> {ont_name.upper()} <---")
    
    for gene_class in ['U', 'S', 'hS']:
        terms_05 = significant_go_terms[ont_key][gene_class].get(0.5, set())
        terms_10 = significant_go_terms[ont_key][gene_class].get(1.0, set()) # Original
        terms_20 = significant_go_terms[ont_key][gene_class].get(2.0, set())
        
        if len(terms_10) > 0:
            overlap_05 = len(terms_10.intersection(terms_05)) / len(terms_10) * 100
            overlap_20 = len(terms_10.intersection(terms_20)) / len(terms_10) * 100
            
            print(f"  Class {gene_class} (Original significant terms: {len(terms_10)}):")
            print(f"    - {overlap_05:.1f}% remain significant at TPM > 0.5")
            print(f"    - {overlap_20:.1f}% remain significant at TPM > 2.0")
        else:
            print(f"  Class {gene_class}: No significant terms found at TPM > 1.0.")
            
            
# ===================================================================================
# 10. Terms lost in the S class
# ===================================================================================

# Extraemos los sets de la ontología Cellular Component para la clase S
terms_10_S_cc = significant_go_terms['cell_comp']['S'].get(1.0, set())
terms_20_S_cc = significant_go_terms['cell_comp']['S'].get(2.0, set())

# Buscamos los que estaban en 1.0 pero NO están en 2.0
dropped_terms = terms_10_S_cc - terms_20_S_cc

print("\n" + "="*60)
print(f" TÉRMINOS GO QUE SE PIERDEN EN LA CLASE 'S' A TPM > 2.0 (CC)")
print("="*60)
print(f"Total de términos perdidos: {len(dropped_terms)}\n")

for i, term in enumerate(dropped_terms, 1):
    print(f"{i}. {term}")

print("\n" + "="*60)
