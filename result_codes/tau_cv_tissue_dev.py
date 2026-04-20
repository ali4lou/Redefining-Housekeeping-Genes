
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import warnings; warnings.simplefilter('ignore')
from scipy.stats import pearsonr
from itertools import product
from scipy.stats import ks_2samp
import seaborn as sns
from scipy.stats import spearmanr
import obonet
import networkx as nx


path_figures='PATH_WHERE_YOU_KEEP_YOUR_DATASETS'
#========================================================================================
#TISSUE
#========================================================================================
path_data='PATH_WHERE_YOU_KEEP_YOUR_DATASETS'
path_save_data='PATH_WHERE_YOU_KEEP_YOUR_DATASETS'
path_bulk_initial='PATH_WHERE_YOU_KEEP_YOUR_DATASETS'

tissue=['brain', 'gill', 'heart', 'intestine', 'kidney', 'liver', 'muscle', 'spleen']
temp = ['10', '18', '28']
# temp=['28']
label= [t + tp for t in tissue for tp in temp]


data=pd.read_csv(path_data+'%s.expression.txt' %label[0], sep='\t')
genes=np.array(list(data['Gene']))
genes, ind_genes=np.unique(genes, return_index=True)
    
genes_tis=genes
del genes

#1.) We are going to analyse temperature 28ºC
tissue=['brain', 'gill', 'heart', 'intestine', 'kidney', 'liver', 'muscle', 'spleen']
temp=['28']
label= [t + tp for t in tissue for tp in temp]

matrix=np.zeros((len(genes_tis), len(tissue)))
for i in range(len(label)):
    data=pd.read_csv(path_data+'%s.expression.txt' %label[i], sep='\t')
    genes=np.array(list(data['Gene']))
    genes=genes[ind_genes]
    counts=np.array(list(data['RPKM']))
    counts_new=counts[ind_genes]
    matrix[:, i]=counts_new

np.sum(matrix[:, 2])

matrix_TPM=np.zeros((len(genes), len(tissue)))
sum_tissue=np.sum(matrix, axis=0)

#2.) Convert matrix from RPKM to TPM
for i in range(len(tissue)):
    for j in range(len(genes)):
        matrix_TPM[j][i]=matrix[j][i]*1000000/sum_tissue[i]

matrix=matrix_TPM
del matrix_TPM

sum_tissue=np.sum(matrix, axis=0)


#2.1.) Total expression level per tissue
n_genes_per_tis=np.zeros(len(tissue))
for i in range(len(tissue)):
    non_null=len(np.where(matrix[:, i]>1)[0])
    n_genes_per_tis[i]=non_null

plt.figure(figsize=(11, 6), dpi=600)
plt.bar(tissue, n_genes_per_tis, color='deepskyblue', edgecolor='black')
plt.xticks(ticks=np.arange(0, len(tissue)), labels=tissue, rotation=90, fontsize=25)
plt.yticks(fontsize=25)
plt.ylabel('# genes\n (TPM>1)', fontsize=30, fontweight='bold')
# plt.ylabel('# genes (TPM>1)', fontsize=30, fontweight='bold')
plt.savefig(path_save_data+'total_genes_per_tissue_tpm_1.png', dpi=600, bbox_inches='tight')
plt.show()



#3.) Active tissues
act_tis_per_gene=np.zeros(len(genes))
unique_tissue=[]
unique_genes=[]
for i in range(len(genes)):
    non_null=len(np.where(matrix[i, :]>1)[0])
    act_tis_per_gene[i]=non_null
    if non_null==1:
        ind_tis=np.where(matrix[i, :]>1)[0]
        unique_tissue.append(tissue[int(ind_tis)])
        unique_genes.append(genes[i])



n_stages, n_times=np.unique(act_tis_per_gene, return_counts=True)
n_stages=np.array(n_stages, dtype=int)
plt.figure(figsize=(11, 6), dpi=600)
plt.bar(n_stages, n_times, color='darkorange', edgecolor='black')
plt.xticks(ticks=np.arange(0, len(n_stages)), labels=n_stages, rotation=90, fontsize=25)
plt.yticks(fontsize=25)
plt.xlabel('# tissues expressed', fontsize=30, fontweight='bold')
plt.ylabel('# genes (TPM>1)', fontsize=30, fontweight='bold')
plt.savefig(path_save_data+'act_tissues_per_gene.png', dpi=600, bbox_inches='tight')
plt.show()


unique_tissues_times, n_times_unique_genes=np.unique(np.array(unique_tissue), return_counts=True)
plt.figure(figsize=(11, 6), dpi=600)
plt.bar(unique_tissues_times, n_times_unique_genes/n_genes_per_tis, color='violet', edgecolor='black')
plt.xticks(ticks=np.arange(0, len(unique_tissues_times)), labels=unique_tissues_times, rotation=90, fontsize=25)
plt.yticks(fontsize=25)
plt.ylabel('% hS genes\n(just in one tissue)', fontsize=30, fontweight='bold')
# plt.ylabel('# genes (TPM>1)', fontsize=30, fontweight='bold')
plt.savefig(path_save_data+'hS_genes_per_tissue.png', dpi=600, bbox_inches='tight')
plt.show()


#4.) We divided the genes in specific, high-specific and ubiqutouly expressed
null_genes=[]
hS_tis=[]
U_tis=[]
S_tis=[]
act_tis_S=[]
act_tis_hS=[]
act_tis_U=[]
for i in range(len(genes)):
    if act_tis_per_gene[i]==0:
        null_genes.append(genes[i])
    if (act_tis_per_gene[i]>0) & (act_tis_per_gene[i]<=4):
        hS_tis.append(genes[i])
        act_tis_hS.append(act_tis_per_gene[i])
    if (act_tis_per_gene[i]>4) & (act_tis_per_gene[i]<=7):
        S_tis.append(genes[i])
        act_tis_S.append(act_tis_per_gene[i])
    if act_tis_per_gene[i]==8:
        U_tis.append(genes[i])
        act_tis_U.append(act_tis_per_gene[i])


frac_genes=np.array([len(null_genes)/len(genes), len(U_tis)/len(genes), len(S_tis)/len(genes), len(hS_tis)/len(genes)])
groups=['NE', 'Ut', 'St', 'hSt']
color=['gray', 'royalblue', 'mediumturquoise', 'tomato']

plt.figure(figsize=(4, 3), dpi=600)
plt.bar(groups, frac_genes, color=color)
plt.xticks(ticks=np.arange(len(groups)), labels=groups, fontsize=20, fontweight='bold')
plt.yticks(fontsize=15)
plt.ylabel('Gene fraction', fontsize=20, fontweight='bold')
plt.savefig(path_save_data+'gene_groups.png', dpi=600, bbox_inches='tight')
plt.show()


frac_tissue=act_tis_per_gene/len(tissue)
plt.hist(frac_tissue)



# Calculate percentiles BEFORE plotting
x1 = np.percentile(frac_tissue, 33.33) # 33rd percentile (Tertile 1)
x2 = np.percentile(frac_tissue, 66.66) # 66th percentile (Tertile 2)
print(f"33rd Percentile: {x1:.2f} (approx {x1*8:.1f} stages)")
print(f"66th Percentile: {x2:.2f} (approx {x2*8:.1f} stages)")
plt.figure(figsize=(4.5, 3), dpi=600)
# Plo the histogram of the continuous variable
plt.hist(frac_tissue, bins=8, color='lightgrey', alpha=0.8)
# Add vertical lines for the percentiles
plt.axvline(x=x1, color='tomato', linestyle='--', linewidth=2, label=f'33rd Percentile')
plt.axvline(x=x2, color='royalblue', linestyle='--', linewidth=2, label=f'66th Percentile')
plt.xlabel('Fraction of expressed tissues', fontsize=14, fontweight='bold')
plt.ylabel('# genes', fontsize=14, fontweight='bold')
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(fontsize=10, loc='upper left')
plt.tight_layout()
plt.savefig(path_figures + 'empirical_cutoff_justification_histogram_tissue.png', dpi=600, bbox_inches='tight')
plt.show()


# 5. Tau and CV just for tissues
log_tis = np.log2(matrix + 1)
x_max_tis = np.max(log_tis, axis=1)

with np.errstate(divide='ignore', invalid='ignore'):
    x_norm_tis = log_tis / x_max_tis[:, None]
    tau_tis = np.sum(1 - x_norm_tis, axis=1) / (log_tis.shape[1] - 1)


log_tpm_tis = np.log1p(matrix)
with np.errstate(divide='ignore', invalid='ignore'):
    cv_tis = np.std(log_tpm_tis, axis=1) / np.mean(log_tpm_tis, axis=1)
cv_tis[np.isnan(cv_tis)] = np.nan


# Asegurarnos de que los infinitos o divisiones por cero se marquen como NaN
tau_tis[np.isnan(tau_tis)] = np.nan
mascara_validos = ~np.isnan(tau_tis) & ~np.isinf(tau_tis) & ~np.isnan(frac_tissue) & ~np.isinf(frac_tissue)
coef_pearson, p_valor = pearsonr(frac_tissue[mascara_validos], tau_tis[mascara_validos])

print(f"Correlación de Pearson: {coef_pearson:.4f}")
print(f"P-valor: {p_valor:.4e}")


cv_tis[np.isnan(cv_tis)] = np.nan
mascara_validos = ~np.isnan(cv_tis) & ~np.isinf(cv_tis) & ~np.isnan(frac_tissue) & ~np.isinf(frac_tissue)
coef_pearson, p_valor = pearsonr(frac_tissue[mascara_validos], cv_tis[mascara_validos])

print(f"Correlación de Pearson: {coef_pearson:.4f}")
print(f"P-valor: {p_valor:.4e}")


df_tis = pd.DataFrame({
    'Gene': genes_tis,
    'Tau_Tis': tau_tis,
    'CV_Tis': cv_tis
})



#Figure!!

set_U_tis, set_S_tis, set_hS_tis = set(U_tis), set(S_tis), set(hS_tis)

df_tis['Class_Tis'] = df_tis['Gene'].apply(
    lambda x: 'U' if x in set_U_tis else ('S' if x in set_S_tis else ('hS' if x in set_hS_tis else 'Other'))
)

classes_tis = ['U', 'S', 'hS']
colors = ['royalblue', 'mediumturquoise', 'tomato']

data_tau_tis = [df_tis[df_tis['Class_Tis'] == c]['Tau_Tis'].dropna().values for c in classes_tis]
data_cv_tis=[df_tis[df_tis['Class_Tis'] == c]['CV_Tis'].dropna().values for c in classes_tis]

fig, axes = plt.subplots(1, 2, figsize=(10, 4.5), dpi=600)
def style_violin(parts, colors):
    for pc, color in zip(parts['bodies'], colors):
        pc.set_facecolor(color)
        pc.set_edgecolor(None)
        pc.set_alpha(1.0) # Transparencia
    
    # Colorear las líneas (bordes, máximos, mínimos y cuartiles)
    for partname in ('cbars', 'cmins', 'cmaxes', 'cquantiles'):
        if partname in parts:
            vp = parts[partname]
            vp.set_edgecolor('black')
            vp.set_linewidth(0.9)

parts_tau = axes[0].violinplot(
    data_tau_tis, 
    positions=[1, 2, 3], 
    showmeans=False, 
    showmedians=False,
    showextrema=True, 
    quantiles=[[0.25, 0.5, 0.75]] * len(classes_tis)
)
style_violin(parts_tau, colors)
axes[0].set_ylabel(r'Spatial specificity ($\tau$)', fontsize=20, fontweight='bold')
axes[0].set_xticks([1, 2, 3])
axes[0].set_xticklabels(['Ut', 'St', 'hSt']) 

parts_cv = axes[1].violinplot(
    data_cv_tis, 
    positions=[1, 2, 3], 
    showmeans=False, 
    showmedians=False,
    showextrema=True, 
    quantiles=[[0.25, 0.5, 0.75]] * len(classes_tis)
)
style_violin(parts_cv, colors)
axes[1].set_ylabel('Spatial variability (CV)', fontsize=20, fontweight='bold')
axes[1].set_xticks([1, 2, 3])
axes[1].set_xticklabels(['Ut', 'St', 'hSt'])

for ax in axes:
    ax.tick_params(axis='x', labelsize=24)
    for label in ax.get_xticklabels():
        label.set_fontweight('bold')
    
    ax.tick_params(axis='y', labelsize=15)

plt.tight_layout()
plt.savefig(path_figures + 'continuos_metrics_tissue_classes.png', dpi=600, bbox_inches='tight')
plt.show()



#DEVELOPMENT
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

genes_dev_all=genes
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

#tau dev all
# ==============================================================================
# 3. Tau and CV for all the genes in development
# ==============================================================================
log_dev_all = np.log2(mean_bulk_matrix + 1)
x_max_dev_all = np.max(log_dev_all, axis=1)

with np.errstate(divide='ignore', invalid='ignore'):
    x_norm_dev_all = log_dev_all / x_max_dev_all[:, None]
    tau_dev_all = np.sum(1 - x_norm_dev_all, axis=1) / (log_dev_all.shape[1] - 1)
tau_dev_all[np.isnan(tau_dev_all)] = np.nan

log_tpm_dev_all = np.log1p(mean_bulk_matrix)
with np.errstate(divide='ignore', invalid='ignore'):
    cv_dev_all = np.std(log_tpm_dev_all, axis=1) / np.mean(log_tpm_dev_all, axis=1)
cv_dev_all[np.isnan(cv_dev_all)] = np.nan

df_dev = pd.DataFrame({
    'Gene': genes_dev_all,
    'Tau_Dev': tau_dev_all,
    'CV_Dev': cv_dev_all
})


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



#=============================================================================================
#Comparison development and tissue!!!!!
#=============================================================================================
# 1. Common genes between development and tissue
genes_dev = genes_dev_all  
genes_tis = genes_tis       
mat_dev = mean_bulk_matrix 
mat_tis = matrix          

common_genes, dev_idx, tis_idx = np.intersect1d(genes_dev, genes_tis, return_indices=True)

mat_dev_common = mat_dev[dev_idx]
mat_tis_common = mat_tis[tis_idx]
frac_tis_common=frac_tissue[tis_idx]

# 2. TAU and CV for development (common genes)
log_dev = np.log2(mat_dev_common + 1)
x_max_dev = np.max(log_dev, axis=1)
with np.errstate(divide='ignore', invalid='ignore'):
    x_norm_dev = log_dev / x_max_dev[:, None]
    tau_dev = np.sum(1 - x_norm_dev, axis=1) / (log_dev.shape[1] - 1)
tau_dev[np.isnan(tau_dev)] = np.nan

log_tpm_dev = np.log1p(mat_dev_common)
with np.errstate(divide='ignore', invalid='ignore'):
    cv_dev = np.std(log_tpm_dev, axis=1) / np.mean(log_tpm_dev, axis=1)
cv_dev[np.isnan(cv_dev)] = np.nan

# 3. TAU and CV for tissues (common genes)
log_tis = np.log2(mat_tis_common + 1)
x_max_tis = np.max(log_tis, axis=1)

with np.errstate(divide='ignore', invalid='ignore'):
    x_norm_tis = log_tis / x_max_tis[:, None]
    tau_tis = np.sum(1 - x_norm_tis, axis=1) / (log_tis.shape[1] - 1)


log_tpm_tis = np.log1p(mat_tis_common)
with np.errstate(divide='ignore', invalid='ignore'):
    cv_tis = np.std(log_tpm_tis, axis=1) / np.mean(log_tpm_tis, axis=1)
cv_tis[np.isnan(cv_tis)] = np.nan

# 4. DataFrame 
df_metrics = pd.DataFrame({
    'Gene': common_genes,
    'Tau_Dev': tau_dev,
    'CV_Dev': cv_dev,
    'Tau_Tis': tau_tis,
    'CV_Tis': cv_tis
})


mascara_validos = ~np.isnan(tau_tis) & ~np.isinf(tau_tis) & ~np.isnan(tau_dev) & ~np.isinf(tau_dev)
coef_pearson, p_valor = pearsonr(tau_tis[mascara_validos], tau_dev[mascara_validos])

print(f"Correlación de Pearson: {coef_pearson:.4f}")
print(f"P-valor: {p_valor:.4e}")



# 5. Map original classifications (U, S, U_tis, S_tis)
set_U_dev, set_S_dev, set_hS_dev = set(U_id), set(S_id), set(hS_id)

set_U_tis, set_S_tis, set_hS_tis = set(U_tis), set(S_tis), set(hS_tis)

df_metrics['Class_Tis'] = df_metrics['Gene'].apply(
    lambda x: 'Ut' if x in set_U_tis else ('St' if x in set_S_tis else ('hSt' if x in set_hS_tis else 'Other'))
)

df_metrics['Class_Dev'] = df_metrics['Gene'].apply(
    lambda x: 'U' if x in set_U_dev else ('S' if x in set_S_dev else ('hS' if x in set_hS_dev else 'Other'))
)

classes_tis = ['U_tis', 'S_tis', 'hS_tis'] 
colors = ['royalblue', 'mediumturquoise', 'tomato']

data_tau_tis = [df_metrics[df_metrics['Class_Tis'] == c]['Tau_Tis'].dropna().values for c in classes_tis]
data_cv_tis = [df_metrics[df_metrics['Class_Tis'] == c]['CV_Tis'].dropna().values for c in classes_tis]

classes_dev = ['hS', 'S', 'U']
colors = ['tomato', 'mediumturquoise', 'royalblue']

data_tau_dev = [df_metrics[df_metrics['Class_Dev'] == c]['Tau_Dev'].dropna().values for c in classes_dev]
data_cv_dev = [df_metrics[df_metrics['Class_Dev'] == c]['CV_Dev'].dropna().values for c in classes_dev]


#6. Figure GENEREAL COMPARISON: TAU DEV VS TAU TISSUE
x = tau_dev[mascara_validos] 
y = tau_tis[mascara_validos] 

# Relation 1:1 (y = x)
y_predicha = x 
residuales = y - y_predicha

# top 5% deviations respect diagonal
umbral = np.percentile(np.abs(residuales), 95)
mascara_outliers = np.abs(residuales) > umbral

genes_desv = common_genes[mascara_validos][mascara_outliers]

# ==========================================
#figure deviations
# ==========================================
colores_custom = [
    "deepskyblue", "blue", "darkblue", # Azules (U)
    "darkgreen", "lightgreen", "forestgreen", # Violetas/Fucsia (S)
    "darkred", "firebrick", "lightcoral"  # Rosas (hS)
]

orden_combinaciones = [
    'U-Ut', 'U-St', 'U-hSt',
    'S-U_tis', 'S-S_tis', 'S-hSt',
    'hS-Ut', 'hS-St', 'hS-hSt'
]
paleta_outliers = dict(zip(orden_combinaciones, colores_custom))

df_outliers = df_metrics[df_metrics['Gene'].isin(genes_desv)].copy()
df_outliers['Combinacion'] = df_outliers['Class_Dev'] + "-" + df_outliers['Class_Tis']


plt.figure(figsize=(12.5, 7), dpi=300)
plt.scatter(x, y, alpha=0.3, color='lightgray', s=90)
plt.plot(x, y_predicha, color='blue', linewidth=1, linestyle='--', label='y=x')
sns.scatterplot(
    data=df_outliers, 
    x='Tau_Dev', 
    y='Tau_Tis', 
    hue='Combinacion',
    hue_order=orden_combinaciones,
    palette=paleta_outliers,
    s=90,              
    alpha=0.7,
    edgecolor='white', 
    zorder=3           
)
plt.xlabel(r'Temporal $\tau$', fontsize=30, fontweight='bold')
plt.ylabel(r'Spatial $\tau$', fontsize=30, fontweight='bold')
plt.yticks(fontsize=22)
plt.xticks(fontsize=22)
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=22)
plt.tight_layout()
plt.savefig(path_figures+'comparison_tau.png', dpi=600, bbox_inches='tight')
plt.show()

print(df_outliers['Combinacion'].value_counts())
print(f"\nGenes with highest deviations: {len(genes_desv)}")

conteos = df_outliers['Combinacion'].value_counts()
porcentajes = df_outliers['Combinacion'].value_counts()/len(genes_desv)*100

tabla_resumen = pd.DataFrame({
    'Conteo': conteos,
    'Porcentaje (%)': porcentajes.round(1) # Redondeamos a 1 decimal
})

print("-" * 40)
print(tabla_resumen)


# ==========================================
# 6.1. Comparison: RESIDUAL DISTRIBUTION
# ==========================================
media = np.mean(residuales)
sigma = np.std(residuales)
total_genes = len(residuales)

mask_1sigma = np.abs(residuales - media) <= sigma
mask_2sigma = (np.abs(residuales - media) > sigma) & (np.abs(residuales - media) <= 2 * sigma)
mask_out_2sigma = np.abs(residuales - media) > 2 * sigma

pct_1sigma = (np.sum(mask_1sigma) / total_genes) * 100
pct_2sigma = (np.sum(mask_2sigma) / total_genes) * 100
pct_out_2sigma = (np.sum(mask_out_2sigma) / total_genes) * 100

# ==========================================
# FIGURE
# ==========================================
plt.figure(figsize=(4, 3), dpi=600) 
plt.hist(residuales, bins=100, color='orangered', alpha=0.7)
plt.axvline(media, color='black', linestyle='-', linewidth=1, label=f'$\mu$ ({media:.2f})')
plt.axvline(media + sigma, color='blue', linestyle='--', linewidth=1.2, label=r'$\pm 1\sigma$')
plt.axvline(media - sigma, color='blue', linestyle='--', linewidth=1.2)
plt.axvline(media + 2*sigma, color='green', linestyle=':', linewidth=1.5, label=r'$\pm 2\sigma$')
plt.axvline(media - 2*sigma, color='green', linestyle=':', linewidth=1.5)
plt.xlabel("Residuals\n"r"(Temporal $\tau$ - Spatial $\tau$ )", fontsize=12, fontweight='bold')
plt.ylabel("# genes", fontsize=13, fontweight='bold')
plt.legend(fontsize=8, loc='upper left')
plt.tight_layout()
plt.savefig(path_figures+'residuals_diagonal_tau.png', dpi=600, bbox_inches='tight')
plt.show()

# ==========================================
# statistics
# ==========================================
print("-" * 50)
print("DISTRIBUCIÓN DE GENES SEGÚN DESVIACIÓN (\u03c3)")
print("-" * 50)
print(f"Total de genes analizados: {total_genes}")
print(f"Media (\u03bc): {media:.4f}")
print(f"Desviación Estándar (\u03c3): {sigma:.4f}\n")

print(f"[Zona Normal] Dentro de 1\u03c3:      {np.sum(mask_1sigma)} genes ({pct_1sigma:.1f}%)")
print(f"[Zona Límite] Entre 1\u03c3 y 2\u03c3:   {np.sum(mask_2sigma)} genes ({pct_2sigma:.1f}%)")
print(f"[Outliers]    Fuera de 2\u03c3 (>2\u03c3): {np.sum(mask_out_2sigma)} genes ({pct_out_2sigma:.1f}%)")
print("-" * 50)


genes_outliers_2sigma = common_genes[mascara_validos][mask_out_2sigma]

df_out_2sigma = df_metrics[df_metrics['Gene'].isin(genes_outliers_2sigma)].copy()
if not df_out_2sigma.empty:
    df_out_2sigma['Combinacion'] = df_out_2sigma['Class_Dev'] + "-" + df_out_2sigma['Class_Tis']
    
    print("\nComposición de los Outliers severos (> 2\u03c3):")
    conteos_2sigma = df_out_2sigma['Combinacion'].value_counts()
    porcentajes_2sigma = (conteos_2sigma / len(genes_outliers_2sigma)) * 100
    
    tabla_2sigma = pd.DataFrame({
        'Conteo': conteos_2sigma,
        'Porcentaje (%)': porcentajes_2sigma.round(1)
    })
    print(tabla_2sigma)
    

    
#=============================================================================================
#COMPARISON TAU TISSUE AND AGE!!!!!!
#=============================================================================================
# 1. Common genes between development and tissue

df_gene_origin = pd.read_csv(path + "gene_origin.csv")
genes_df = np.array(list(df_gene_origin['ensembl_gene_id']))
age_df = np.array(list(df_gene_origin['gene_age']))

print("\n" + "="*50)
print(" CORR: TAU TEMP VS GENE ORIGIN")
print("="*50)

# Creamos un diccionario rápido de Gen -> Edad
age_dict = dict(zip(df_gene_origin['ensembl_gene_id'], df_gene_origin['gene_age']))

common_genes, ind_age, ind_dev=np.intersect1d(df_gene_origin['ensembl_gene_id'], genes_dev_all, return_indices=True)
age_compare=list(df_gene_origin['gene_age'][ind_age])
ages_class = np.array(['4290' if x == '>4290' else x for x in age_compare], dtype=int)
tau_dev_compare=np.array(tau_dev_all[ind_dev])

# Calculamos la correlación de Spearman (mejor para edades evolutivas/filoestratigrafía)
corr, p_val = pearsonr(ages_class, tau_dev_compare)

print(f"# common genes: {len(ages_class)}")
print(f"Corr Pearson (rho): {corr:.4f}")
print(f"P-value: {p_val:.2e}")



#sliding window
sw=400
#genes
tau_sorted=np.sort(tau_dev_compare)
index_sort=np.argsort(tau_dev_compare)
age_per_gen_sorted=np.zeros(len(tau_sorted))
for i in range(len(index_sort)):
    age_per_gen_sorted[i]=ages_class[int(index_sort[i])]
    

serie1 = pd.Series(tau_sorted)
serie2 = pd.Series(age_per_gen_sorted)

sw_tau = serie1.rolling(window=sw, center=False).mean()
sw_age = serie2.rolling(window=sw, center=False).mean()

# # figure
plt.figure(figsize=(4, 3), dpi=600)
plt.scatter(sw_tau, sw_age, s=0.8, color='violet')
plt.xlabel(r'Temporal $\tau$', fontsize=13, fontweight='bold')
plt.ylabel('Age (mry)', fontsize=13, fontweight='bold')
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)
plt.legend(markerscale=5)
plt.savefig(path_figures+'age_vs_tau_sw_%d.png' %sw, dpi=600, bbox_inches='tight')
plt.show()



#=================================================================================================
#GO FUNCTIONS FOR TAU!!!!!!!
#=================================================================================================

def enrichment_go_ks_continuous(
    big_matrix: np.ndarray, 
    go_descrip_list: list, 
    genes_matrix: np.ndarray, 
    genes_tau: np.ndarray, 
    tau_values: np.ndarray, 
    min_genes: int = 10, 
    p_threshold: float = 0.001
) -> pd.DataFrame:
    """
    Performs GO enrichment analysis on a continuous variable (e.g., Tau)
    using the Two-Sample Kolmogorov-Smirnov test.
    
    Parameters:
    -----------
    big_matrix : np.ndarray
        Binary association matrix (Genes x GO_Terms). 1 if gene belongs to GO term, 0 otherwise.
    go_descrip_list : list
        List of descriptions for each GO term (corresponds to columns of big_matrix).
    genes_matrix : np.ndarray
        Array of gene IDs corresponding to the rows of big_matrix.
    genes_tau : np.ndarray
        Array of gene IDs corresponding to the tau_values.
    tau_values : np.ndarray
        Continuous metric values (e.g., tau_dev) corresponding to genes_tau.
    min_genes : int
        Minimum number of genes required in a GO term to perform the test (reduces noise).
    p_threshold : float
        P-value threshold for significance.
        
    Returns:
    --------
    pd.DataFrame
        DataFrame containing significant GO terms sorted by P-value.
    """
    # 1. Intersect genes from the GO matrix with genes that have valid Tau values
    common_genes, ind_mat, ind_tau = np.intersect1d(
        genes_matrix, genes_tau, return_indices=True
    )
    
    matrix_matched = big_matrix[ind_mat, :]
    tau_matched = tau_values[ind_tau]
    
    resultados = []
    
    # 2. Iterate through each GO term (columns in the matrix)
    for fen in range(matrix_matched.shape[1]):
        # Boolean mask: genes inside this specific GO term
        in_term = matrix_matched[:, fen] == 1
        n_in = np.sum(in_term)
        
        # Ignore terms with too few genes to avoid false positives from noise
        if n_in < min_genes:
            continue
            
        # Tau values of genes IN the term vs OUTSIDE the term
        tau_in = tau_matched[in_term]
        tau_out = tau_matched[~in_term]
        
        # 3. Kolmogorov-Smirnov Test
        stat, p_val = ks_2samp(tau_in, tau_out)
        
        if p_val < p_threshold:
            # Calculate shift direction using medians
            diff_median = np.median(tau_in) - np.median(tau_out)
            
            if diff_median < 0:
                direccion = "Ubiquitous (Low Tau)"
            else:
                direccion = "Specific (High Tau)"
                
            resultados.append({
                'GO_Term': go_descrip_list[fen],
                'P_value': p_val,
                'KS_Stat': stat,
                'Diff_Median': diff_median,
                'Direction': direccion,
                'N_Genes': n_in
            })
            
    # 4. Convert to DataFrame and sort by significance
    df_res = pd.DataFrame(resultados)
    if not df_res.empty:
        df_res = df_res.sort_values('P_value').reset_index(drop=True)
        
    return df_res

def print_top_ks_results(df_res: pd.DataFrame, ontology_name: str, top_n: int = 5):
    """
    Prints a formatted summary of the top enriched GO terms separated by direction.
    """
    print(f"\n---> {ontology_name.upper()} <---")
    
    if df_res.empty:
        print("  No terms passed the significance filter.")
        return
        
    # Separate into ubiquitous vs specific tendencies
    df_ubicuo = df_res[df_res['Direction'] == "Ubiquitous (Low Tau)"]
    df_especifico = df_res[df_res['Direction'] == "Specific (High Tau)"]
    
    print(f"  [Enriched towards UBIQUITOUS (Low \u03c4)] - Total: {len(df_ubicuo)}")
    for _, row in df_ubicuo.head(top_n).iterrows():
        print(f"    - {row['GO_Term']} (n={row['N_Genes']}, p={row['P_value']:.2e})")
        
    print(f"\n  [Enriched towards SPECIFIC (High \u03c4)] - Total: {len(df_especifico)}")
    for _, row in df_especifico.head(top_n).iterrows():
        print(f"    - {row['GO_Term']} (n={row['N_Genes']}, p={row['P_value']:.2e})")


# ===================================================================================
# EXECUTE: CONTINUOUS GO ENRICHMENT (KS TEST ON TAU_DEV)
# ===================================================================================
print("\n" + "="*60)
print(" CONTINUOUS GO ENRICHMENT (KS TEST ON TAU)")
print("="*60)

# Assuming common_genes, mascara_validos, and tau_dev are defined earlier in your script
genes_validos_tau = genes_dev_all
valores_tau_dev = tau_dev_all


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

big_matrix_bio_process = create_association_matrix_gene_go_term(genes_id_type, GO_bio_process, pathlist_bio_process_unique, GO_bio_process_descrip)
big_matrix_molecular = create_association_matrix_gene_go_term(genes_id_type, GO_molecular, pathlist_molecular_unique, GO_molecular_descrip)
big_matrix_cell_comp = create_association_matrix_gene_go_term(genes_id_type, GO_cell_comp, pathlist_cell_comp_unique, GO_cell_comp_descrip)



# 1. Biological Process
print("\nCalculating KS for Biological Process...")
df_ks_bp = enrichment_go_ks_continuous(
    big_matrix=big_matrix_bio_process, 
    go_descrip_list=GO_bio_process_descrip, 
    genes_matrix=genes_id_type, 
    genes_tau=genes_validos_tau, 
    tau_values=valores_tau_dev
)

# 2. Molecular Function
print("Calculating KS for Molecular Function...")
df_ks_mf = enrichment_go_ks_continuous(
    big_matrix=big_matrix_molecular, 
    go_descrip_list=GO_molecular_descrip, 
    genes_matrix=genes_id_type, 
    genes_tau=genes_validos_tau, 
    tau_values=valores_tau_dev
)

# 3. Cellular Component
print("Calculating KS for Cellular Component...")
df_ks_cc = enrichment_go_ks_continuous(
    big_matrix=big_matrix_cell_comp, 
    go_descrip_list=GO_cell_comp_descrip, 
    genes_matrix=genes_id_type, 
    genes_tau=genes_validos_tau, 
    tau_values=valores_tau_dev
)

# Print Summaries
print_top_ks_results(df_ks_bp, "Biological Process", top_n=5)
print_top_ks_results(df_ks_mf, "Molecular Function", top_n=5)
print_top_ks_results(df_ks_cc, "Cellular Component", top_n=5)



# 1. Add an 'Ontology' label to track where each term came from
# (Assuming df_ks_bp, df_ks_mf, and df_ks_cc have already been generated)
if not df_ks_bp.empty: df_ks_bp['Ontology'] = 'Biological Process'
if not df_ks_mf.empty: df_ks_mf['Ontology'] = 'Molecular Function'
if not df_ks_cc.empty: df_ks_cc['Ontology'] = 'Cellular Component'

# 2. Combine all dataframes into a single unified DataFrame
list_of_dfs = [df for df in [df_ks_bp, df_ks_mf, df_ks_cc] if not df.empty]
df_ks_all = pd.concat(list_of_dfs, ignore_index=True)

# 3. Sort the unified dataframe by P-value
df_ks_all = df_ks_all.sort_values(by='P_value').reset_index(drop=True)

# 4. Separate into High Tau (Specific) and Low Tau (Ubiquitous) dataframes
df_high_tau = df_ks_all[df_ks_all['Direction'] == "Specific (High Tau)"].reset_index(drop=True)
df_low_tau  = df_ks_all[df_ks_all['Direction'] == "Ubiquitous (Low Tau)"].reset_index(drop=True)

# ==========================================
# VERIFICATION PRINTS
# ==========================================
print(f"Total Significant GO Terms: {len(df_ks_all)}")
print(f"Total High Tau (Specific) Terms: {len(df_high_tau)}")
print(f"Total Low Tau (Ubiquitous) Terms: {len(df_low_tau)}\n")

print("--> TOP 5 HIGH TAU (SPECIFIC) TERMS <--")
print(df_high_tau[['Ontology', 'GO_Term', 'P_value', 'N_Genes']].head())

print("\n--> TOP 5 LOW TAU (UBIQUITOUS) TERMS <--")
print(df_low_tau[['Ontology', 'GO_Term', 'P_value', 'N_Genes']].head())

nombre_archivo_excel = path_figures + 'Tau_GO_Enrichment.xlsx'

with pd.ExcelWriter(nombre_archivo_excel, engine='openpyxl') as writer:
    df_ks_all.to_excel(writer, sheet_name='Todos_los_Resultados', index=False)
    
    df_high_tau.to_excel(writer, sheet_name='High_Tau_Específicos', index=False)
    
    df_low_tau.to_excel(writer, sheet_name='Low_Tau_Ubicuos', index=False)
