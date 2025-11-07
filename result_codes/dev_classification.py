# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 15:59:37 2025

Bulk zebrafish

@author: Alicia
"""
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import warnings; warnings.simplefilter('ignore')
from itertools import product
from scipy import stats


path_save_data='PATH_WHERE_YOU_KEEP_THE_REQUIRED_DATASETS'


#read the data
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


#1.) We plot the distribution of the #genes expressed in each stage.
#We take all the embryos in each stage
#In each stage, i count the number of active genes (in each embryo)
mean_per_stage=np.zeros(len(bulk_matrix[0, :]))
act_genes_per_stage=np.zeros(len(bulk_matrix[0, :]))
for i in range(len(bulk_matrix[0, :])):
    non_null=len(np.where(bulk_matrix[:, i]>0)[0])
    mean_per_stage[i]=np.mean(bulk_matrix[:, i])
    act_genes_per_stage[i]=non_null/len(bulk_matrix[:, 0])

data_reshaped = act_genes_per_stage.reshape(18,5)
data_reshaped=data_reshaped.transpose()

#1.1.) Figure of gene divesity per stage
plt.figure(figsize=(11, 6), dpi=600)
sns.boxplot(data=data_reshaped)
plt.xticks(ticks=np.arange(len(time_unique)), labels=time_unique, rotation=90, fontsize=25)
plt.yticks(fontsize=25)
plt.xlabel('Stages', fontsize=30, fontweight='bold')
plt.ylabel('% different genes', fontsize=26, fontweight='bold')
plt.savefig(path_save_data+'frac_diff_genes_per_stage.png', dpi=600, bbox_inches='tight')
plt.show()


data_reshaped = mean_per_stage.reshape(18,5)
data_reshaped=data_reshaped.transpose()

stage_1=data_reshaped[:, 14]
stage_2=data_reshaped[:, 16]
# Perform the paired t-test
t_stat, p_value = stats.ttest_rel(stage_1, stage_2)

# Output the result
print(f"T-statistic: {t_stat}")
print(f"P-value: {p_value}")

stage_1=data_reshaped[:, 0]
stage_2=data_reshaped[:, 2]
# Perform the paired t-test
t_stat, p_value = stats.ttest_rel(stage_1, stage_2)

# Output the result
print(f"T-statistic: {t_stat}")
print(f"P-value: {p_value}")

#1.2.) Figure of average TPM (gene expression) per stage
plt.figure(figsize=(11, 6), dpi=600)
sns.boxplot(data=data_reshaped)
plt.xticks(ticks=np.arange(len(time_unique)), labels=time_unique, rotation=90, fontsize=25)
plt.yticks(fontsize=25)
plt.xlabel('Stages', fontsize=30, fontweight='bold')
plt.ylabel('<TPM>', fontsize=26, fontweight='bold')
plt.savefig(path_save_data+'average_TPM_per_stage.png', dpi=600, bbox_inches='tight')
plt.show()

# #2.) Thereshold: 1 TPM
# bulk_matrix[bulk_matrix<1]=0


#3.) We create a mean_bulk_matrix averaging the embryos in the same times 
mean_bulk_matrix=np.zeros((len(genes), len(time_unique)))
for k in range(len(genes)):
    count=0
    for i in range(len(time_unique)):
        sum_embryo=0
        for j in range(len(embryo)):
            sum_embryo=sum_embryo+bulk_matrix[k][int(count)]
            count=count+1
        mean_bulk_matrix[k][i]=sum_embryo/len(embryo)

total_tpm_per_stage=np.sum(mean_bulk_matrix, axis=0)


#3.1.) TPM mean expression
mean_bulk_matrix_tpm=np.zeros((len(genes), len(time_unique)))
#We TPM the mean_bulk_matrix
for i in range(len(time_unique)):
    for j in range(len(genes)):
        mean_bulk_matrix_tpm[j][i]=mean_bulk_matrix[j][i]*1000000/(np.sum(mean_bulk_matrix[:, i]))

total_tpm_per_stage2=np.sum(mean_bulk_matrix_tpm, axis=0)
mean_bulk_matrix=mean_bulk_matrix_tpm
del mean_bulk_matrix_tpm


# 3.2.) Thereshold: 1 TPM
mean_bulk_matrix_clean = np.where(mean_bulk_matrix < 1, 0, mean_bulk_matrix)

np.savetxt(path_save_data+'mean_bulk_matrix.txt', mean_bulk_matrix, fmt='%s', delimiter=',')
np.savetxt(path_save_data+'mean_bulk_matrix_clean.txt', mean_bulk_matrix_clean, fmt='%s', delimiter=',')


#3.1.) We plot the histogram of the total mean gene counts per gene in each stage
plt.figure(figsize=(4, 3),dpi=600)
plt.hist(mean_bulk_matrix.flatten(), bins=100, color='lightseagreen', log=True)
plt.xlabel('TPM per gene and stage', fontsize=14, fontweight='bold')
plt.ylabel('freq', fontsize=14, fontweight='bold')
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.grid(False)
# plt.axvline(x=np.median(genes_n_cell), color='darkslategrey', linestyle='--', lw=1)
plt.savefig(path_save_data+'total_mean_gene_counts_per_gene_stage.png', dpi=600, bbox_inches='tight')
plt.show()  

#3.2.) We plot the distribution of stages in which genes are active
act_stage_per_gene=np.zeros(len(genes))
for i in range(len(genes)):
    non_null=len(np.where(mean_bulk_matrix_clean[i, :]>0)[0])
    act_stage_per_gene[i]=non_null

null_genes=len(np.where(act_stage_per_gene==0)[0])

n_stages, n_times=np.unique(act_stage_per_gene, return_counts=True)
n_stages=np.array(n_stages, dtype=int)
plt.figure(figsize=(11, 6), dpi=600)
plt.bar(n_stages, n_times, color='darkorange', edgecolor='black')
plt.xticks(ticks=np.arange(0, len(n_stages)), labels=n_stages, rotation=90, fontsize=25)
plt.yticks(fontsize=25)
plt.xlabel('# stages expressed', fontsize=30, fontweight='bold')
plt.ylabel('# genes (TPM>1)', fontsize=30, fontweight='bold')
plt.savefig(path_save_data+'act_stages_per_gene.png', dpi=600, bbox_inches='tight')
plt.show()



#4.) We divided the genes in specific, high-specific and ubiqutouly expressed
specific=[]
high_specific=[]
ubiq_genes=[]
null_genes=[]
null_genes_id=[]
hS_id=[]
U_id=[]
S_id=[]
act_stages_S=[]
act_stages_hS=[]
act_stages_U=[]
for i in range(len(genes)):
    if act_stage_per_gene[i]==0:
        null_genes.append(genes_name[i])
        null_genes_id.append(genes[i])
    if (act_stage_per_gene[i]>0) & (act_stage_per_gene[i]<=8):
        high_specific.append(genes_name[i])
        hS_id.append(genes[i])
        act_stages_hS.append(act_stage_per_gene[i])
    if (act_stage_per_gene[i]>8) & (act_stage_per_gene[i]<=17):
        specific.append(genes_name[i])
        S_id.append(genes[i])
        act_stages_S.append(act_stage_per_gene[i])
    if act_stage_per_gene[i]==18:
        ubiq_genes.append(genes_name[i])
        U_id.append(genes[i])
        act_stages_U.append(act_stage_per_gene[i])

        
print(len(U_id)+len(S_id)+len(hS_id)+len(null_genes_id))

np.savetxt(path_save_data+'hS_id_bulk.txt', hS_id, fmt='%s', delimiter=',')
np.savetxt(path_save_data+'S_id_bulk.txt', S_id, fmt='%s', delimiter=',')
np.savetxt(path_save_data+'U_id_bulk.txt', U_id, fmt='%s', delimiter=',')
np.savetxt(path_save_data+'null_genes.txt', null_genes_id, fmt='%s', delimiter=',')

np.median(act_stages_S)
np.median(act_stages_hS)
from matplotlib.ticker import MaxNLocator
data=[act_stages_U, act_stages_S, act_stages_hS]
fig, ax = plt.subplots(figsize=(4,3), dpi=600)
box = ax.boxplot(data, labels=["U", "S", "hS"],  widths=0.45, patch_artist=True, 
                 flierprops=dict(marker='o', color='gray', markersize=0.5), 
                 boxprops=dict(edgecolor="none"))
colors = ["royalblue", "mediumturquoise", "tomato"]
for patch, color in zip(box["boxes"], colors):
    patch.set_facecolor(color)
plt.setp(box["medians"], color="black", linewidth=1)
plt.setp(box["whiskers"], color="black", linestyle="--", linewidth=0.5)
plt.setp(box["caps"], color="black", linewidth=0.5)
plt.ylabel("# expressed stages", fontsize=15, fontweight='bold')
# plt.xlabel("Gene classes", fontsize=15)
plt.xticks(fontsize=18, fontweight='bold')
plt.yticks(fontsize=16)
ax = plt.gca()
ax.yaxis.set_major_locator(MaxNLocator(integer=True))
plt.savefig(path_save_data+'act_stages_per_class_dev.png', dpi=600, bbox_inches='tight')
plt.show()


#===================================================================================
#4.0.) Classification without the threshold

#4.0.1.) We plot the distribution of stages in which genes are active
act_stage_per_gene=np.zeros(len(genes))
for i in range(len(genes)):
    non_null=len(np.where(mean_bulk_matrix[i, :]>0)[0])
    act_stage_per_gene[i]=non_null

null_genes=len(np.where(act_stage_per_gene==0)[0])

n_stages, n_times=np.unique(act_stage_per_gene, return_counts=True)
n_stages=np.array(n_stages, dtype=int)
plt.figure(figsize=(11, 6), dpi=600)
plt.bar(n_stages, n_times, color='darkorange', edgecolor='black')
plt.xticks(ticks=np.arange(0, len(n_stages)), labels=n_stages, rotation=90, fontsize=25)
plt.yticks(fontsize=25)
plt.xlabel('# stages expressed', fontsize=30, fontweight='bold')
plt.ylabel('# genes (TPM>0)', fontsize=30, fontweight='bold')
plt.savefig(path_save_data+'act_stages_per_gene_without_thres.png', dpi=600, bbox_inches='tight')
plt.show()

#4.0.2.) we search for the classification without the threshold
null_genes_id_no_thres=[]
hS_id_no_thres=[]
U_id_no_thres=[]
S_id_no_thres=[]
for i in range(len(genes)):
    if act_stage_per_gene[i]==0:
        null_genes_id_no_thres.append(genes[i])
    if (act_stage_per_gene[i]>0) & (act_stage_per_gene[i]<=8):
        hS_id_no_thres.append(genes[i])
    if (act_stage_per_gene[i]>8) & (act_stage_per_gene[i]<=17):
        S_id_no_thres.append(genes[i])
    if act_stage_per_gene[i]==18:
        U_id_no_thres.append(genes[i])
        
print(len(U_id_no_thres)+len(S_id_no_thres)+len(hS_id_no_thres)+len(null_genes_id_no_thres))

U_no_thres_U=np.intersect1d(U_id_no_thres, U_id)
U_no_thres_S=np.intersect1d(U_id_no_thres, S_id)
U_no_thres_hS=np.intersect1d(U_id_no_thres, hS_id)
U_no_thres_null_genes=np.intersect1d(U_id_no_thres, null_genes_id)

frac_genes_U=np.array([len(U_no_thres_U)/len(U_id_no_thres), 
                     len(U_no_thres_S)/len(U_id_no_thres), 
                     len(U_no_thres_hS)/len(U_id_no_thres), 
                     len(U_no_thres_null_genes)/len(U_id_no_thres)])

S_no_thres_U=np.intersect1d(S_id_no_thres, U_id)
S_no_thres_S=np.intersect1d(S_id_no_thres, S_id)
S_no_thres_hS=np.intersect1d(S_id_no_thres, hS_id)
S_no_thres_null_genes=np.intersect1d(S_id_no_thres, null_genes_id)

frac_genes_S=np.array([len(S_no_thres_U)/len(S_id_no_thres), 
                     len(S_no_thres_S)/len(S_id_no_thres), 
                     len(S_no_thres_hS)/len(S_id_no_thres), 
                     len(S_no_thres_null_genes)/len(S_id_no_thres)])

hS_no_thres_U=np.intersect1d(hS_id_no_thres, U_id)
hS_no_thres_S=np.intersect1d(hS_id_no_thres, S_id)
hS_no_thres_hS=np.intersect1d(hS_id_no_thres, hS_id)
hS_no_thres_null_genes=np.intersect1d(hS_id_no_thres, null_genes_id)

frac_genes_hS=np.array([len(hS_no_thres_U)/len(hS_id_no_thres), 
                     len(hS_no_thres_S)/len(hS_id_no_thres), 
                     len(hS_no_thres_hS)/len(hS_id_no_thres), 
                     len(hS_no_thres_null_genes)/len(hS_id_no_thres)])


gene_names = ["U", "S", "hS", 'NE']
categories = ["hS", "S", "U"]
frac_genes = np.array([frac_genes_hS, frac_genes_S, frac_genes_U])
n_genes_plot = np.array([frac_genes_hS*len(hS_id_no_thres), frac_genes_S*len(S_id_no_thres), frac_genes_U*len(U_id_no_thres)])

color=['royalblue', 'mediumturquoise', 'tomato', 'gray', 
       'royalblue', 'mediumturquoise', 'tomato', 'gray', 
       'royalblue', 'mediumturquoise', 'tomato', 'gray']

fig, ax = plt.subplots(figsize=(5, 4))

count=0
for i, category in enumerate(categories):
    for j, gene in enumerate(gene_names):
        # plt.scatter(j+0.5, i+0.5, s=frac_genes[i][j] * 1000, color=color[count])
        plt.text(j + 0.5, i + 0.5, str(int(n_genes_plot[i][j])), color=color[count],
                 fontsize=15, ha='center', va='center', fontweight='bold')
        count=count+1

# Add dashed guide lines at midpoints
for i in range(len(categories)-1):
    plt.axhline(y=i+1, color='gray', linestyle='--', alpha=0.5, lw=1)
for j in range(len(gene_names)-1):
    plt.axvline(x=j+1, color='gray', linestyle='--', alpha=0.5, lw=1)
ax.set_xticks(np.arange(len(gene_names)) + 0.5)
ax.set_xticklabels(gene_names, fontsize=20)
ax.set_yticks(np.arange(len(categories)) + 0.5)
ax.set_yticklabels(categories, fontsize=20)

# Add row totals on the right
row_totals = [len(hS_id_no_thres), len(S_id_no_thres), len(U_id_no_thres)]
for i, total in enumerate(row_totals):
    plt.text(len(gene_names) + 0.5, i + 0.5, str(int(total)),
             va='center', ha='center', fontsize=15, fontweight='bold')

col_totals = [len(U_id), len(S_id), len(hS_id), len(null_genes_id)]  # Asegúrate de que el orden coincida con gene_names
for j, total in enumerate(col_totals):
    plt.text(j + 0.5, -0.2, str(int(total)),  # y = 0 lo coloca debajo
             va='center', ha='center', fontsize=15, fontweight='bold')

plt.ylabel("TPM > 0", fontsize=19, fontweight='bold')
plt.xlim(0, 4)
plt.ylim(0, 3)
plt.xlabel("TPM > 1", fontsize=19, fontweight='bold')
ax.xaxis.set_label_position('top') 
ax.xaxis.tick_top()
plt.savefig(path_save_data+'threshold_comparison.png', dpi=600, bbox_inches='tight')
plt.show()
#===================================================================================






frac_genes=np.array([len(null_genes_id)/len(genes), len(ubiq_genes)/len(genes), len(specific)/len(genes), len(high_specific)/len(genes)])
groups=['NE', 'U', 'S', 'hS']
color=['gray', 'royalblue', 'mediumturquoise', 'tomato']

plt.figure(figsize=(4, 3), dpi=600)
plt.bar(groups, frac_genes, color=color)
plt.xticks(ticks=np.arange(len(groups)), labels=groups, fontsize=20, fontweight='bold')
plt.yticks(fontsize=15)
plt.ylabel('Gene fraction', fontsize=20, fontweight='bold')
plt.savefig(path_save_data+'gene_groups.png', dpi=600, bbox_inches='tight')
plt.show()


#4.1.) We plot the expression of two of them in each stage
gene_list=['rpl24', 'cnppd1', 'cyyr1', 'mbnl2', 'urp2', 'adamts5']
label_title=['U', 'U', 'S', 'S', 'hS', 'hS']
color=['royalblue', 'royalblue', 'mediumturquoise', 'mediumturquoise','tomato', 'tomato']
for i in range(len(gene_list)):
    ind=np.where(genes_name==gene_list[i])[0]
    plt.figure(figsize=(6, 4), dpi=600)
    plt.plot(time_unique, mean_bulk_matrix[int(ind), :], color=color[i])
    plt.suptitle(label_title[i], fontsize=12)
    plt.title(gene_list[i], fontsize=15, fontweight='bold')
    plt.ylabel('<TPM>', fontsize=20)
    plt.yticks(fontsize=15)
    plt.xticks(ticks=np.arange(len(time_unique)), labels=time_unique, rotation=90, fontsize=13)
    plt.savefig(path_save_data+'gene_exp_%s.png' %gene_list[i], dpi=600, bbox_inches='tight')
    plt.show()
    

#5.) We search for the highest expressed stage in each gene
U_id=np.array(U_id, dtype=str)
S_id=np.array(S_id, dtype=str)
hS_id=np.array(hS_id, dtype=str)

high_stage_U=np.zeros(len(time_unique))
for i in range(len(ubiq_genes)):
    ind=np.where(genes==U_id[i])[0]
    max_stage=np.argmax(mean_bulk_matrix[int(ind), :])
    high_stage_U[int(max_stage)]=high_stage_U[int(max_stage)]+1
    
plt.figure(figsize=(11, 6), dpi=600)
plt.bar(time_unique, high_stage_U/len(U_id)*100, color='royalblue')
plt.xticks(ticks=np.arange(0, len(time_unique)), labels=time_unique, rotation=90, fontsize=25)
plt.yticks(fontsize=25)
plt.title('Ubiquitously expressed genes', fontweight='bold', fontsize=30)
plt.xlabel('Stages', fontsize=30, fontweight='bold')
plt.ylabel('% genes with max expression', fontsize=25)
plt.savefig(path_save_data+'max_expressed_stages_U_genes.png', dpi=600, bbox_inches='tight')
plt.show()
    
high_stage_S=np.zeros(len(time_unique))
for i in range(len(S_id)):
    ind=np.where(genes==S_id[i])[0]
    max_stage=np.argmax(mean_bulk_matrix[int(ind), :])
    high_stage_S[int(max_stage)]=high_stage_S[int(max_stage)]+1
    
plt.figure(figsize=(11, 6), dpi=600)
plt.bar(time_unique, high_stage_S/len(S_id)*100, color='mediumturquoise')
plt.xticks(ticks=np.arange(0, len(time_unique)), labels=time_unique, rotation=90, fontsize=25)
plt.yticks(fontsize=25)
plt.title('Specific genes', fontweight='bold', fontsize=30)
plt.xlabel('Stages', fontsize=30, fontweight='bold')
plt.ylabel('% genes with max expression', fontsize=25)
plt.savefig(path_save_data+'max_expressed_stages_S_genes.png', dpi=600, bbox_inches='tight')
plt.show()

high_stage_hS=np.zeros(len(time_unique))
for i in range(len(hS_id)):
    ind=np.where(genes==hS_id[i])[0]
    max_stage=np.argmax(mean_bulk_matrix[int(ind), :])
    high_stage_hS[int(max_stage)]=high_stage_hS[int(max_stage)]+1
            
plt.figure(figsize=(11, 6), dpi=600)
plt.bar(time_unique, high_stage_hS/len(hS_id)*100, color='tomato')
plt.xticks(ticks=np.arange(0, len(time_unique)), labels=time_unique, rotation=90, fontsize=25)
plt.yticks(fontsize=25)
plt.title('Highly Specific genes', fontweight='bold', fontsize=30)
plt.xlabel('Stages', fontsize=30, fontweight='bold')
plt.ylabel('% genes with max expression', fontsize=25)
plt.savefig(path_save_data+'max_expressed_stages_hS_genes.png', dpi=600, bbox_inches='tight')
plt.show()

#mixed fig
width = 0.25  # ancho de cada barra
x = np.arange(len(time_unique))
plt.figure(figsize=(11, 6), dpi=600)
plt.bar(x - width, high_stage_U / len(U_id) * 100, width=width, label='U', color='royalblue')
plt.bar(x,         high_stage_S / len(S_id) * 100, width=width, label='S', color='mediumturquoise')
plt.bar(x + width, high_stage_hS / len(hS_id) * 100, width=width, label='hS', color='tomato')

plt.xticks(ticks=x, labels=time_unique, rotation=90, fontsize=25)
plt.yticks(fontsize=25)
plt.xlabel('Stages', fontsize=30)
plt.ylabel('% genes with\nmax expression', fontsize=28, fontweight='bold')
plt.legend(fontsize=20)

plt.tight_layout()
plt.savefig(path_save_data + 'max_expressed_stages_all_classes.png', dpi=600, bbox_inches='tight')
plt.show()


#6.) We find the fraction of genes expressed in each stage
frac_U_genes_per_stage=np.zeros(len(time_unique))
frac_hS_genes_per_stage=np.zeros(len(time_unique))
frac_S_genes_per_stage=np.zeros(len(time_unique))

for i in range(len(time_unique)):
    ind_genes=np.where(mean_bulk_matrix_clean[:, i]>0)[0]
    genes_set=genes[ind_genes]
    frac_U_genes_per_stage[i]=len(np.intersect1d(genes_set, U_id))/len(genes_set)
    frac_hS_genes_per_stage[i]=len(np.intersect1d(genes_set, hS_id))/len(genes_set)
    frac_S_genes_per_stage[i]=len(np.intersect1d(genes_set, S_id))/len(genes_set)    
    
#figure
plt.figure(figsize=(11, 6), dpi=600)
plt.bar(time_unique, frac_U_genes_per_stage, 0.7, label='U', color='royalblue')
plt.bar(time_unique, frac_S_genes_per_stage, 0.7, bottom=frac_U_genes_per_stage, label='S', color='mediumturquoise')
plt.bar(time_unique, frac_hS_genes_per_stage, 0.7, bottom=frac_U_genes_per_stage + frac_S_genes_per_stage, label='hS', color='tomato')
plt.xticks(ticks=np.arange(0, len(time_unique)), labels=time_unique, rotation=90, fontsize=25)
plt.ylabel('% expressed genes', fontsize=25, fontweight='bold')
plt.xlabel('Stages', fontsize=30)
plt.yticks(fontsize=20)
plt.legend(
    fontsize=25,  
    loc='center left', 
    bbox_to_anchor=(1, 0.7)  
)
plt.tight_layout(rect=[0, 0, 1, 1]) 
plt.savefig(path_save_data+'frac_genes_each_group_per_stage.png', dpi=600, bbox_inches='tight')
plt.show()


print('Initial and final gene fraction')
print('hS', frac_hS_genes_per_stage[0]*100, frac_hS_genes_per_stage[int(len(time_unique)-1)]*100)
print('U', frac_U_genes_per_stage[0]*100, frac_U_genes_per_stage[int(len(time_unique)-1)]*100)
print('S', frac_S_genes_per_stage[0]*100, frac_S_genes_per_stage[int(len(time_unique)-1)]*100)

#7.) We find the contribution of each group to the total TPM in each stage
total_TPM_stage=np.sum(mean_bulk_matrix, axis=0)

frac_TPM_U_genes_per_stage=np.zeros(len(time_unique))
frac_TPM_hS_genes_per_stage=np.zeros(len(time_unique))
frac_TPM_S_genes_per_stage=np.zeros(len(time_unique))

for i in range(len(time_unique)):
    ind_genes=np.where(mean_bulk_matrix[:, i]>0)[0]
    genes_set=genes[ind_genes]
    U_genes_per_stage=np.intersect1d(genes_set, U_id)
    for j in range(len(U_genes_per_stage)):
        ind_gene=np.where(genes==U_genes_per_stage[j])[0]
        frac_TPM_U_genes_per_stage[i]=frac_TPM_U_genes_per_stage[i]+mean_bulk_matrix[int(ind_gene), i]
    frac_TPM_U_genes_per_stage[i]=frac_TPM_U_genes_per_stage[i]/total_TPM_stage[i]
    S_genes_per_stage=np.intersect1d(genes_set, S_id)
    for j in range(len(S_genes_per_stage)):
        ind_gene=np.where(genes==S_genes_per_stage[j])[0]
        frac_TPM_S_genes_per_stage[i]=frac_TPM_S_genes_per_stage[i]+mean_bulk_matrix[int(ind_gene), i]
    frac_TPM_S_genes_per_stage[i]=frac_TPM_S_genes_per_stage[i]/total_TPM_stage[i]
    hS_genes_per_stage=np.intersect1d(genes_set, hS_id)
    for j in range(len(hS_genes_per_stage)):
        ind_gene=np.where(genes==hS_genes_per_stage[j])[0]
        frac_TPM_hS_genes_per_stage[i]=frac_TPM_hS_genes_per_stage[i]+mean_bulk_matrix[int(ind_gene), i]
    frac_TPM_hS_genes_per_stage[i]=frac_TPM_hS_genes_per_stage[i]/total_TPM_stage[i]


#figure
plt.figure(figsize=(11, 6), dpi=600)
plt.bar(time_unique, frac_TPM_U_genes_per_stage, 0.7, label='U', color='royalblue')
plt.bar(time_unique, frac_TPM_S_genes_per_stage, 0.7, bottom=frac_TPM_U_genes_per_stage, label='S', color='mediumturquoise')
plt.bar(time_unique, frac_TPM_hS_genes_per_stage, 0.7, bottom=frac_TPM_U_genes_per_stage + frac_TPM_S_genes_per_stage, label='hS', color='tomato')
plt.xticks(ticks=np.arange(0, len(time_unique)), labels=time_unique, rotation=90, fontsize=25)
plt.ylabel('% <TPM>', fontsize=25, fontweight='bold')
plt.xlabel('Stages', fontsize=30)
plt.yticks(fontsize=20)
plt.legend(
    fontsize=25,  
    loc='center left', 
    bbox_to_anchor=(1, 0.7)  
)
plt.tight_layout(rect=[0, 0, 1, 1]) 
plt.savefig(path_save_data+'frac_TPM_each_group_per_stage.png', dpi=600, bbox_inches='tight')
plt.show()


print('Initial and final <TPM> fraction')
print('hS', frac_TPM_hS_genes_per_stage[0]*100, frac_TPM_hS_genes_per_stage[int(len(time_unique)-1)]*100)
print('U', frac_TPM_U_genes_per_stage[0]*100, frac_TPM_U_genes_per_stage[int(len(time_unique)-1)]*100)
print('S', frac_TPM_S_genes_per_stage[0]*100, frac_TPM_S_genes_per_stage[int(len(time_unique)-1)]*100)


#====================================
#Scatter plot (6 vs 7)
#Frac <TPM> vs frac #genes

#U genes
plt.figure(figsize=(4,3))
plt.scatter(frac_U_genes_per_stage, frac_TPM_U_genes_per_stage, color='royalblue', label='U', s=800)
for i in range(len(time_unique)):
    plt.text(frac_U_genes_per_stage[i]-0.017, frac_TPM_U_genes_per_stage[i]-0.005, time_unique[i], fontsize=8, color='white', fontweight='bold')
plt.ylabel('% <TPM>', fontsize=15)
plt.xlabel('% different genes', fontsize=15)
plt.title('U genes', fontsize=15, fontweight='bold')
plt.yticks(fontsize=13)
plt.xticks(fontsize=13)
plt.savefig(path_save_data+'frac_U_genes_vs_frac_tpm_U.png', dpi=600, bbox_inches='tight')
plt.show()

#S genes
plt.figure(figsize=(4,3))
plt.scatter(frac_S_genes_per_stage, frac_TPM_S_genes_per_stage, color='mediumturquoise', label='S', s=800)
for i in range(len(time_unique)):
    plt.text(frac_S_genes_per_stage[i]-0.01, frac_TPM_S_genes_per_stage[i]-0.002, time_unique[i], fontsize=8, color='blue', fontweight='bold')
plt.ylabel('% <TPM>', fontsize=15)
plt.xlabel('% different genes', fontsize=15)
plt.title('S genes', fontsize=15, fontweight='bold')
plt.yticks(fontsize=13)
plt.xticks(fontsize=13)
plt.savefig(path_save_data+'frac_S_genes_vs_frac_tpm_S.png', dpi=600, bbox_inches='tight')
plt.show()

#hS genes
plt.figure(figsize=(4,3))
plt.scatter(frac_hS_genes_per_stage, frac_TPM_hS_genes_per_stage, color='tomato', label='hS', s=800)
for i in range(len(time_unique)):
    plt.text(frac_hS_genes_per_stage[i]-0.01, frac_TPM_hS_genes_per_stage[i]-0.002, time_unique[i], fontsize=8, color='white', fontweight='bold')
plt.ylabel('% <TPM>', fontsize=15)
plt.xlabel('% different genes', fontsize=15)
plt.title('hS genes', fontsize=15, fontweight='bold')
plt.yticks(fontsize=13)
plt.xticks(fontsize=13)
plt.savefig(path_save_data+'frac_hS_genes_vs_frac_tpm_hS.png', dpi=600, bbox_inches='tight')
plt.show()


#====================================



#7.1.) Mean expression per gene in each class
i, ind_U, g=np.intersect1d(genes, U_id, return_indices=True)
i, ind_S, g=np.intersect1d(genes, S_id, return_indices=True)
i, ind_hS, g=np.intersect1d(genes, hS_id, return_indices=True)
bulk_U=mean_bulk_matrix[ind_U]
bulk_hS=mean_bulk_matrix[ind_hS]
bulk_S=mean_bulk_matrix[ind_S]

mean_exp_U=np.zeros(len(U_id))
for i in range(len(U_id)):
    ind_non_null=np.where(bulk_U[i, :]>0)[0]
    exp=bulk_U[i, ind_non_null]
    mean_exp_U[i]=np.mean(exp)
    
mean_exp_S=np.zeros(len(S_id))
for i in range(len(S_id)):
    ind_non_null=np.where(bulk_S[i, :]>0)[0]
    exp=bulk_S[i, ind_non_null]
    mean_exp_S[i]=np.mean(exp)
        
mean_exp_hS=np.zeros(len(hS_id))
for i in range(len(hS_id)):
    ind_non_null=np.where(bulk_hS[i, :]>0)[0]
    exp=bulk_hS[i, ind_non_null]
    mean_exp_hS[i]=np.mean(exp)

#7.1.1.) #figure hist
plt.figure(figsize=(4, 3), dpi=600)
plt.hist(mean_exp_U, color='royalblue', label='U', log=True, bins=100)
plt.hist(mean_exp_S, color='mediumturquoise', label='S', log=True, bins=100)
plt.hist(mean_exp_hS, color='tomato', label='hS', log=True, bins=100)
plt.xlabel('<TPM>', fontsize=18)
plt.ylabel('# genes', fontsize=18)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.legend(fontsize=13)
plt.grid(False)
plt.tight_layout()
# plt.xlim(0, 100)
plt.savefig(path_save_data+'non_zero_mean_gene_expression_bulk_histogram.png', dpi=600, bbox_inches='tight')
plt.show()



#7.1.2.) figure boxplot
data=[mean_exp_U, mean_exp_S, mean_exp_hS]
fig, ax = plt.subplots(figsize=(4,3), dpi=600)
box = ax.boxplot(data, labels=["U", "S", "hS"],  widths=0.45, patch_artist=True, 
                 flierprops=dict(marker='o', color='gray', markersize=0.5), 
                 boxprops=dict(edgecolor="none"))
colors = ["royalblue", "mediumturquoise", "tomato"]
for patch, color in zip(box["boxes"], colors):
    patch.set_facecolor(color)
plt.setp(box["medians"], color="black", linewidth=1)
plt.setp(box["whiskers"], color="black", linestyle="--", linewidth=0.5)
plt.setp(box["caps"], color="black", linewidth=0.5)
plt.ylabel("<TPM>", fontsize=15)
plt.xlabel("Gene classes", fontsize=15)
plt.xticks(fontsize=18, fontweight='bold')
plt.yticks(fontsize=13)
plt.yscale('log')
plt.savefig(path_save_data+'non_zero_mean_exp_bulk.png', dpi=600, bbox_inches='tight')
plt.show()


print('KS test between <TPM> distributions')
ks_U_S, pvalue_U_S = stats.ks_2samp(mean_exp_U, mean_exp_S)
print('U vs S:', ks_U_S, pvalue_U_S)
ks_U_hS, pvalue_U_hS = stats.ks_2samp(mean_exp_U, mean_exp_hS)
print('U vs hS:', ks_U_hS, pvalue_U_hS)
ks_S_hS, pvalue_S_hS = stats.ks_2samp(mean_exp_S, mean_exp_hS)
print('S vs hS:', ks_S_hS, pvalue_S_hS)
    




#8.) U, hS and S matrix
common_U, ind_all_U, ind_U=np.intersect1d(genes, U_id, return_indices=True)
U_bulk_matrix=mean_bulk_matrix[ind_all_U, :]
del common_U, ind_U
common_S, ind_all_S, ind_S=np.intersect1d(genes, S_id, return_indices=True)
S_bulk_matrix=mean_bulk_matrix[ind_all_S, :]
del common_S, ind_S
common_hS, ind_all_hS, ind_hS=np.intersect1d(genes, hS_id, return_indices=True)
hS_bulk_matrix=mean_bulk_matrix[ind_all_hS, :]
del common_hS, ind_hS


# 0 - 0.25 → Low inequality
# 0.25 - 0.40 → Moderate inequality
# 0.40 - 0.55 → High inequality
# 0.55 - 1.00 → Extreme inequality


#10.) Analysis of the std and mean of <TPM> per stage per gene class to undertand the cv

#10.1.) std of <TPM> per gene class
#10.1.1.) std hS genes
std_bulk_matrix_hS=np.zeros((len(hS_id), len(time_unique)))
common_hS, ind_all_hS, ind_hS=np.intersect1d(genes, hS_id, return_indices=True)

hs_all=bulk_matrix[ind_all_hS]

for k in range(len(hS_id)):
    count=0
    for i in range(len(time_unique)):
        sum_embryo=0
        for j in range(len(embryo)):
            sum_embryo=sum_embryo+(hs_all[k][int(count)]-hS_bulk_matrix[k][i])**2
            count=count+1
        std_bulk_matrix_hS[k][i]=np.sqrt(sum_embryo/(len(embryo)-1))

plt.hist(std_bulk_matrix_hS.flatten(), bins=100, log=True)
plt.show()

#10.1.2.) std S genes
std_bulk_matrix_S=np.zeros((len(S_id), len(time_unique)))
common_S, ind_all_S, ind_S=np.intersect1d(genes, S_id, return_indices=True)

S_all=bulk_matrix[ind_all_S]

for k in range(len(S_id)):
    count=0
    for i in range(len(time_unique)):
        sum_embryo=0
        for j in range(len(embryo)):
            sum_embryo=sum_embryo+(S_all[k][int(count)]-S_bulk_matrix[k][i])**2
            count=count+1
        std_bulk_matrix_S[k][i]=np.sqrt(sum_embryo/(len(embryo)-1))

plt.hist(std_bulk_matrix_S.flatten(), bins=100, log=True)
plt.show()

#10.1.3.) std U genes
std_bulk_matrix_U=np.zeros((len(U_id), len(time_unique)))
common_U, ind_all_U, ind_U=np.intersect1d(genes, U_id, return_indices=True)

U_all=bulk_matrix[ind_all_U]

for k in range(len(U_id)):
    count=0
    for i in range(len(time_unique)):
        sum_embryo=0
        for j in range(len(embryo)):
            sum_embryo=sum_embryo+(U_all[k][int(count)]-U_bulk_matrix[k][i])**2
            count=count+1
        std_bulk_matrix_U[k][i]=np.sqrt(sum_embryo/(len(embryo)-1))

plt.hist(std_bulk_matrix_U.flatten(), bins=100, log=True)
plt.show()

#10.1.4.) We plot the median of the std <TPM> per stage in each gene class
median_U_std=np.nanmedian(std_bulk_matrix_U, axis=0)
median_S_std=np.nanmedian(std_bulk_matrix_S, axis=0)
median_hS_std=np.nanmedian(std_bulk_matrix_hS, axis=0)

plt.figure(figsize=(11, 6), dpi=600)
plt.plot(np.linspace(0, len(time_unique)-1, len(time_unique)), median_U_std, marker='o', color='royalblue', label='U')
plt.plot(np.linspace(0, len(time_unique)-1, len(time_unique)), median_S_std, marker='o', color='mediumturquoise', label='S')
plt.plot(np.linspace(0, len(time_unique)-1, len(time_unique)), median_hS_std, marker='o', color='tomato', label='hS')
plt.xticks(ticks=np.arange(len(time_unique)), labels=time_unique, rotation=90, fontsize=25)
plt.yticks(fontsize=25)
plt.legend(fontsize=20)
plt.xlabel('Stages', fontsize=30, fontweight='bold')
plt.ylabel('std (median)', fontsize=26, fontweight='bold')
plt.savefig(path_save_data+'std_per_stage_comparing_embryo.png', dpi=600, bbox_inches='tight')
plt.show()


#10.2.) We plot the median of the mean expression <TPM> per stage in each gene class
median_U_mean=np.nanmedian(U_bulk_matrix, axis=0)
median_S_mean=np.nanmedian(S_bulk_matrix, axis=0)
median_hS_mean=np.nanmedian(hS_bulk_matrix, axis=0)


plt.figure(figsize=(11, 6), dpi=600)
plt.plot(np.linspace(0, len(time_unique)-1, len(time_unique)), median_U_mean, marker='o', color='royalblue', label='U')
plt.plot(np.linspace(0, len(time_unique)-1, len(time_unique)), median_S_mean, marker='o', color='mediumturquoise', label='S')
plt.plot(np.linspace(0, len(time_unique)-1, len(time_unique)), median_hS_mean, marker='o', color='tomato', label='hS')
plt.xticks(ticks=np.arange(len(time_unique)), labels=time_unique, rotation=90, fontsize=25)
plt.yticks(fontsize=25)
plt.legend(fontsize=20)
plt.xlabel('Stages', fontsize=30, fontweight='bold')
plt.ylabel('mean (median)', fontsize=26, fontweight='bold')
plt.savefig(path_save_data+'mean_per_stage_comparing_embryo.png', dpi=600, bbox_inches='tight')
plt.show()


#11.) Varialibility per gene between embryos in each stage
#ADJUSTED CV

std_bulk_matrix=np.zeros((len(genes), len(time_unique)))

for k in range(len(genes)):
    count=0
    for i in range(len(time_unique)):
        sum_embryo=0
        for j in range(len(embryo)):
            sum_embryo=sum_embryo+(bulk_matrix[k][int(count)]-mean_bulk_matrix[k][i])**2
            count=count+1
        std_bulk_matrix[k][i]=np.sqrt(sum_embryo/(len(embryo)-1))

# cv_genes=std_bulk_matrix/mean_bulk_matrix

cv_genes=np.zeros((len(genes), len(time_unique)))
for k in range(len(genes)):
    for i in range(len(time_unique)):    
        if mean_bulk_matrix[k][i]>0:
            cv_genes[k][i]=std_bulk_matrix[k][i]/mean_bulk_matrix[k][i]
        else:
            cv_genes[k][i]=np.log(0)

cv_genes = np.where(np.isinf(cv_genes), np.nan, cv_genes)



#cv log per gene and stage
log_cv2 = np.log10(cv_genes**2)

#normalize using rolling median
mean_genes = np.nanmean(mean_bulk_matrix, axis=1)

order_genes = np.argsort(mean_genes)

log_cv2_ordered = log_cv2[order_genes, :]

log_cv2_df = pd.DataFrame(log_cv2_ordered)
roll_median = log_cv2_df.rolling(window=50, min_periods=1, center=True).median()

#final norm
cv_norm = log_cv2_ordered - roll_median.values

# Reordenar a orden original de genes
reorder_index = np.argsort(order_genes)
cv_norm_original_order = cv_norm[reorder_index, :]

# cv_norm_abs = 10**cv_norm_original_order  # vuelve a escala CV²
# cv_norm_abs = np.sqrt(cv_norm_abs) 

common, ind_U, ind_i=np.intersect1d(genes, U_id, return_indices=True)
common, ind_S, ind_i=np.intersect1d(genes, S_id, return_indices=True)
common, ind_hS, ind_i=np.intersect1d(genes, hS_id, return_indices=True)

del common, ind_i

cv_genes_U=cv_norm_original_order[ind_U, :]
cv_genes_S=cv_norm_original_order[ind_S, :]
cv_genes_hS=cv_norm_original_order[ind_hS, :]


#11.1.) We compute the median in each stage (without taking into account null )
median_U=np.nanmedian(cv_genes_U, axis=0)
median_S=np.nanmedian(cv_genes_S, axis=0)
median_hS=np.nanmedian(cv_genes_hS, axis=0)

plt.figure(figsize=(11, 6), dpi=600)
plt.plot(np.linspace(0, len(time_unique)-1, len(time_unique)), median_U, marker='o', color='royalblue', label='U')
plt.plot(np.linspace(0, len(time_unique)-1, len(time_unique)), median_S, marker='o', color='mediumturquoise', label='S')
plt.plot(np.linspace(0, len(time_unique)-1, len(time_unique)), median_hS, marker='o', color='tomato', label='hS')
plt.xticks(ticks=np.arange(len(time_unique)), labels=time_unique, rotation=90, fontsize=25)
plt.yticks(fontsize=25)
plt.legend(fontsize=20)
plt.xlabel('Stages', fontsize=34)
plt.ylabel('CV adjusted (median)', fontsize=26, fontweight='bold')
plt.savefig(path_save_data+'cv_per_stage_comparing_embryo_norm.png', dpi=600, bbox_inches='tight')
plt.show()



#12.) Distrbution of protein coding genes in each class
gene_type_U=gene_type[ind_U]
n_prot_cod_U=np.where(gene_type_U=='protein_coding')[0]
gene_type_S=gene_type[ind_S]
n_prot_cod_S=np.where(gene_type_S=='protein_coding')[0]
gene_type_hS=gene_type[ind_hS]
n_prot_cod_hS=np.where(gene_type_hS=='protein_coding')[0]

frac_genes=np.array([len(n_prot_cod_U)/len(U_id), len(n_prot_cod_S)/len(S_id), len(n_prot_cod_hS)/len(hS_id)])
groups=['U', 'S', 'hS']
color=['royalblue', 'mediumturquoise', 'tomato']

plt.figure(figsize=(4, 3), dpi=600)
plt.bar(groups, frac_genes, color=color)
for i in range(len(groups)):
    plt.text(i-0.18, frac_genes[i]-0.1, np.round(frac_genes[i], 2), fontsize=12, fontweight='bold', color='white')
plt.xticks(ticks=np.arange(len(groups)), labels=groups, fontsize=16, fontweight='bold')
plt.yticks(fontsize=14)
plt.ylabel('Frac of protein\ncoding genes', fontsize=18, fontweight='bold')
plt.savefig(path_save_data+'frac_prot_coding_per_gene_class.png', dpi=600, bbox_inches='tight')
plt.show()

print(frac_genes)













