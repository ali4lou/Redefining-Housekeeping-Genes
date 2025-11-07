# -*- coding: utf-8 -*-
"""
Created on Mon Mar 17 15:15:52 2025

bulk analysis tissue - temperatures

@author: Alicia
"""


import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import warnings; warnings.simplefilter('ignore')
from scipy import stats
from scipy.stats import ranksums



path_save_data='PATH_WHERE_YOU_KEEP_YOUR_DATASETS'


tissue=['brain', 'gill', 'heart', 'intestine', 'kidney', 'liver', 'muscle', 'spleen']
temp = ['10', '18', '28']
# temp=['28']
label= [t + tp for t in tissue for tp in temp]


data=pd.read_csv(path_save_data+'%s.expression.txt' %label[0], sep='\t')
genes=np.array(list(data['Gene']))
genes, ind_genes=np.unique(genes, return_index=True)

#1.) We are going to analyse temperature 28ºC
tissue=['brain', 'gill', 'heart', 'intestine', 'kidney', 'liver', 'muscle', 'spleen']
temp=['28']
label= [t + tp for t in tissue for tp in temp]

matrix=np.zeros((len(genes), len(tissue)))
for i in range(len(label)):
    data=pd.read_csv(path_save_data+'%s.expression.txt' %label[i], sep='\t')
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


np.savetxt(path_save_data+'hS_id_tissue.txt', hS_tis, fmt='%s', delimiter=',')
np.savetxt(path_save_data+'S_id_tissue.txt', S_tis, fmt='%s', delimiter=',')
np.savetxt(path_save_data+'U_id_tissue.txt', U_tis, fmt='%s', delimiter=',')
np.savetxt(path_save_data+'null_genes_tissue.txt', null_genes, fmt='%s', delimiter=',')



np.median(act_tis_S)
np.median(act_tis_hS)
from matplotlib.ticker import MaxNLocator
data=[act_tis_U, act_tis_S, act_tis_hS]
fig, ax = plt.subplots(figsize=(4,3), dpi=600)
box = ax.boxplot(data, labels=["Ut", "St", "hSt"],  widths=0.45, patch_artist=True, 
                 flierprops=dict(marker='o', color='gray', markersize=0.5), 
                 boxprops=dict(edgecolor="none"))
colors = ["royalblue", "mediumturquoise", "tomato"]
for patch, color in zip(box["boxes"], colors):
    patch.set_facecolor(color)
plt.setp(box["medians"], color="black", linewidth=1)
plt.setp(box["whiskers"], color="black", linestyle="--", linewidth=0.5)
plt.setp(box["caps"], color="black", linewidth=0.5)
plt.ylabel("# expressed tissues", fontsize=16, fontweight='bold')
plt.xticks(fontsize=20, fontweight='bold')
plt.yticks(fontsize=13)
ax = plt.gca()
ax.yaxis.set_major_locator(MaxNLocator(integer=True))
plt.savefig(path_save_data+'act_tissue_per_class_tis.png', dpi=600, bbox_inches='tight')
plt.show()



#5.) We check if the classification match with ours

#5.1.) We read the bulk data
f=open(path_save_data+'genes_bulk.txt', 'r')
txt = f.read()
genes_bulk = txt.split('\n')
del txt, f
genes_bulk=np.delete(genes_bulk, len(genes_bulk)-1)
genes_bulk=np.array(genes_bulk)

f=open(path_save_data+'hS_id_bulk.txt', 'r')
txt = f.read()
hS_id_bulk = txt.split('\n')
del txt, f
hS_id_bulk=np.delete(hS_id_bulk, len(hS_id_bulk)-1)
hS_id_bulk=np.array(hS_id_bulk)

f=open(path_save_data+'S_id_bulk.txt', 'r')
txt = f.read()
S_id_bulk = txt.split('\n')
del txt, f
S_id_bulk=np.delete(S_id_bulk, len(S_id_bulk)-1)
S_id_bulk=np.array(S_id_bulk)

f=open(path_save_data+'U_id_bulk.txt', 'r')
txt = f.read()
U_id_bulk = txt.split('\n')
del txt, f
U_id_bulk=np.delete(U_id_bulk, len(U_id_bulk)-1)
U_id_bulk=np.array(U_id_bulk)

f=open(path_save_data+'null_genes.txt', 'r')
txt = f.read()
null_genes_id_bulk = txt.split('\n')
del txt, f
null_genes_id_bulk=np.delete(null_genes_id_bulk, len(null_genes_id_bulk)-1)
null_genes_id_bulk=np.array(null_genes_id_bulk)

f=open(path_save_data+'time_bulk.txt', 'r')
txt = f.read()
time_bulk = txt.split('\n')
del txt, f
time_bulk=np.delete(time_bulk, len(time_bulk)-1)
time_bulk=np.array(time_bulk)


mean_bulk_matrix= np.genfromtxt(path_save_data+'mean_bulk_matrix.txt', delimiter=',', dtype=None, encoding=None)


#5.2.) Common genes
common_genes=np.intersect1d(genes, genes_bulk)

common_genes_hS=np.intersect1d(genes, hS_id_bulk)
print(len(hS_id_bulk), len(common_genes_hS))
hS_ini_U_tis, ind_all_hS, ind_bulk_hS=np.intersect1d(U_tis, hS_id_bulk, return_indices=True)
hS_ini_S_tis, ind_all_hS, ind_bulk_hS=np.intersect1d(S_tis, hS_id_bulk, return_indices=True)
hS_ini_hS_tis, ind_all_hS, ind_bulk_hS=np.intersect1d(hS_tis, hS_id_bulk, return_indices=True)
hS_ini_null_tis, ind_all_hS, ind_bulk_hS=np.intersect1d(null_genes, hS_id_bulk, return_indices=True)

print(len(hS_ini_U_tis)+len(hS_ini_S_tis)+len(hS_ini_hS_tis)+len(hS_ini_null_tis))

np.savetxt(path_save_data+"hS_dev_U_tis.txt", hS_ini_U_tis, fmt="%s")
np.savetxt(path_save_data+"hS_dev_S_tis.txt", hS_ini_S_tis, fmt="%s")
np.savetxt(path_save_data+"hS_dev_hS_tis.txt", hS_ini_hS_tis, fmt="%s")
np.savetxt(path_save_data+"hS_dev_null_tis.txt", hS_ini_null_tis, fmt="%s")

frac_genes_hS=np.array([len(hS_ini_U_tis)/len(common_genes_hS), 
                     len(hS_ini_S_tis)/len(common_genes_hS), 
                     len(hS_ini_hS_tis)/len(common_genes_hS), 
                     len(hS_ini_null_tis)/len(common_genes_hS)])
groups=['U', 'S', 'hS', 'NE']
color=['gray', 'royalblue', 'mediumturquoise', 'tomato']

plt.figure(figsize=(4, 3), dpi=600)
plt.bar(groups, frac_genes_hS, color=color)
plt.xticks(ticks=np.arange(len(groups)), labels=groups, fontsize=15, fontweight='bold')
plt.yticks(fontsize=14)
plt.title('hS genes development', fontsize=16, fontweight='bold')
plt.ylabel('Fraction of genes', fontsize=16, fontweight='bold')
plt.xlabel('Tissue gene classes', fontsize=16, fontweight='bold')
for i in range(len(groups)):
    plt.text(i-0.25, frac_genes_hS[i], np.round(frac_genes_hS[i], 2), fontsize=12, fontweight='bold', color='black')
plt.savefig(path_save_data+'hS_genes_class.png', dpi=600, bbox_inches='tight')
plt.show()


common_genes_S=np.intersect1d(genes, S_id_bulk)
print(len(S_id_bulk), len(common_genes_S))
S_ini_U_tis, ind_all_hS, ind_bulk_hS=np.intersect1d(U_tis, S_id_bulk, return_indices=True)
S_ini_S_tis, ind_all_hS, ind_bulk_hS=np.intersect1d(S_tis, S_id_bulk, return_indices=True)
S_ini_hS_tis, ind_all_hS, ind_bulk_hS=np.intersect1d(hS_tis, S_id_bulk, return_indices=True)
S_ini_null_tis, ind_all_hS, ind_bulk_hS=np.intersect1d(null_genes, S_id_bulk, return_indices=True)

print(len(S_ini_U_tis)+len(S_ini_S_tis)+len(S_ini_hS_tis)+len(S_ini_null_tis))

np.savetxt(path_save_data+"S_dev_U_tis.txt", S_ini_U_tis, fmt="%s")
np.savetxt(path_save_data+"S_dev_S_tis.txt", S_ini_S_tis, fmt="%s")
np.savetxt(path_save_data+"S_dev_hS_tis.txt", S_ini_hS_tis, fmt="%s")
np.savetxt(path_save_data+"S_dev_null_tis.txt", S_ini_null_tis, fmt="%s")



frac_genes_S=np.array([len(S_ini_U_tis)/len(common_genes_S), 
                     len(S_ini_S_tis)/len(common_genes_S), 
                     len(S_ini_hS_tis)/len(common_genes_S), 
                     len(S_ini_null_tis)/len(common_genes_S)])
groups=['U', 'S', 'hS', 'NE']
color=['gray', 'royalblue', 'mediumturquoise', 'tomato']

plt.figure(figsize=(4, 3), dpi=600)
plt.bar(groups, frac_genes_S, color=color)
plt.xticks(ticks=np.arange(len(groups)), labels=groups, fontsize=15, fontweight='bold')
plt.yticks(fontsize=14)
plt.title('S genes development', fontsize=16, fontweight='bold')
plt.ylabel('Fraction of genes', fontsize=16, fontweight='bold')
plt.xlabel('Tissue gene classes', fontsize=16, fontweight='bold')
for i in range(len(groups)):
    plt.text(i-0.25, frac_genes_S[i], np.round(frac_genes_S[i], 2), fontsize=12, fontweight='bold', color='black')
plt.savefig(path_save_data+'S_genes_class.png', dpi=600, bbox_inches='tight')
plt.show()


common_genes_U=np.intersect1d(genes, U_id_bulk)
print(len(U_id_bulk), len(common_genes_U))
U_ini_U_tis, ind_all_hS, ind_bulk_hS=np.intersect1d(U_tis, U_id_bulk, return_indices=True)
U_ini_S_tis, ind_all_hS, ind_bulk_hS=np.intersect1d(S_tis, U_id_bulk, return_indices=True)
U_ini_hS_tis, ind_all_hS, ind_bulk_hS=np.intersect1d(hS_tis, U_id_bulk, return_indices=True)
U_ini_null_tis, ind_all_hS, ind_bulk_hS=np.intersect1d(null_genes, U_id_bulk, return_indices=True)

print(len(U_ini_U_tis)+len(U_ini_S_tis)+len(U_ini_hS_tis)+len(U_ini_null_tis))

np.savetxt(path_save_data+"U_dev_U_tis.txt", U_ini_U_tis, fmt="%s")
np.savetxt(path_save_data+"U_dev_S_tis.txt", U_ini_S_tis, fmt="%s")
np.savetxt(path_save_data+"U_dev_hS_tis.txt", U_ini_hS_tis, fmt="%s")
np.savetxt(path_save_data+"U_dev_null_tis.txt", U_ini_null_tis, fmt="%s")



frac_genes_U=np.array([len(U_ini_U_tis)/len(common_genes_U), 
                     len(U_ini_S_tis)/len(common_genes_U), 
                     len(U_ini_hS_tis)/len(common_genes_U), 
                     len(U_ini_null_tis)/len(common_genes_U)])
groups=['U', 'S', 'hS', 'NE']
color=['gray', 'royalblue', 'mediumturquoise', 'tomato']

plt.figure(figsize=(4, 3), dpi=600)
plt.bar(groups, frac_genes_U, color=color)
plt.xticks(ticks=np.arange(len(groups)), labels=groups, fontsize=15, fontweight='bold')
plt.yticks(fontsize=14)
plt.title('U genes development', fontsize=16, fontweight='bold')
plt.ylabel('Fraction of genes', fontsize=16, fontweight='bold')
plt.xlabel('Tissue gene classes', fontsize=16, fontweight='bold')
for i in range(len(groups)):
    plt.text(i-0.34, frac_genes_U[i], np.round(frac_genes_U[i], 3), fontsize=12, fontweight='bold', color='black')
plt.savefig(path_save_data+'U_genes_class.png', dpi=600, bbox_inches='tight')
plt.show()



common_genes_null=np.intersect1d(genes, null_genes_id_bulk)

len(U_tis)+len(S_tis)+len(hS_tis)+len(null_genes)

null_ini_U_tis, ind_all_hS, ind_bulk_hS=np.intersect1d(U_tis, null_genes_id_bulk, return_indices=True)
null_ini_S_tis, ind_all_hS, ind_bulk_hS=np.intersect1d(S_tis, null_genes_id_bulk, return_indices=True)
null_ini_hS_tis, ind_all_hS, ind_bulk_hS=np.intersect1d(hS_tis, null_genes_id_bulk, return_indices=True)
null_ini_null_tis, ind_all_hS, ind_bulk_hS=np.intersect1d(null_genes, null_genes_id_bulk, return_indices=True)

np.savetxt(path_save_data+"null_ini_U_tis.txt", null_ini_U_tis, fmt="%s")
np.savetxt(path_save_data+"null_ini_S_tis.txt", null_ini_S_tis, fmt="%s")
np.savetxt(path_save_data+"null_ini_hS_tis.txt", null_ini_hS_tis, fmt="%s")
np.savetxt(path_save_data+"null_ini_null_tis.txt", null_ini_null_tis, fmt="%s")

frac_genes_null=np.array([len(null_ini_U_tis)/len(common_genes_null), 
                     len(null_ini_S_tis)/len(common_genes_null), 
                     len(null_ini_hS_tis)/len(common_genes_null), 
                     len(null_ini_null_tis)/len(common_genes_null)])


#Matrix dotplot 
gene_names = ["Ut", "St", "hSt", 'NE']
categories = ['NE', "hS", "S", "U"]
frac_genes = np.array([frac_genes_null, frac_genes_hS, frac_genes_S, frac_genes_U])
n_genes_plot=np.array([frac_genes_null*len(common_genes_null), frac_genes_hS*len(common_genes_hS), frac_genes_S*len(common_genes_S), frac_genes_U*len(common_genes_U)])

fig, ax = plt.subplots(figsize=(4, 3))
color=['gray', 'gray', 'gray', 'gray',
       'tomato', 'tomato', 'tomato', 'tomato', 
       'mediumturquoise', 'mediumturquoise', 'mediumturquoise', 'mediumturquoise', 
       'royalblue', 'royalblue', 'royalblue', 'royalblue']
count=0
for i, category in enumerate(categories):
    for j, gene in enumerate(gene_names):
        plt.scatter(j+0.5, i+0.5, s=frac_genes[i][j] * 2000, color=color[count])
        plt.text(j + 0.5, i + 0.5, int(np.round(n_genes_plot[i][j], 0)), color='black',
                 fontsize=10, ha='center', va='center', fontweight='bold')
        print(n_genes_plot[i][j])
        count=count+1
# Add dashed guide lines at midpoints
for i in range(len(categories)-1):
    plt.axhline(y=i+1, color='gray', linestyle='--', alpha=0.5, lw=1)
for j in range(len(gene_names)-1):
    plt.axvline(x=j+1, color='gray', linestyle='--', alpha=0.5, lw=1)
ax.set_xticks(np.arange(len(gene_names)) + 0.5)
ax.set_xticklabels(gene_names, fontsize=18, fontweight='bold')
ax.set_yticks(np.arange(len(categories)) + 0.5)
ax.set_yticklabels(categories, fontsize=18, fontweight='bold')
plt.ylabel("Developmental\nclassification", fontsize=20, fontweight='bold', color='gray')
plt.xlim(0, 4)
plt.ylim(0, 4)

# Add row totals on the right
row_totals = [len(common_genes_null), len(common_genes_hS), len(common_genes_S), len(common_genes_U)]
for i, total in enumerate(row_totals):
    plt.text(len(gene_names) + 0.5, i + 0.5, str(int(total)),
             va='center', ha='center', fontsize=15)

col_totals = np.sum(n_genes_plot, axis=0)  # Asegúrate de que el orden coincida con gene_names
for j, total in enumerate(col_totals):
    plt.text(j + 0.5, -0.2, str(int(total)),  # y = 0 lo coloca debajo
             va='center', ha='center', fontsize=15)


plt.xlabel("Adult tissue", fontsize=20, fontweight='bold', color='gray')
ax.xaxis.set_label_position('top') 
ax.xaxis.tick_top()
plt.savefig(path_save_data+'dev_class_to_tissue.png', dpi=600, bbox_inches='tight')
plt.show()




#6.) We search for differences along development
def development(gene_subset, title):
    g, ind_dev, j=np.intersect1d(genes_bulk, gene_subset, return_indices=True)
    del g, j
    matrix_subset=mean_bulk_matrix[ind_dev, :]
    act_dev_stages=np.zeros(len(gene_subset))
    mean_exp=np.zeros(len(gene_subset))
    frac_genes_per_stage=np.zeros(len(time_bulk))
    mean_tpm_stage=np.zeros(len(time_bulk))
    for k in range(len(gene_subset)):
        act_dev_stages[k]=len(np.where(matrix_subset[k, :]>1)[0])
        ind_non_zero=np.where(matrix_subset[k, :]>0)[0]
        mean_exp[k]=np.mean(matrix_subset[k, ind_non_zero])
    for k in range(len(time_bulk)):
        ind_genes_non_zero=np.where(matrix_subset[:, k]>0)[0]
        frac_genes_per_stage[k]=len(np.where(matrix_subset[:, k]>1)[0])/len(gene_subset)
        mean_tpm_stage[k]=np.mean(matrix_subset[ind_genes_non_zero, k])
    
    high_stage=np.zeros(len(time_bulk))
    for i in range(len(gene_subset)):
        max_stage=np.argmax(matrix_subset[i, :])
        high_stage[int(max_stage)]=high_stage[int(max_stage)]+1
        
    
    #Figure of m_types
    matrix_subset_represent=np.where(matrix_subset < 1, np.log(0), matrix_subset)
    plt.figure(figsize=(5, 4), dpi=600)
    plt.imshow(matrix_subset_represent , cmap='viridis', aspect='auto')
    cbar=plt.colorbar(shrink=0.5, aspect=15)
    cbar.set_label('TPM', size=20)  # Aquí puedes ajustar el tamaño como desees
    plt.grid(False)
    plt.ylabel('Genes', fontsize=20)
    plt.xlabel('Stage', fontsize=20)
    plt.xticks(np.arange(len(time_bulk)), time_bulk, fontsize=15, rotation=90)  # Set text labels.
    plt.title(title)
    # plt.yticks(np.arange(len(cell_type_unique)), cell_type_unique, fontsize=15)
    # plt.savefig(path_save_data+'percentile_hS_genes_pseudo_bulk.png', dpi=600, bbox_inches='tight')  
    plt.show()   
    
    return act_dev_stages, mean_exp, frac_genes_per_stage, mean_tpm_stage, high_stage/len(gene_subset)



#6.1.) hS genes
act_dev_stages_hS_ini_U_tis, mean_exp_hS_ini_U_tis, frac_stage_hS_ini_U_tis, mean_tpm_stage_hS_ini_U_tis, max_stage_hS_ini_U_tis=development(hS_ini_U_tis, 'hS dev - U tis')
act_dev_stages_hS_ini_S_tis, mean_exp_hS_ini_S_tis, frac_stage_hS_ini_S_tis, mean_tpm_stage_hS_ini_S_tis, max_stage_hS_ini_S_tis=development(hS_ini_S_tis, 'hS dev - S tis')
act_dev_stages_hS_ini_hS_tis, mean_exp_hS_ini_hS_tis, frac_stage_hS_ini_hS_tis, mean_tpm_stage_hS_ini_hS_tis, max_stage_hS_ini_hS_tis=development(hS_ini_hS_tis, 'hS dev - hS tis')
act_dev_stages_hS_ini_null_tis, mean_exp_hS_ini_null_tis, frac_stage_hS_ini_null_tis, mean_tpm_stage_hS_ini_null_tis, max_stage_hS_ini_null_tis=development(hS_ini_null_tis, 'hS dev - zero tis')

#maximal stage
time_unique=np.array(['0hpf', '0.75hpf', '2.25hpf', '3hpf', '4.3hpf', '5.25hpf', '6hpf', '8hpf', '10.3hpf', '16hpf', '19hpf', '24hpf', '30hpf', '36hpf', '2dpf', '3dpf', '4dpf', '5dpf'])

width = 0.25  # ancho de cada barra
x = np.arange(len(time_unique))
plt.figure(figsize=(11, 6), dpi=600)
plt.bar(x - width, max_stage_hS_ini_U_tis*100, width=width, label='hS-Ut', color="mediumvioletred")
plt.bar(x,         max_stage_hS_ini_S_tis*100, width=width, label='hS-St', color="hotpink")
plt.bar(x + width, max_stage_hS_ini_hS_tis*100, width=width, label='hS-hSt', color="pink")

plt.xticks(ticks=x, labels=time_unique, rotation=90, fontsize=25)
plt.yticks(fontsize=25)
plt.xlabel('Stages', fontsize=30)
plt.ylabel('% genes\nwith max exp', fontsize=28, fontweight='bold')
plt.legend(fontsize=20)

plt.tight_layout()
plt.savefig(path_save_data + 'max_stage_gene_subclass_hS.png', dpi=600, bbox_inches='tight')
plt.show()


plt.figure(figsize=(11, 6), dpi=600)
plt.plot(time_bulk, frac_stage_hS_ini_U_tis, marker='o', label='U tissue', color='royalblue')
plt.plot(time_bulk, frac_stage_hS_ini_S_tis, marker='o', label='S tissue', color='mediumturquoise')
plt.plot(time_bulk, frac_stage_hS_ini_hS_tis, marker='o', label='hS tissue', color='tomato')
plt.plot(time_bulk, frac_stage_hS_ini_null_tis, marker='o', label='null tissue', color='gray')
plt.xticks(ticks=np.arange(0, len(time_bulk)), labels=time_bulk, rotation=90, fontsize=25)
plt.title('hS genes development', fontsize=30, fontweight='bold')
plt.ylabel('% genes', fontsize=25, fontweight='bold')
plt.xlabel('Stages', fontsize=30)
plt.yticks(fontsize=20)
plt.legend(
    fontsize=25,  
    loc='center left', 
    bbox_to_anchor=(1, 0.7)  
)
plt.tight_layout(rect=[0, 0, 1, 1]) 
plt.savefig(path_save_data+'frac_gene_per_stage_hS.png', dpi=600, bbox_inches='tight')
plt.show()

plt.figure(figsize=(11, 6), dpi=600)
plt.plot(time_bulk, mean_tpm_stage_hS_ini_U_tis, label='U tissue',  marker='o', color='royalblue')
plt.plot(time_bulk, mean_tpm_stage_hS_ini_S_tis, label='S tissue', marker='o', color='mediumturquoise')
plt.plot(time_bulk, mean_tpm_stage_hS_ini_hS_tis, label='hS tissue', marker='o', color='tomato')
plt.plot(time_bulk, mean_tpm_stage_hS_ini_null_tis, label='null tissue', marker='o', color='gray')
plt.xticks(ticks=np.arange(0, len(time_bulk)), labels=time_bulk, rotation=90, fontsize=25)
plt.ylabel('<TPM> (non zero)', fontsize=25, fontweight='bold')
plt.xlabel('Stages', fontsize=30)
plt.title('hS genes development', fontsize=30, fontweight='bold')
plt.yticks(fontsize=20)
plt.legend(
    fontsize=25,  
    loc='center left', 
    bbox_to_anchor=(1, 0.7)  
)
plt.tight_layout(rect=[0, 0, 1, 1]) 
plt.savefig(path_save_data+'mean_gene_exp_per_stage_hS.png', dpi=600, bbox_inches='tight')
plt.show()

data=[mean_exp_hS_ini_U_tis, mean_exp_hS_ini_S_tis, mean_exp_hS_ini_hS_tis, mean_exp_hS_ini_null_tis]
fig, ax = plt.subplots(figsize=(4,3), dpi=600)
box = ax.boxplot(data, labels=["U", "S", "hS", 'null'],  widths=0.45, patch_artist=True, 
                 flierprops=dict(marker='o', color='gray', markersize=0.5), 
                 boxprops=dict(edgecolor="none"))
colors = ["royalblue", "mediumturquoise", "tomato", 'gray']
for patch, color in zip(box["boxes"], colors):
    patch.set_facecolor(color)
plt.setp(box["medians"], color="black", linewidth=1)
plt.setp(box["whiskers"], color="black", linestyle="--", linewidth=0.5)
plt.setp(box["caps"], color="black", linewidth=0.5)
plt.ylabel("<TPM> (non zero)", fontsize=15)
plt.xlabel("Tissue gene classes", fontsize=15)
plt.xticks(fontsize=18, fontweight='bold')
plt.title('hS genes development', fontsize=15, fontweight='bold')
plt.yticks(fontsize=13)
plt.yscale('log')
plt.savefig(path_save_data+'non_zero_mean_exp_bulk_hS.png', dpi=600, bbox_inches='tight')
plt.show()


data=[act_dev_stages_hS_ini_U_tis, act_dev_stages_hS_ini_S_tis, act_dev_stages_hS_ini_hS_tis, act_dev_stages_hS_ini_null_tis]
fig, ax = plt.subplots(figsize=(4,3), dpi=600)
box = ax.boxplot(data, labels=["U", "S", "hS", 'null'],  widths=0.45, patch_artist=True, 
                 flierprops=dict(marker='o', color='gray', markersize=0.5), 
                 boxprops=dict(edgecolor="none"))
colors = ["royalblue", "mediumturquoise", "tomato", 'gray']
for patch, color in zip(box["boxes"], colors):
    patch.set_facecolor(color)
plt.setp(box["medians"], color="black", linewidth=1)
plt.setp(box["whiskers"], color="black", linestyle="--", linewidth=0.5)
plt.setp(box["caps"], color="black", linewidth=0.5)
plt.ylabel('% genes with max exp', fontsize=15)
plt.xlabel("Tissue gene classes", fontsize=15)
plt.xticks(fontsize=18, fontweight='bold')
plt.yticks(fontsize=13)
plt.title('hS genes development', fontsize=15, fontweight='bold')
# plt.yscale('log')
plt.savefig(path_save_data+'act_dev_stages_bulk_hS.png.png', dpi=600, bbox_inches='tight')
plt.show()

from scipy import stats

print('KS test between <TPM> distributions')
ks_U_S, pvalue_U_S = stats.ks_2samp(mean_exp_hS_ini_U_tis, mean_exp_hS_ini_S_tis)
print('U vs S:', ks_U_S, pvalue_U_S)
ks_U_hS, pvalue_U_hS = stats.ks_2samp(mean_exp_hS_ini_U_tis, mean_exp_hS_ini_hS_tis)
print('U vs hS:', ks_U_hS, pvalue_U_hS)
ks_S_hS, pvalue_S_hS = stats.ks_2samp(mean_exp_hS_ini_S_tis, mean_exp_hS_ini_hS_tis)
print('S vs hS:', ks_S_hS, pvalue_S_hS)
print('median hS tissue:', np.median(mean_exp_hS_ini_hS_tis))
print('median U tissue:', np.median(mean_exp_hS_ini_U_tis))


#6.2.) S genes
act_dev_stages_S_ini_U_tis, mean_exp_S_ini_U_tis, frac_stage_S_ini_U_tis, mean_tpm_stage_S_ini_U_tis, max_stage_S_ini_U_tis=development(S_ini_U_tis, 'S dev - U tis')
act_dev_stages_S_ini_S_tis, mean_exp_S_ini_S_tis, frac_stage_S_ini_S_tis, mean_tpm_stage_S_ini_S_tis, max_stage_S_ini_S_tis=development(S_ini_S_tis, 'S dev - S tis')
act_dev_stages_S_ini_hS_tis, mean_exp_S_ini_hS_tis, frac_stage_S_ini_hS_tis, mean_tpm_stage_S_ini_hS_tis, max_stage_S_ini_hS_tis=development(S_ini_hS_tis, 'S dev - hS tis')
act_dev_stages_S_ini_null_tis, mean_exp_S_ini_null_tis, frac_stage_S_ini_null_tis, mean_tpm_stage_S_ini_null_tis, max_stage_S_ini_null_tis=development(S_ini_null_tis, 'S dev - zero tis')


#maximal stage
time_unique=np.array(['0hpf', '0.75hpf', '2.25hpf', '3hpf', '4.3hpf', '5.25hpf', '6hpf', '8hpf', '10.3hpf', '16hpf', '19hpf', '24hpf', '30hpf', '36hpf', '2dpf', '3dpf', '4dpf', '5dpf'])

width = 0.25  # ancho de cada barra
x = np.arange(len(time_unique))
plt.figure(figsize=(11, 6), dpi=600)
plt.bar(x - width, max_stage_S_ini_U_tis * 100, width=width, label='S-Ut', color="#7E57C2")
plt.bar(x,         max_stage_S_ini_S_tis * 100, width=width, label='S-St', color="#AB47BC")
plt.bar(x + width, max_stage_S_ini_hS_tis * 100, width=width, label='S-hSt', color="#EC407A")

plt.xticks(ticks=x, labels=time_unique, rotation=90, fontsize=25)
plt.yticks(fontsize=25)
plt.xlabel('Stages', fontsize=30)
plt.ylabel('% genes\nwith max exp', fontsize=28, fontweight='bold')
plt.legend(fontsize=20)

plt.tight_layout()
plt.savefig(path_save_data + 'max_stage_gene_subclass_S.png', dpi=600, bbox_inches='tight')
plt.show()

plt.figure(figsize=(11, 6), dpi=600)
plt.plot(time_bulk, frac_stage_S_ini_U_tis, marker='o', label='U tissue', color='royalblue')
plt.plot(time_bulk, frac_stage_S_ini_S_tis, marker='o', label='S tissue', color='mediumturquoise')
plt.plot(time_bulk, frac_stage_S_ini_hS_tis, marker='o', label='hS tissue', color='tomato')
plt.plot(time_bulk, frac_stage_S_ini_null_tis, marker='o', label='null tissue', color='gray')
plt.xticks(ticks=np.arange(0, len(time_bulk)), labels=time_bulk, rotation=90, fontsize=25)
plt.title('S genes development', fontsize=30, fontweight='bold')
plt.ylabel('% genes', fontsize=25, fontweight='bold')
plt.xlabel('Stages', fontsize=30)
plt.yticks(fontsize=20)
plt.legend(
    fontsize=25,  
    loc='center left', 
    bbox_to_anchor=(1, 0.7)  
)
plt.tight_layout(rect=[0, 0, 1, 1]) 
plt.savefig(path_save_data+'frac_gene_per_stage_S.png', dpi=600, bbox_inches='tight')
plt.show()

plt.figure(figsize=(11, 6), dpi=600)
plt.plot(time_bulk, mean_tpm_stage_S_ini_U_tis, label='U tissue',  marker='o', color='royalblue')
plt.plot(time_bulk, mean_tpm_stage_S_ini_S_tis, label='S tissue', marker='o', color='mediumturquoise')
plt.plot(time_bulk, mean_tpm_stage_S_ini_hS_tis, label='hS tissue', marker='o', color='tomato')
plt.plot(time_bulk, mean_tpm_stage_S_ini_null_tis, label='null tissue', marker='o', color='gray')
plt.xticks(ticks=np.arange(0, len(time_bulk)), labels=time_bulk, rotation=90, fontsize=25)
plt.ylabel('<TPM> (non zero)', fontsize=25, fontweight='bold')
plt.xlabel('Stages', fontsize=30)
plt.title('S genes development', fontsize=30, fontweight='bold')
plt.yticks(fontsize=20)
plt.legend(
    fontsize=25,  
    loc='center left', 
    bbox_to_anchor=(1, 0.7)  
)
plt.tight_layout(rect=[0, 0, 1, 1]) 
plt.savefig(path_save_data+'mean_gene_exp_per_stage_S.png', dpi=600, bbox_inches='tight')
plt.show()

data=[mean_exp_S_ini_U_tis, mean_exp_S_ini_S_tis, mean_exp_S_ini_hS_tis, mean_exp_S_ini_null_tis]
fig, ax = plt.subplots(figsize=(4,3), dpi=600)
box = ax.boxplot(data, labels=["U", "S", "hS", 'null'],  widths=0.45, patch_artist=True, 
                 flierprops=dict(marker='o', color='gray', markersize=0.5), 
                 boxprops=dict(edgecolor="none"))
colors = ["royalblue", "mediumturquoise", "tomato", 'gray']
for patch, color in zip(box["boxes"], colors):
    patch.set_facecolor(color)
plt.setp(box["medians"], color="black", linewidth=1)
plt.setp(box["whiskers"], color="black", linestyle="--", linewidth=0.5)
plt.setp(box["caps"], color="black", linewidth=0.5)
plt.ylabel("<TPM> (non zero)", fontsize=15)
plt.xlabel("Tissue gene classes", fontsize=15)
plt.xticks(fontsize=18, fontweight='bold')
plt.title('S genes development', fontsize=15, fontweight='bold')
plt.yticks(fontsize=13)
plt.yscale('log')
plt.savefig(path_save_data+'non_zero_mean_exp_bulk_S.png', dpi=600, bbox_inches='tight')
plt.show()


data=[act_dev_stages_S_ini_U_tis, act_dev_stages_S_ini_S_tis, act_dev_stages_S_ini_hS_tis, act_dev_stages_S_ini_null_tis]
fig, ax = plt.subplots(figsize=(4,3), dpi=600)
box = ax.boxplot(data, labels=["U", "S", "hS", 'null'],  widths=0.45, patch_artist=True, 
                 flierprops=dict(marker='o', color='gray', markersize=0.5), 
                 boxprops=dict(edgecolor="none"))
colors = ["royalblue", "mediumturquoise", "tomato", 'gray']
for patch, color in zip(box["boxes"], colors):
    patch.set_facecolor(color)
plt.setp(box["medians"], color="black", linewidth=1)
plt.setp(box["whiskers"], color="black", linestyle="--", linewidth=0.5)
plt.setp(box["caps"], color="black", linewidth=0.5)
plt.ylabel("# active stages\n(<TPM> > 1)", fontsize=15)
plt.xlabel("Tissue gene classes", fontsize=15)
plt.xticks(fontsize=18, fontweight='bold')
plt.yticks(fontsize=13)
plt.title('S genes development', fontsize=15, fontweight='bold')
# plt.yscale('log')
plt.savefig(path_save_data+'act_dev_stages_bulk_S.png', dpi=600, bbox_inches='tight')
plt.show()


print('KS test between <TPM> distributions')
ks_U_S, pvalue_U_S = stats.ks_2samp(mean_exp_S_ini_U_tis, mean_exp_S_ini_S_tis)
print('U vs S:', ks_U_S, pvalue_U_S)
ks_U_hS, pvalue_U_hS = stats.ks_2samp(mean_exp_S_ini_U_tis, mean_exp_S_ini_hS_tis)
print('U vs hS:', ks_U_hS, pvalue_U_hS)
ks_S_hS, pvalue_S_hS = stats.ks_2samp(mean_exp_S_ini_S_tis, mean_exp_S_ini_hS_tis)
print('S vs hS:', ks_S_hS, pvalue_S_hS)
print('median hS tissue:', np.median(mean_exp_S_ini_hS_tis))
print('median U tissue:', np.median(mean_exp_S_ini_U_tis))



#6.2.) U genes
act_dev_stages_U_ini_U_tis, mean_exp_U_ini_U_tis, frac_stage_U_ini_U_tis, mean_tpm_stage_U_ini_U_tis, max_stage_U_ini_U_tis=development(U_ini_U_tis, 'U dev - U tis')
act_dev_stages_U_ini_S_tis, mean_exp_U_ini_S_tis, frac_stage_U_ini_S_tis, mean_tpm_stage_U_ini_S_tis, max_stage_U_ini_S_tis=development(U_ini_S_tis, 'U dev - S tis')
act_dev_stages_U_ini_hS_tis, mean_exp_U_ini_hS_tis, frac_stage_U_ini_hS_tis, mean_tpm_stage_U_ini_hS_tis, max_stage_U_ini_hS_tis=development(U_ini_hS_tis, 'U dev - hS tis')
act_dev_stages_U_ini_null_tis, mean_exp_U_ini_null_tis, frac_stage_U_ini_null_tis, mean_tpm_stage_U_ini_null_tis, max_stage_U_ini_null_tis=development(U_ini_null_tis, 'U dev - zero tis')


#maximal stage
time_unique=np.array(['0hpf', '0.75hpf', '2.25hpf', '3hpf', '4.3hpf', '5.25hpf', '6hpf', '8hpf', '10.3hpf', '16hpf', '19hpf', '24hpf', '30hpf', '36hpf', '2dpf', '3dpf', '4dpf', '5dpf'])

width = 0.25  # ancho de cada barra
x = np.arange(len(time_unique))
plt.figure(figsize=(11, 6), dpi=600)
plt.bar(x - width, max_stage_U_ini_U_tis * 100, width=width, label='U-Ut', color="#0D47A1")
plt.bar(x,         max_stage_U_ini_S_tis * 100, width=width, label='U-St', color="#1976D2")
plt.bar(x + width, max_stage_U_ini_hS_tis * 100, width=width, label='U-hSt', color="#5C6BC0")

plt.xticks(ticks=x, labels=time_unique, rotation=90, fontsize=25)
plt.yticks(fontsize=25)
plt.xlabel('Stages', fontsize=30)
plt.ylabel('% genes\nwith max exp', fontsize=28, fontweight='bold')
plt.legend(fontsize=20)

plt.tight_layout()
plt.savefig(path_save_data + 'max_stage_gene_subclass_U.png', dpi=600, bbox_inches='tight')
plt.show()

plt.figure(figsize=(11, 6), dpi=600)
plt.plot(time_bulk, frac_stage_U_ini_U_tis, marker='o', label='U tissue', color='royalblue')
plt.plot(time_bulk, frac_stage_U_ini_S_tis, marker='o', label='S tissue', color='mediumturquoise')
plt.plot(time_bulk, frac_stage_U_ini_hS_tis, marker='o', label='hS tissue', color='tomato')
plt.plot(time_bulk, frac_stage_U_ini_null_tis, marker='o', label='null tissue', color='gray')
plt.xticks(ticks=np.arange(0, len(time_bulk)), labels=time_bulk, rotation=90, fontsize=25)
plt.title('U genes development', fontsize=30, fontweight='bold')
plt.ylabel('% genes', fontsize=25, fontweight='bold')
plt.xlabel('Stages', fontsize=30)
plt.yticks(fontsize=20)
plt.legend(
    fontsize=25,  
    loc='center left', 
    bbox_to_anchor=(1, 0.7)  
)
plt.tight_layout(rect=[0, 0, 1, 1]) 
plt.savefig(path_save_data+'frac_gene_per_stage_U.png', dpi=600, bbox_inches='tight')
plt.show()

plt.figure(figsize=(11, 6), dpi=600)
plt.plot(time_bulk, mean_tpm_stage_U_ini_U_tis, label='U tissue',  marker='o', color='royalblue')
plt.plot(time_bulk, mean_tpm_stage_U_ini_S_tis, label='S tissue', marker='o', color='mediumturquoise')
plt.plot(time_bulk, mean_tpm_stage_U_ini_hS_tis, label='hS tissue', marker='o', color='tomato')
plt.plot(time_bulk, mean_tpm_stage_U_ini_null_tis, label='null tissue', marker='o', color='gray')
plt.xticks(ticks=np.arange(0, len(time_bulk)), labels=time_bulk, rotation=90, fontsize=25)
plt.ylabel('<TPM> (non zero)', fontsize=25, fontweight='bold')
plt.xlabel('Stages', fontsize=30)
plt.title('U genes development', fontsize=30, fontweight='bold')
plt.yticks(fontsize=20)
plt.legend(
    fontsize=25,  
    loc='center left', 
    bbox_to_anchor=(1, 0.7)  
)
plt.tight_layout(rect=[0, 0, 1, 1]) 
plt.savefig(path_save_data+'mean_gene_exp_per_stage_U.png', dpi=600, bbox_inches='tight')
plt.show()

data=[mean_exp_U_ini_U_tis, mean_exp_U_ini_S_tis, mean_exp_U_ini_hS_tis, mean_exp_U_ini_null_tis]
fig, ax = plt.subplots(figsize=(4,3), dpi=600)
box = ax.boxplot(data, labels=["U", "S", "hS", 'null'],  widths=0.45, patch_artist=True, 
                 flierprops=dict(marker='o', color='gray', markersize=0.5), 
                 boxprops=dict(edgecolor="none"))
colors = ["royalblue", "mediumturquoise", "tomato", 'gray']
for patch, color in zip(box["boxes"], colors):
    patch.set_facecolor(color)
plt.setp(box["medians"], color="black", linewidth=1)
plt.setp(box["whiskers"], color="black", linestyle="--", linewidth=0.5)
plt.setp(box["caps"], color="black", linewidth=0.5)
plt.ylabel("<TPM> (non zero)", fontsize=15)
plt.xlabel("Tissue gene classes", fontsize=15)
plt.xticks(fontsize=18, fontweight='bold')
plt.title('U genes development', fontsize=15, fontweight='bold')
plt.yticks(fontsize=13)
plt.yscale('log')
plt.savefig(path_save_data+'non_zero_mean_exp_bulk_U.png', dpi=600, bbox_inches='tight')
plt.show()


print('KS test between <TPM> distributions')
ks_U_S, pvalue_U_S = stats.ks_2samp(mean_exp_U_ini_U_tis, mean_exp_U_ini_S_tis)
print('U vs S:', ks_U_S, pvalue_U_S)
ks_U_hS, pvalue_U_hS = stats.ks_2samp(mean_exp_U_ini_U_tis, mean_exp_U_ini_hS_tis)
print('U vs hS:', ks_U_hS, pvalue_U_hS)
ks_S_hS, pvalue_S_hS = stats.ks_2samp(mean_exp_U_ini_S_tis, mean_exp_U_ini_hS_tis)
print('S vs hS:', ks_S_hS, pvalue_S_hS)
print('median hS tissue:', np.median(mean_exp_U_ini_hS_tis))
print('median U tissue:', np.median(mean_exp_U_ini_U_tis))




#big figure of bulk mean expression dev vs tissue

#figure of gene age distrib
data=[mean_exp_U_ini_U_tis, mean_exp_U_ini_S_tis, mean_exp_U_ini_hS_tis, 
      mean_exp_S_ini_U_tis, mean_exp_S_ini_S_tis, mean_exp_S_ini_hS_tis, 
      mean_exp_hS_ini_U_tis, mean_exp_hS_ini_S_tis, mean_exp_hS_ini_hS_tis]
fig, ax = plt.subplots(figsize=(4,3), dpi=600)
box = ax.boxplot(data, labels=["U-Ut", "U-St", "U-hSt", 
                               "S-Ut", "S-St", "S-hSt",
                               "hS-Ut", "hS-St", "hS-hSt"],  widths=0.45, patch_artist=True, 
                 flierprops=dict(marker='o', color='gray', markersize=0.5), 
                 boxprops=dict(edgecolor="none"))
colors = [
    "#0D47A1",  # Azul profundo
    "#1976D2",  # Azul medio
    "#5C6BC0",  # Azul violáceo
    "#7E57C2",  # Violeta suave
    "#AB47BC",  # Púrpura
    "#EC407A",  # Fucsia fuerte
    "mediumvioletred",  # Rosa medio
    "hotpink",  # Rosa claro
    "pink" ]  # Rosa pálido

for patch, color in zip(box["boxes"], colors):
    patch.set_facecolor(color)

plt.setp(box["medians"], color="black", linewidth=1)
plt.setp(box["whiskers"], color="black", linestyle="--", linewidth=0.5)
plt.setp(box["caps"], color="black", linewidth=0.5)

# plt.yscale(log)
plt.ylabel("<TPM>", fontsize=15)
plt.xlabel("Gene classes", fontsize=15)
plt.xticks(fontsize=18, fontweight='bold', rotation=90)
plt.yticks(fontsize=13)

comparisons = [(0,1), (0,2), (1,2)]  
labels = ['****', '****', '***']   
y_max = max([max(d) for d in data])
h = 0.4 * y_max  
line_heights = [y_max + 1.0*y_max, y_max + 10.2*y_max, y_max + 50.4*y_max]
for i, ((x1, x2), label) in enumerate(zip(comparisons, labels)):
    y = line_heights[i]
    ax.plot([x1+1, x1+1, x2+1, x2+1], [y, y+h, y+h, y], lw=0.5, c='black')
    ax.text((x1 + x2 + 2) / 2, y + 0.01 * y_max, label, ha='center', fontsize=10, c='red')
    
comparisons = [(3,4), (3,5), (4,5)]  
labels = ['****', '****', '****']   
y_max = max([max(d) for d in data])
h = 0.4 * y_max  
line_heights = [y_max + 1.0*y_max, y_max + 10.2*y_max, y_max + 50.4*y_max]
for i, ((x1, x2), label) in enumerate(zip(comparisons, labels)):
    y = line_heights[i]
    ax.plot([x1+1, x1+1, x2+1, x2+1], [y, y+h, y+h, y], lw=0.5, c='black')
    ax.text((x1 + x2 + 2) / 2, y + 0.01 * y_max, label, ha='center', fontsize=10, c='red')
    
comparisons = [(6,7), (6,8), (7,8)]  
labels = ['', '***', '']   
y_max = max([max(d) for d in data])
h = 0.4 * y_max  
line_heights = [y_max + 1.0*y_max, y_max + 10.2*y_max, y_max + 50.4*y_max]
for i, ((x1, x2), label) in enumerate(zip(comparisons, labels)):
    y = line_heights[i]
    ax.plot([x1+1, x1+1, x2+1, x2+1], [y, y+h, y+h, y], lw=0.5, c='black')
    ax.text((x1 + x2 + 2) / 2, y + 0.01 * y_max, label, ha='center', fontsize=10, c='red')
    
plt.yscale('log')
plt.savefig(path_save_data+'mean_exp_subclasses.png', dpi=600, bbox_inches='tight')

plt.show()



#KS test between maximal exp stages
print('KS test between maximal exp distributions: U dev genes')
ks_U_S, pvalue_U_S = stats.ks_2samp(max_stage_U_ini_U_tis, max_stage_U_ini_S_tis)
print('U vs S:', ks_U_S, pvalue_U_S)
ks_U_hS, pvalue_U_hS = stats.ks_2samp(max_stage_U_ini_U_tis, max_stage_U_ini_hS_tis)
print('U vs hS:', ks_U_hS, pvalue_U_hS)
ks_S_hS, pvalue_S_hS = stats.ks_2samp(max_stage_U_ini_S_tis, max_stage_U_ini_hS_tis)
print('S vs hS:', ks_S_hS, pvalue_S_hS)

print('KS test between maximal exp distributions: S dev genes')
ks_U_S, pvalue_U_S = stats.ks_2samp(max_stage_S_ini_U_tis, max_stage_S_ini_S_tis)
print('U vs S:', ks_U_S, pvalue_U_S)
ks_U_hS, pvalue_U_hS = stats.ks_2samp(max_stage_S_ini_U_tis, max_stage_S_ini_hS_tis)
print('U vs hS:', ks_U_hS, pvalue_U_hS)
ks_S_hS, pvalue_S_hS = stats.ks_2samp(max_stage_S_ini_S_tis, max_stage_S_ini_hS_tis)
print('S vs hS:', ks_S_hS, pvalue_S_hS)

print('KS test between maximal exp distributions: hS dev genes')
ks_U_S, pvalue_U_S = stats.ks_2samp(max_stage_hS_ini_U_tis, max_stage_hS_ini_S_tis)
print('U vs S:', ks_U_S, pvalue_U_S)
ks_U_hS, pvalue_U_hS = stats.ks_2samp(max_stage_hS_ini_U_tis, max_stage_hS_ini_hS_tis)
print('U vs hS:', ks_U_hS, pvalue_U_hS)
ks_S_hS, pvalue_S_hS = stats.ks_2samp(max_stage_hS_ini_S_tis, max_stage_hS_ini_hS_tis)
print('S vs hS:', ks_S_hS, pvalue_S_hS)


#7.) Differences in tissue
def tissue_differences(gene_subset, gene_subset_title, gene_subset_label, color_subset):
    frac_genes_per_tissue=np.zeros(len(tissue))
    g, ind_tis, j=np.intersect1d(genes, gene_subset, return_indices=True)
    del g, j
    matrix_subset=matrix[ind_tis, :]
    print(gene_subset_title, len(gene_subset))
    for i in range(len(tissue)):
        frac_genes_per_tissue[i]=len(np.where(matrix_subset[:, i]>1)[0])/len(gene_subset)

    plt.figure(figsize=(9, 4), dpi=600)
    plt.bar(tissue, frac_genes_per_tissue, color=color_subset, edgecolor='black')
    plt.xticks(ticks=np.arange(0, len(tissue)), labels=tissue, rotation=90, fontsize=25)
    plt.yticks(fontsize=25)
    plt.ylabel('% genes', fontsize=30, fontweight='bold')
    plt.title(gene_subset_title, fontsize=30, fontweight='bold')
    # plt.ylabel('# genes (TPM>1)', fontsize=30, fontweight='bold')
    plt.savefig(path_save_data+'%s_per_tissue.png' %gene_subset_label, dpi=600, bbox_inches='tight')
    plt.show()

tissue_differences(U_ini_S_tis, 'U dev - S tissue', 'U_ini_S_tis', 'mediumturquoise')
tissue_differences(U_ini_hS_tis, 'U dev - hS tissue', 'U_ini_hS_tis', 'tomato')

tissue_differences(S_ini_S_tis, 'S dev - S tissue', 'S_ini_S_tis', 'mediumturquoise')
tissue_differences(S_ini_hS_tis, 'S dev - hS tissue', 'S_ini_hS_tis', 'tomato')

tissue_differences(hS_ini_S_tis, 'hS dev - S tissue', 'hS_ini_S_tis', 'mediumturquoise')
tissue_differences(hS_ini_hS_tis, 'hS dev - hS tissue', 'hS_ini_hS_tis', 'tomato')



def stat_test_tissue(gene_subset):
    wilcoxon_results = {}
    
    g, ind_tis, j=np.intersect1d(genes, gene_subset, return_indices=True)
    del g, j
    matrix_subset=matrix[ind_tis, :]
    matrix_rest=np.delete(matrix, ind_tis, axis=0)
    
    for i in range(len(tissue)):
        
        # stat, p_value=stats.ks_2samp(matrix_subset[:, i], matrix_rest[:, i])
        # stat, p_value = ttest_ind(matrix_subset[:, i], matrix_rest[:, i], equal_var=False)  # Welch's t-test
        stat, p_value = ranksums(matrix_subset[:, i], matrix_rest[:, i])
        
    
        wilcoxon_results[i] = {"Wilcoxon Statistic": stat, "P-Value": p_value}

        if (p_value<0.0001) & (np.median(matrix_subset[:, i])>np.median(matrix_rest[:, i])):
            print(tissue[i], i)
            print('my subset:', np.median(matrix_subset[:, i]))
            print('the rest subset:', np.median(matrix_rest[:, i]), '\n')

    
    

print('U_ini_S_tis')
stat_test_tissue(U_ini_S_tis)
print('U_ini_hS_tis')
stat_test_tissue(U_ini_hS_tis)

print('S_ini_S_tis')
stat_test_tissue(S_ini_S_tis)
print('S_ini_hS_tis')
stat_test_tissue(S_ini_hS_tis)

print('hS_ini_hS_tis')
stat_test_tissue(hS_ini_S_tis)
print('hS_ini_hS_tis')
stat_test_tissue(hS_ini_hS_tis)











