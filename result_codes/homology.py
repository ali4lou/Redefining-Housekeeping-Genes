# -*- coding: utf-8 -*-
"""
Created on Mon May 26 13:11:07 2025

Homology
"""

import numpy as np
import os
import matplotlib.pyplot as plt 
import pandas as pd
import matplotlib as mpl
from scipy.stats import mstats, kstest, ttest_ind, fisher_exact
from scipy import stats
from scipy.stats import hypergeom

def test_hipergeometrico(total_poblacion, exitos_totales, tamaño_muestra, exitos_observados, verbose=True):

    N = total_poblacion
    K = exitos_totales
    n = tamaño_muestra
    k = exitos_observados

    # P-valor (cola superior)
    p_valor = hypergeom.sf(k-1, N, K, n)
    
    print(p_valor)
    return p_valor



path='PATH_WHERE_YOU_KEEP_YOUR_DATASETS'


#0.) Zfin gene related to ENS gene
f=open(path+'zfin_to_ENS.txt', 'r')
txt = f.read()
zfin_to_ENS = txt.split('\n')
del txt, f
zfin_to_ENS=np.delete(zfin_to_ENS, len(zfin_to_ENS)-1)
zfin_to_ENS=np.array(zfin_to_ENS)

zfin_all = []
ENS_all=[]

for elemento in zfin_to_ENS:
    partes = elemento.split("\t")

    for i in range(len(partes)):
    
        if partes[i].startswith("ZDB"):
            zfin_all.append(partes[i])
        
        if partes[i].startswith("ENS"):
            ENS_all.append(partes[i])

ENS_all=np.array(ENS_all)
zfin_all=np.array(zfin_all)

del zfin_to_ENS, partes, elemento



#1.) Read the data
f=open(path+'genes_bulk.txt', 'r')
txt = f.read()
genes_bulk = txt.split('\n')
del txt, f
genes_bulk=np.delete(genes_bulk, len(genes_bulk)-1)
genes_bulk=np.array(genes_bulk)

f=open(path+'hS_id_bulk.txt', 'r')
txt = f.read()
hS_id_bulk = txt.split('\n')
del txt, f
hS_id_bulk=np.delete(hS_id_bulk, len(hS_id_bulk)-1)
hS_id_bulk=np.array(hS_id_bulk)

f=open(path+'S_id_bulk.txt', 'r')
txt = f.read()
S_id_bulk = txt.split('\n')
del txt, f
S_id_bulk=np.delete(S_id_bulk, len(S_id_bulk)-1)
S_id_bulk=np.array(S_id_bulk)

f=open(path+'U_id_bulk.txt', 'r')
txt = f.read()
U_id_bulk = txt.split('\n')
del txt, f
U_id_bulk=np.delete(U_id_bulk, len(U_id_bulk)-1)
U_id_bulk=np.array(U_id_bulk)



#bulk genes and matrix
f=open(path+'time_bulk.txt', 'r')
txt = f.read()
time_bulk = txt.split('\n')
del txt, f
time_bulk=np.delete(time_bulk, len(time_bulk)-1)
time_bulk=np.array(time_bulk)


mean_bulk_matrix= np.genfromtxt(path+'mean_bulk_matrix.txt', delimiter=',', dtype=None, encoding=None)


f=open(path+'genes_bulk.txt', 'r')
txt = f.read()
genes_bulk = txt.split('\n')
del txt, f
genes_bulk=np.delete(genes_bulk, len(genes_bulk)-1)
genes_bulk=np.array(genes_bulk)



#9.) Paralog 

#9.1.) Ortology association
f=open(path+'paralog.txt', 'r')
txt = f.read()
paralog_data = txt.split('\n')
del txt, f
paralog_data=np.delete(paralog_data, len(paralog_data)-1)
paralog_data=np.array(paralog_data)

gene_with_paralog = []
paralog_mate=[]

for elemento in paralog_data:
    partes = elemento.split(",")
    if partes[2]=='within_species_paralog':
        gene_with_paralog.append(partes[0])
        paralog_mate.append(partes[1])

gene_with_paralog=np.array(gene_with_paralog)
paralog_mate=np.array(paralog_mate)
gene_with_paralog_unique, n_paralog_times=np.unique(gene_with_paralog, return_counts=True)


print('total paralog: ', len(gene_with_paralog))


#9.1.) We search the diseases associated to each gene class
common_genes_U_par, ind_U_par, ind_U_bulk=np.intersect1d(gene_with_paralog_unique, U_id_bulk, return_indices=True)
print('U genes:', len(U_id_bulk), ' - ', 'par U:', len(common_genes_U_par))
n_times_U_par=np.sum(n_paralog_times[ind_U_par])
n_paralog_all_U=n_paralog_times[ind_U_par]
print(n_times_U_par)

common_genes_S_par, ind_S_par, ind_S_bulk=np.intersect1d(gene_with_paralog_unique, S_id_bulk, return_indices=True)
print('S genes:', len(S_id_bulk), ' - ', 'par S:', len(common_genes_S_par))
n_times_S_par=np.sum(n_paralog_times[ind_S_par])
n_paralog_all_S=n_paralog_times[ind_S_par]
print(n_times_S_par)

common_genes_hS_par, ind_hS_par, ind_hS_bulk=np.intersect1d(gene_with_paralog_unique, hS_id_bulk, return_indices=True)
print('hS genes:', len(hS_id_bulk), ' - ', 'par hS:', len(common_genes_hS_par))
n_times_hS_par=np.sum(n_paralog_times[ind_hS_par])
n_paralog_all_hS=n_paralog_times[ind_hS_par]
print(n_times_hS_par)

print('p-value')
print('U genes', len(common_genes_U_par), len(U_id_bulk))
p_valor_U = test_hipergeometrico(
    total_poblacion=len(genes_bulk),
    exitos_totales=len(gene_with_paralog_unique),
    tamaño_muestra=len(U_id_bulk),
    exitos_observados=len(common_genes_U_par)
)

print('S genes', len(common_genes_S_par), len(S_id_bulk))
p_valor_S = test_hipergeometrico(
    total_poblacion=len(genes_bulk),
    exitos_totales=len(gene_with_paralog_unique),
    tamaño_muestra=len(S_id_bulk),
    exitos_observados=len(common_genes_S_par)
)

print('hS genes', len(common_genes_hS_par), len(hS_id_bulk))
p_valor_hS = test_hipergeometrico(
    total_poblacion=len(genes_bulk),
    exitos_totales=len(gene_with_paralog_unique),
    tamaño_muestra=len(hS_id_bulk),
    exitos_observados=len(common_genes_hS_par)
)




#==========================================================
#==========================================================
#we check if each couple are in the same class
same_class = 0
other = 0
s_class_par = 0
hs_class_par = 0
for g, p in zip(gene_with_paralog, paralog_mate):
    if np.any(U_id_bulk == g):
        if np.any(U_id_bulk == p):
            same_class += 1
        else:
            other += 1
            if np.any(S_id_bulk == p):
                s_class_par += 1
            elif np.any(hS_id_bulk == p):
                hs_class_par += 1
 
    
print('same class:', same_class/n_times_U_par, same_class)
print('other class:', other/n_times_U_par, other)
print('S class', s_class_par/other, s_class_par)
print('hs class', hs_class_par/other, hs_class_par)
 
#=======================================================================
print('U genes')
U_set = set(U_id_bulk)  # convierte a set para búsquedas rápidas

genes_bulk_set = set(genes_bulk) 

found=[g for g in gene_with_paralog if g in U_set]
print(len(found))

in_U_and_any = [(g, p) for g, p in zip(gene_with_paralog, paralog_mate) if g in U_set and p in genes_bulk_set]
print('U genes with paralogs in the dataset', len(in_U_and_any), len(in_U_and_any)/len(found))


both_in_U = [(g, p) for g, p in zip(gene_with_paralog, paralog_mate) if g in U_set and p in U_set]
print('U-U', len(both_in_U), len(both_in_U)/len(in_U_and_any))

S_set = set(S_id_bulk)  # convierte a set para búsquedas rápidas
in_U_and_S = [(g, p) for g, p in zip(gene_with_paralog, paralog_mate) if g in U_set and p in S_set]
print('U-S', len(in_U_and_S), len(in_U_and_S)/len(in_U_and_any))

hS_set = set(hS_id_bulk)  # convierte a set para búsquedas rápidas
in_U_and_hS = [(g, p) for g, p in zip(gene_with_paralog, paralog_mate) if g in U_set and p in hS_set]
print('U-hS', len(in_U_and_hS), len(in_U_and_hS)/len(in_U_and_any))

all_genes = U_set | S_set | hS_set
not_in_classes = genes_bulk_set - all_genes

in_U_and_non_exp = [(g, p) for g, p in zip(gene_with_paralog, paralog_mate) if g in U_set and p in not_in_classes]
print('non exp', len(in_U_and_non_exp), len(in_U_and_non_exp)/len(in_U_and_any))
#=======================================================================


            

same_class=0
other=0
u_class_par=0
hs_class_par=0
for i in range(len(gene_with_paralog)):
    ind_g=np.where(S_id_bulk==gene_with_paralog[i])[0]
    if len(ind_g)>0:
        ind_couple=np.where(S_id_bulk==paralog_mate[i])[0]
        if len(ind_couple)>0:
            same_class=same_class+1
        else:
            other=other+1
            if len(np.where(U_id_bulk==paralog_mate[i])[0])>0:
                u_class_par=u_class_par+1
            if len(np.where(hS_id_bulk==paralog_mate[i])[0])>0:
                hs_class_par=hs_class_par+1
            
print('same class:', same_class/n_times_S_par, same_class)
print('other class:', other/n_times_S_par, other)
print('u class', u_class_par/other, u_class_par)
print('hs class', hs_class_par/other, hs_class_par)

#=======================================================================
print('S genes')
S_set = set(S_id_bulk)  # convierte a set para búsquedas rápidas

found=[g for g in gene_with_paralog if g in S_set]
print(len(found))

in_S_and_any = [(g, p) for g, p in zip(gene_with_paralog, paralog_mate) if g in S_set and p in genes_bulk_set]
print('S genes with paralog in the dataset', len(in_S_and_any), len(in_S_and_any)/len(found))

both_in_S = [(g, p) for g, p in zip(gene_with_paralog, paralog_mate) if g in S_set and p in S_set]
print('S-S', len(both_in_S), len(both_in_S)/len(in_S_and_any))

U_set = set(U_id_bulk)  # convierte a set para búsquedas rápidas
in_S_and_U = [(g, p) for g, p in zip(gene_with_paralog, paralog_mate) if g in S_set and p in U_set]
print('S-U', len(in_S_and_U), len(in_S_and_U)/len(in_S_and_any))

hS_set = set(hS_id_bulk)  # convierte a set para búsquedas rápidas
in_S_and_hS = [(g, p) for g, p in zip(gene_with_paralog, paralog_mate) if g in S_set and p in hS_set]
print('S-hS', len(in_S_and_hS), len(in_S_and_hS)/len(in_S_and_any))

all_genes = U_set | S_set | hS_set
genes_bulk_set = set(genes_bulk) 
not_in_classes = genes_bulk_set - all_genes

in_S_and_non_exp = [(g, p) for g, p in zip(gene_with_paralog, paralog_mate) if g in S_set and p in not_in_classes]
print('non exp', len(in_S_and_non_exp), len(in_S_and_non_exp)/len(in_S_and_any))
#=======================================================================



same_class=0
other=0
u_class_par=0
s_class_par=0
for i in range(len(gene_with_paralog)):
    ind_g=np.where(hS_id_bulk==gene_with_paralog[i])[0]
    if len(ind_g)>0:
        ind_couple=np.where(hS_id_bulk==paralog_mate[i])[0]
        if len(ind_couple)>0:
            same_class=same_class+1
        else:
            other=other+1
            if len(np.where(U_id_bulk==paralog_mate[i])[0])>0:
                u_class_par=u_class_par+1
            if len(np.where(S_id_bulk==paralog_mate[i])[0])>0:
                s_class_par=s_class_par+1
            
print('same class:', same_class/n_times_hS_par, same_class)
print('other class:', other/n_times_hS_par, other)
print('S class', s_class_par/other, s_class_par)
print('u class', u_class_par/other, u_class_par)


#=======================================================================
print('hS genes')
hS_set = set(hS_id_bulk)  # convierte a set para búsquedas rápidas

found=[g for g in gene_with_paralog if g in hS_set]
print(len(found))

in_hS_and_any = [(g, p) for g, p in zip(gene_with_paralog, paralog_mate) if g in hS_set and p in genes_bulk_set]
print('hS genes with paralog in the dataset', len(in_hS_and_any), len(in_hS_and_any)/len(found))

both_in_hS = [(g, p) for g, p in zip(gene_with_paralog, paralog_mate) if g in hS_set and p in hS_set]
print('hS-hS', len(both_in_hS), len(both_in_hS)/len(in_hS_and_any))

U_set = set(U_id_bulk)  # convierte a set para búsquedas rápidas
in_hS_and_U = [(g, p) for g, p in zip(gene_with_paralog, paralog_mate) if g in hS_set and p in U_set]
print('hS-U', len(in_hS_and_U), len(in_hS_and_U)/len(in_hS_and_any))

S_set = set(S_id_bulk)  # convierte a set para búsquedas rápidas
in_hS_and_S = [(g, p) for g, p in zip(gene_with_paralog, paralog_mate) if g in hS_set and p in S_set]
print('hS-S', len(in_hS_and_S), len(in_hS_and_S)/len(in_hS_and_any))

all_genes = U_set | S_set | hS_set
genes_bulk_set = set(genes_bulk) 
not_in_classes = genes_bulk_set - all_genes

in_hS_and_non_exp = [(g, p) for g, p in zip(gene_with_paralog, paralog_mate) if g in hS_set and p in not_in_classes]
print('non exp', len(in_hS_and_non_exp), len(in_hS_and_non_exp)/len(in_hS_and_any))
#=======================================================================



#==========================================================
#==========================================================


#figure
frac_genes=np.array([len(common_genes_U_par)/len(U_id_bulk), len(common_genes_S_par)/len(S_id_bulk), len(common_genes_hS_par)/len(hS_id_bulk)])
groups=['U', 'S', 'hS']
color=['royalblue', 'mediumturquoise', 'tomato']

plt.figure(figsize=(4, 3), dpi=600)
plt.bar(groups, frac_genes, color=color)
plt.xticks(ticks=np.arange(len(groups)), labels=groups, fontsize=16, fontweight='bold')
plt.yticks(fontsize=14)
for i in range(len(groups)):
    plt.text(i-0.18, frac_genes[i]-0.1, np.round(frac_genes[i], 2), fontsize=12, fontweight='bold', color='white')
plt.ylabel('Fraction of\nparalogs', fontsize=18, fontweight='bold')
plt.savefig(path+'frac_paralogs_per_gene_class.png', dpi=600, bbox_inches='tight')
plt.show()

print(frac_genes)


#===========================================================
#stage of maximal expression
#index of bulk genes
common_genes_U_par, ind_U_par, ind_U_bulk=np.intersect1d(gene_with_paralog_unique, U_id_bulk, return_indices=True)


g, ind_U_par_with_dev, ind_U_let=np.intersect1d(common_genes_U_par, genes_bulk, return_indices=True)
g, ind_S_par_with_dev, ind_S_let=np.intersect1d(common_genes_S_par, genes_bulk, return_indices=True)
g, ind_hS_par_with_dev, ind_hS_let=np.intersect1d(common_genes_hS_par, genes_bulk, return_indices=True)
   
U_matrix_par=mean_bulk_matrix[ind_U_let, :]
S_matrix_par=mean_bulk_matrix[ind_S_let, :]
hS_matrix_par=mean_bulk_matrix[ind_hS_let, :]

#n paralog
n_par_U_with_dev=n_paralog_all_U[ind_U_par_with_dev]
n_par_S_with_dev=n_paralog_all_S[ind_S_par_with_dev]
n_par_hS_with_dev=n_paralog_all_hS[ind_hS_par_with_dev]


#For each gene we find the stage of maximal exp 
stage_of_max_exp_U=np.zeros(len(U_matrix_par[:, 0]))
for i in range(len(U_matrix_par[:, 0])):
    ind_max=np.argmax(U_matrix_par[i, :])
    stage_of_max_exp_U[i]=ind_max
    
stage_of_max_exp_S=np.zeros(len(S_matrix_par[:, 0]))
for i in range(len(S_matrix_par[:, 0])):
    ind_max=np.argmax(S_matrix_par[i, :])
    stage_of_max_exp_S[i]=ind_max
    
stage_of_max_exp_hS=np.zeros(len(hS_matrix_par[:, 0]))
for i in range(len(hS_matrix_par[:, 0])):
    ind_max=np.argmax(hS_matrix_par[i, :])
    stage_of_max_exp_hS[i]=ind_max



#mixed fig
time_unique=np.array(['0hpf', '0.75hpf', '2.25hpf', '3hpf', '4.3hpf', '5.25hpf', '6hpf', '8hpf', '10.3hpf', '16hpf', '19hpf', '24hpf', '30hpf', '36hpf', '2dpf', '3dpf', '4dpf', '5dpf'])

act_stages_par_U=np.zeros(len(time_unique))
act_stages_par_S=np.zeros(len(time_unique))
act_stages_par_hS=np.zeros(len(time_unique))
for j in (stage_of_max_exp_U):
    act_stages_par_U[int(j)]=act_stages_par_U[int(j)]+1
for j in (stage_of_max_exp_S):
    act_stages_par_S[int(j)]=act_stages_par_S[int(j)]+1
for j in (stage_of_max_exp_hS):
    act_stages_par_hS[int(j)]=act_stages_par_hS[int(j)]+1

act_stages_par_U=act_stages_par_U/len(U_matrix_par[:, 0])
act_stages_par_S=act_stages_par_S/len(S_matrix_par[:, 0])
act_stages_par_hS=act_stages_par_hS/len(hS_matrix_par[:, 0])



#mixed fig
time_unique=np.array(['0hpf', '0.75hpf', '2.25hpf', '3hpf', '4.3hpf', '5.25hpf', '6hpf', '8hpf', '10.3hpf', '16hpf', '19hpf', '24hpf', '30hpf', '36hpf', '2dpf', '3dpf', '4dpf', '5dpf'])

width = 0.25  # ancho de cada barra
x = np.arange(len(time_unique))
plt.figure(figsize=(11, 6), dpi=600)
plt.bar(x - width, act_stages_par_U * 100, width=width, label='U', color='royalblue')
plt.bar(x,         act_stages_par_S * 100, width=width, label='S', color='mediumturquoise')
plt.bar(x + width, act_stages_par_hS * 100, width=width, label='hS', color='tomato')

plt.xticks(ticks=x, labels=time_unique, rotation=90, fontsize=25)
plt.yticks(fontsize=25)
plt.xlabel('Stages', fontsize=30)
plt.ylabel('% paralog genes\nwith max exp', fontsize=25, fontweight='bold')
plt.legend(fontsize=20)

plt.tight_layout()
plt.savefig(path + 'frac_paralog_genes_stage_maximal_exp.png', dpi=600, bbox_inches='tight')
plt.show()




#boxplot paraglog distrib
data = [n_paralog_all_U, n_paralog_all_S, n_paralog_all_hS]
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
plt.ylabel("# of associated\nparalogs per gene", fontsize=15)
plt.xlabel("Gene classes", fontsize=15)
plt.xticks(fontsize=18, fontweight='bold')
plt.yticks(fontsize=13)
plt.yscale('log')
comparisons = [(0,1), (0,2), (1,2)]  
labels = ['****', '****', '****']   
y_max = max([max(d) for d in data])
h = 0.4 * y_max  
line_heights = [y_max + 1.0*y_max, y_max + 2.2*y_max, y_max + 3.4*y_max]
for i, ((x1, x2), label) in enumerate(zip(comparisons, labels)):
    y = line_heights[i]
    ax.plot([x1+1, x1+1, x2+1, x2+1], [y, y+h, y+h, y], lw=0.5, c='black')
    ax.text((x1 + x2 + 2) / 2, y + 0.1 * y_max, label, ha='center', fontsize=10, c='red')
plt.savefig(path + 'paralog_per_gene_distrib_significance.png', dpi=600, bbox_inches='tight')
plt.show()




ks_U_S, pvalue_U_S = stats.ks_2samp(n_paralog_times[ind_U_par], n_paralog_times[ind_S_par])
print('U vs S:', ks_U_S, pvalue_U_S)
ks_U_hS, pvalue_U_hS = stats.ks_2samp(n_paralog_times[ind_U_par], n_paralog_times[ind_hS_par])
print('U vs hS:', ks_U_hS, pvalue_U_hS)
ks_S_hS, pvalue_S_hS = stats.ks_2samp(n_paralog_times[ind_S_par], n_paralog_times[ind_hS_par])
print('S vs hS:', ks_S_hS, pvalue_S_hS)
    
print(np.median(n_paralog_times[ind_U_par]), np.median(n_paralog_times[ind_S_par]), np.median(n_paralog_times[ind_hS_par]))


#sliding window n paralogs per gene class vs peak gene expression
sw=400
#U genes
stage_of_max_exp_sorted_U=np.sort(stage_of_max_exp_U)
index_stage_sort_U=np.argsort(stage_of_max_exp_U)
n_par_U_with_dev_sorted=np.zeros(len(stage_of_max_exp_U))
for i in range(len(index_stage_sort_U)):
    n_par_U_with_dev_sorted[i]=n_par_U_with_dev[int(index_stage_sort_U[i])]

serie1U = pd.Series(stage_of_max_exp_sorted_U)
serie2U = pd.Series(n_par_U_with_dev_sorted)

sw_stage_U = serie1U.rolling(window=sw, center=False).mean()
sw_pleio_U = serie2U.rolling(window=sw, center=False).mean()

#S genes
stage_of_max_exp_sorted_S=np.sort(stage_of_max_exp_S)
index_stage_sort_S=np.argsort(stage_of_max_exp_S)
n_par_S_with_dev_sorted=np.zeros(len(stage_of_max_exp_S))
for i in range(len(index_stage_sort_S)):
    n_par_S_with_dev_sorted[i]=n_par_S_with_dev[int(index_stage_sort_S[i])]

serie1S = pd.Series(stage_of_max_exp_sorted_S)
serie2S = pd.Series(n_par_S_with_dev_sorted)

sw_stage_S = serie1S.rolling(window=sw, center=False).mean()
sw_pleio_S = serie2S.rolling(window=sw, center=False).mean()

#hS genes
stage_of_max_exp_sorted_hS=np.sort(stage_of_max_exp_hS)
index_stage_sort_hS=np.argsort(stage_of_max_exp_hS)
n_par_hS_with_dev_sorted=np.zeros(len(stage_of_max_exp_hS))
for i in range(len(index_stage_sort_hS)):
    n_par_hS_with_dev_sorted[i]=n_par_hS_with_dev[int(index_stage_sort_hS[i])]

serie1hS = pd.Series(stage_of_max_exp_sorted_hS)
serie2hS = pd.Series(n_par_hS_with_dev_sorted)

sw_stage_hS = serie1hS.rolling(window=sw, center=False).mean()
sw_pleio_hS = serie2hS.rolling(window=sw, center=False).mean()



# # figure
plt.figure(figsize=(5, 3), dpi=600)
plt.scatter(sw_stage_U, sw_pleio_U, s=1, color='royalblue', label='U')
plt.scatter(sw_stage_S, sw_pleio_S, s=1, color='mediumturquoise', label='S')
plt.scatter(sw_stage_hS, sw_pleio_hS, s=1, color='tomato', label='hS')
plt.xlabel('Stage max. expression', fontsize=13, fontweight='bold')
plt.ylabel('# paralogs', fontsize=13, fontweight='bold')
plt.xticks(np.arange(len(time_bulk)), time_bulk, fontsize=13, rotation=90)
plt.yticks(fontsize=13)
plt.legend(markerscale=5)
plt.savefig(path+'n_paralog_vs_max_stage_exp_sw_%d.png' %sw, dpi=600, bbox_inches='tight')
plt.show()




# #========================================================================
# #========================================================================
#10.) Human ortology


#10.1.) Ortology association
f=open(path+'human_orthos.txt', 'r')
txt = f.read()
human_orthos = txt.split('\n')
del txt, f
human_orthos=np.delete(human_orthos, len(human_orthos)-1)
human_orthos=np.array(human_orthos)

zfin_ort_ini = []
human_ort=[]

for elemento in human_orthos:
    partes = elemento.split("\t")

    zfin_ort_ini.append(partes[0])
    human_ort.append(partes[3])

zfin_ort_ini=np.array(zfin_ort_ini)
human_ort=np.array(human_ort)


unique_zfin_ort=np.unique(zfin_ort_ini)
human_ort_per_unique_zfin=[]
n_associated_human_ort=[]
#10.2.) Unique associations
for i in range(len(unique_zfin_ort)):
    ind_ort=np.where(zfin_ort_ini==unique_zfin_ort[i])[0]
    id_ort_hum=human_ort[ind_ort]
    human_ort_per_unique_zfin.append(np.unique(id_ort_hum))
    n_associated_human_ort.append(len(np.unique(id_ort_hum)))

n_associated_human_ort=np.array(n_associated_human_ort)



zfin_ort_common, ind_all, ind_zfin_to_human_ort=np.intersect1d(zfin_all, unique_zfin_ort, return_indices=True)

print('ort ENS:', len(zfin_ort_common))

final_genes_ort=ENS_all[ind_all]
n_associated_human_ort_final=n_associated_human_ort[ind_zfin_to_human_ort]


print('n human genes with ort:', np.sum(n_associated_human_ort_final))

final_human_ort_per_unique_zfin=[]
for i in range(len(ind_zfin_to_human_ort)):
    final_human_ort_per_unique_zfin.append(human_ort_per_unique_zfin[int(ind_zfin_to_human_ort[i])])

del zfin_ort_common, ind_all, ind_zfin_to_human_ort
del n_associated_human_ort

total_human_genes_ort=np.sum(n_associated_human_ort_final)


#9.1.) We search the frac of ortho
common_genes_U_ort, ind_U_ort, ind_U_bulk=np.intersect1d(final_genes_ort, U_id_bulk, return_indices=True)
print('U genes:', len(U_id_bulk), ' - ', 'ort U:', len(common_genes_U_ort))
n_U_ort=np.sum(n_associated_human_ort_final[ind_U_ort])
print(n_U_ort/total_human_genes_ort)

common_genes_S_ort, ind_S_ort, ind_S_bulk=np.intersect1d(final_genes_ort, S_id_bulk, return_indices=True)
print('S genes:', len(S_id_bulk), ' - ', 'ort S:', len(common_genes_S_ort))
n_S_ort=np.sum(n_associated_human_ort_final[ind_S_ort])
print(n_S_ort/total_human_genes_ort)

common_genes_hS_ort, ind_hS_ort, ind_hS_bulk=np.intersect1d(final_genes_ort, hS_id_bulk, return_indices=True)
print('hS genes:', len(hS_id_bulk), ' - ', 'ort hS:', len(common_genes_hS_ort))
n_hS_ort=np.sum(n_associated_human_ort_final[ind_hS_ort])
print(n_hS_ort/total_human_genes_ort)
        
print((n_U_ort + n_S_ort + n_hS_ort)/(len(U_id_bulk)+len(S_id_bulk)+len(hS_id_bulk)))


print(len(common_genes_U_ort)/len(U_id_bulk), len(common_genes_S_ort)/len(S_id_bulk), len(common_genes_hS_ort)/len(hS_id_bulk))
#figure
frac_genes=np.array([len(common_genes_U_ort)/len(U_id_bulk), len(common_genes_S_ort)/len(S_id_bulk), len(common_genes_hS_ort)/len(hS_id_bulk)])
groups=['U', 'S', 'hS']
color=['royalblue', 'mediumturquoise', 'tomato']

plt.figure(figsize=(4, 3), dpi=600)
plt.bar(groups, frac_genes, color=color)
plt.xticks(ticks=np.arange(len(groups)), labels=groups, fontsize=16, fontweight='bold')
plt.yticks(fontsize=14)
for i in range(len(groups)):
    plt.text(i-0.18, frac_genes[i]-0.1, np.round(frac_genes[i], 2), fontsize=12, fontweight='bold', color='white')
if "U" in groups:
    idx_U = groups.index("U")  # posición de la barra 'S'
    plt.text(idx_U, frac_genes[idx_U]-0.05, "*",
             ha='center', va='bottom', fontsize=18, color='red')
plt.ylabel('Fraction of\nhuman orthologs', fontsize=16, fontweight='bold')
plt.savefig(path+'frac_human_ortho_per_gene_class.png', dpi=600, bbox_inches='tight')
plt.show()


print('p-value')
print('U genes', len(common_genes_U_ort), len(U_id_bulk))
p_valor_U = test_hipergeometrico(
    total_poblacion=len(genes_bulk),
    exitos_totales=len(final_genes_ort),
    tamaño_muestra=len(U_id_bulk),
    exitos_observados=len(common_genes_U_ort)
)

print('S genes', len(common_genes_S_ort), len(S_id_bulk))
p_valor_S = test_hipergeometrico(
    total_poblacion=len(genes_bulk),
    exitos_totales=len(final_genes_ort),
    tamaño_muestra=len(S_id_bulk),
    exitos_observados=len(common_genes_S_ort)
)

print('hS genes', len(common_genes_hS_ort), len(hS_id_bulk))
p_valor_hS = test_hipergeometrico(
    total_poblacion=len(genes_bulk),
    exitos_totales=len(final_genes_ort),
    tamaño_muestra=len(hS_id_bulk),
    exitos_observados=len(common_genes_hS_ort)
)

#================================================================
#================================================================
#stage of maximal expression
#index of bulk genes
g, i, ind_U_let=np.intersect1d(common_genes_U_ort, genes_bulk, return_indices=True)
g, i, ind_S_let=np.intersect1d(common_genes_S_ort, genes_bulk, return_indices=True)
g, i, ind_hS_let=np.intersect1d(common_genes_hS_ort, genes_bulk, return_indices=True)
   
U_matrix_ort=mean_bulk_matrix[ind_U_let, :]
S_matrix_ort=mean_bulk_matrix[ind_S_let, :]
hS_matrix_ort=mean_bulk_matrix[ind_hS_let, :]



#For each gene we find the stage of maximal exp 
stage_of_max_exp_U=np.zeros(len(U_matrix_ort[:, 0]))
for i in range(len(U_matrix_ort[:, 0])):
    ind_max=np.argmax(U_matrix_ort[i, :])
    stage_of_max_exp_U[i]=ind_max
    
stage_of_max_exp_S=np.zeros(len(S_matrix_ort[:, 0]))
for i in range(len(S_matrix_ort[:, 0])):
    ind_max=np.argmax(S_matrix_ort[i, :])
    stage_of_max_exp_S[i]=ind_max
    
stage_of_max_exp_hS=np.zeros(len(hS_matrix_ort[:, 0]))
for i in range(len(hS_matrix_ort[:, 0])):
    ind_max=np.argmax(hS_matrix_ort[i, :])
    stage_of_max_exp_hS[i]=ind_max



#mixed fig
time_unique=np.array(['0hpf', '0.75hpf', '2.25hpf', '3hpf', '4.3hpf', '5.25hpf', '6hpf', '8hpf', '10.3hpf', '16hpf', '19hpf', '24hpf', '30hpf', '36hpf', '2dpf', '3dpf', '4dpf', '5dpf'])

act_stages_par_U=np.zeros(len(time_unique))
act_stages_par_S=np.zeros(len(time_unique))
act_stages_par_hS=np.zeros(len(time_unique))
for j in (stage_of_max_exp_U):
    act_stages_par_U[int(j)]=act_stages_par_U[int(j)]+1
for j in (stage_of_max_exp_S):
    act_stages_par_S[int(j)]=act_stages_par_S[int(j)]+1
for j in (stage_of_max_exp_hS):
    act_stages_par_hS[int(j)]=act_stages_par_hS[int(j)]+1

act_stages_par_U=act_stages_par_U/len(U_matrix_par[:, 0])
act_stages_par_S=act_stages_par_S/len(S_matrix_par[:, 0])
act_stages_par_hS=act_stages_par_hS/len(hS_matrix_par[:, 0])



#mixed fig
time_unique=np.array(['0hpf', '0.75hpf', '2.25hpf', '3hpf', '4.3hpf', '5.25hpf', '6hpf', '8hpf', '10.3hpf', '16hpf', '19hpf', '24hpf', '30hpf', '36hpf', '2dpf', '3dpf', '4dpf', '5dpf'])

width = 0.25  # ancho de cada barra
x = np.arange(len(time_unique))
plt.figure(figsize=(11, 6), dpi=600)
plt.bar(x - width, act_stages_par_U * 100, width=width, label='U', color='royalblue')
plt.bar(x,         act_stages_par_S * 100, width=width, label='S', color='mediumturquoise')
plt.bar(x + width, act_stages_par_hS * 100, width=width, label='hS', color='tomato')

plt.xticks(ticks=x, labels=time_unique, rotation=90, fontsize=25)
plt.yticks(fontsize=25)
plt.xlabel('Stages', fontsize=30)
plt.ylabel('% orthologs genes\nwith max exp', fontsize=25, fontweight='bold')
plt.legend(fontsize=20)

plt.tight_layout()
plt.savefig(path + 'frac_orthologs_genes_stage_maximal_exp.png', dpi=600, bbox_inches='tight')
plt.show()



#=========================================================================================================
#=========================================================================================================
#Subclasses Analysis
#9.1.) Paralogs
f=open(path+'paralog.txt', 'r')
txt = f.read()
paralog_data = txt.split('\n')
del txt, f
paralog_data=np.delete(paralog_data, len(paralog_data)-1)
paralog_data=np.array(paralog_data)

gene_with_paralog = []

for elemento in paralog_data:
    partes = elemento.split(",")
    if partes[2]=='within_species_paralog':
        gene_with_paralog.append(partes[0])

gene_with_paralog=np.array(gene_with_paralog)
gene_with_paralog_unique, n_paralog_times=np.unique(gene_with_paralog, return_counts=True)


print('total paralog: ', len(gene_with_paralog))



def paralogs(gene_analyse):
    common_genes_par, ind_par, ind_bulk=np.intersect1d(gene_with_paralog_unique, gene_analyse, return_indices=True)
    print(np.median(n_paralog_times[ind_par]))
    
    #index of bulk genes
    g, i, ind_let=np.intersect1d(common_genes_par, genes_bulk, return_indices=True)
    par_matrix=mean_bulk_matrix[ind_let, :]
    
    act_stages_par=np.zeros(len(mean_bulk_matrix[0, :]))
    
    if len(ind_let)>0:
        for i in range(len(mean_bulk_matrix[0, :])):
            act_stages_par[i]=len(np.where(par_matrix[:, i]>1)[0])/len(ind_let)
       
    return len(common_genes_par)/len(gene_analyse), n_paralog_times[ind_par], act_stages_par


#3.) We read the subclasses
label=['U_dev_U_tis', 'U_dev_S_tis', 'U_dev_hS_tis', 
       'S_dev_U_tis', 'S_dev_S_tis', 'S_dev_hS_tis',
       'hS_dev_U_tis', 'hS_dev_S_tis', 'hS_dev_hS_tis']
for k in (label):
    genes_analyze=np.loadtxt(path+"%s.txt" %k, dtype=str)
    genes_label=k
    globals()[f"{k}_frac_paralogs"], globals()[f"{k}_n_paralogs"], globals()[f"{k}_act_stages_par"]=paralogs(genes_analyze)
    


#one common figure
frac_genes_subclass=np.array([U_dev_U_tis_frac_paralogs, U_dev_S_tis_frac_paralogs, U_dev_hS_tis_frac_paralogs, 
                              S_dev_U_tis_frac_paralogs, S_dev_S_tis_frac_paralogs, S_dev_hS_tis_frac_paralogs, 
                              hS_dev_U_tis_frac_paralogs, hS_dev_S_tis_frac_paralogs, hS_dev_hS_tis_frac_paralogs])
groups=['U-Ut', 'U-St', 'U-hSt',
        'S-Ut', 'S-St', 'S-hSt', 
        'hS-Ut', 'hS-St', 'hS-hSt']
color=["#0D47A1", "#1976D2", "#5C6BC0", 
       "#7E57C2", "#AB47BC", "#EC407A", 
       "#F06292", "#F48FB1", "#F8BBD0"]
plt.figure(figsize=(4, 3), dpi=600)
plt.bar(groups, frac_genes_subclass, color=color)
plt.xticks(ticks=np.arange(len(groups)), labels=groups, fontsize=16, fontweight='bold', rotation=90)
plt.yticks(fontsize=14)
for i in range(len(groups)):
    plt.text(i-0.4, frac_genes_subclass[i], np.round(frac_genes_subclass[i], 2), fontsize=8, fontweight='bold', color='black')
plt.ylabel('% paralogs', fontsize=18, fontweight='bold')
plt.savefig(path+'frac_paralogs_per_gene_subclass.png', dpi=600, bbox_inches='tight')
plt.show()


#paralog distrib
data=[U_dev_U_tis_n_paralogs, U_dev_S_tis_n_paralogs, U_dev_hS_tis_n_paralogs, 
      S_dev_U_tis_n_paralogs, S_dev_S_tis_n_paralogs, S_dev_hS_tis_n_paralogs, 
      hS_dev_U_tis_n_paralogs, hS_dev_S_tis_n_paralogs, hS_dev_hS_tis_n_paralogs]
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
    "#F06292",  # Rosa medio
    "#F48FB1",  # Rosa claro
    "#F8BBD0" ]  # Rosa pálido

for patch, color in zip(box["boxes"], colors):
    patch.set_facecolor(color)

plt.setp(box["medians"], color="black", linewidth=1)
plt.setp(box["whiskers"], color="black", linestyle="--", linewidth=0.5)
plt.setp(box["caps"], color="black", linewidth=0.5)

# plt.yscale(log)
plt.ylabel("# paralogs", fontsize=15)
plt.xlabel("Gene classes", fontsize=15)
plt.xticks(fontsize=18, fontweight='bold', rotation=90)
plt.yticks(fontsize=13)

plt.yscale('log')
plt.savefig(path+'n_paralogs_boxplot_subclasses.png', dpi=600, bbox_inches='tight')

plt.show()




#4.) KS test
from scipy import stats
ks_U_S, pvalue_U_S = stats.ks_2samp(U_dev_U_tis_n_paralogs, U_dev_S_tis_n_paralogs)
print('U vs S:', ks_U_S, pvalue_U_S)
ks_U_hS, pvalue_U_hS = stats.ks_2samp(U_dev_U_tis_n_paralogs, U_dev_hS_tis_n_paralogs)
print('U vs hS:', ks_U_hS, pvalue_U_hS)
ks_S_hS, pvalue_S_hS = stats.ks_2samp(U_dev_S_tis_n_paralogs, U_dev_hS_tis_n_paralogs)
print('S vs hS:', ks_S_hS, pvalue_S_hS)
print('significant difference between all the age distrib')
    
from scipy import stats
ks_U_S, pvalue_U_S = stats.ks_2samp(S_dev_U_tis_n_paralogs, S_dev_S_tis_n_paralogs)
print('U vs S:', ks_U_S, pvalue_U_S)
ks_U_hS, pvalue_U_hS = stats.ks_2samp(S_dev_U_tis_n_paralogs, S_dev_hS_tis_n_paralogs)
print('U vs hS:', ks_U_hS, pvalue_U_hS)
ks_S_hS, pvalue_S_hS = stats.ks_2samp(S_dev_S_tis_n_paralogs, S_dev_hS_tis_n_paralogs)
print('S vs hS:', ks_S_hS, pvalue_S_hS)
print('significant difference between all the age distrib')

from scipy import stats
ks_U_S, pvalue_U_S = stats.ks_2samp(hS_dev_U_tis_n_paralogs, hS_dev_S_tis_n_paralogs)
print('U vs S:', ks_U_S, pvalue_U_S)
ks_U_hS, pvalue_U_hS = stats.ks_2samp(hS_dev_U_tis_n_paralogs, hS_dev_hS_tis_n_paralogs)
print('U vs hS:', ks_U_hS, pvalue_U_hS)
ks_S_hS, pvalue_S_hS = stats.ks_2samp(hS_dev_S_tis_n_paralogs, hS_dev_hS_tis_n_paralogs)
print('S vs hS:', ks_S_hS, pvalue_S_hS)
print('significant difference between all the age distrib')



#figures percentage paralog genes per gene subclass in dev

#mixed fig
time_unique=np.array(['0hpf', '0.75hpf', '2.25hpf', '3hpf', '4.3hpf', '5.25hpf', '6hpf', '8hpf', '10.3hpf', '16hpf', '19hpf', '24hpf', '30hpf', '36hpf', '2dpf', '3dpf', '4dpf', '5dpf'])

width = 0.25  # ancho de cada barra
x = np.arange(len(time_unique))
plt.figure(figsize=(11, 6), dpi=600)
plt.bar(x - width, S_dev_U_tis_act_stages_par * 100, width=width, label='S-Ut', color="#7E57C2")
plt.bar(x,         S_dev_S_tis_act_stages_par * 100, width=width, label='S-St', color= "#AB47BC")
plt.bar(x + width, S_dev_hS_tis_act_stages_par * 100, width=width, label='S-hSt', color="#EC407A")

plt.xticks(ticks=x, labels=time_unique, rotation=90, fontsize=25)
plt.yticks(fontsize=25)
plt.xlabel('Stages', fontsize=30)
plt.ylabel('% paralog genes\n(<TPM> > 1)', fontsize=28, fontweight='bold')
plt.legend(fontsize=24)
plt.title('S genes', fontsize=28)

plt.tight_layout()
plt.savefig(path + 'frac_paralog_genes_expression_per_stage_gene_subclass_S.png', dpi=600, bbox_inches='tight')
plt.show()

time_unique=np.array(['0hpf', '0.75hpf', '2.25hpf', '3hpf', '4.3hpf', '5.25hpf', '6hpf', '8hpf', '10.3hpf', '16hpf', '19hpf', '24hpf', '30hpf', '36hpf', '2dpf', '3dpf', '4dpf', '5dpf'])

width = 0.25  # ancho de cada barra
x = np.arange(len(time_unique))
plt.figure(figsize=(11, 6), dpi=600)
plt.bar(x - width, hS_dev_U_tis_act_stages_par * 100, width=width, label='hS-Ut', color="#F06292")
plt.bar(x,         hS_dev_S_tis_act_stages_par * 100, width=width, label='hS-St', color="#F48FB1")
plt.bar(x + width, hS_dev_hS_tis_act_stages_par * 100, width=width, label='hS-hSt', color="#F8BBD0")

plt.xticks(ticks=x, labels=time_unique, rotation=90, fontsize=25)
plt.yticks(fontsize=25)
plt.xlabel('Stages', fontsize=30)
plt.ylabel('% paralog genes\n(<TPM> > 1)', fontsize=28, fontweight='bold')
plt.legend(fontsize=24)
plt.title('hS genes', fontsize=28)

plt.tight_layout()
plt.savefig(path + 'frac_paralog_genes_expression_per_stage_gene_subclass_hS.png', dpi=600, bbox_inches='tight')
plt.show()






#====================================================================================
#====================================================================================
#Human ort
#10.1.) Ortology association
f=open(path+'human_orthos.txt', 'r')
txt = f.read()
human_orthos = txt.split('\n')
del txt, f
human_orthos=np.delete(human_orthos, len(human_orthos)-1)
human_orthos=np.array(human_orthos)

zfin_ort_ini = []
human_ort=[]

for elemento in human_orthos:
    partes = elemento.split("\t")

    zfin_ort_ini.append(partes[0])
    human_ort.append(partes[3])

zfin_ort_ini=np.array(zfin_ort_ini)
human_ort=np.array(human_ort)


unique_zfin_ort=np.unique(zfin_ort_ini)
human_ort_per_unique_zfin=[]
n_associated_human_ort=[]
#10.2.) Unique associations
for i in range(len(unique_zfin_ort)):
    ind_ort=np.where(zfin_ort_ini==unique_zfin_ort[i])[0]
    id_ort_hum=human_ort[ind_ort]
    human_ort_per_unique_zfin.append(np.unique(id_ort_hum))
    n_associated_human_ort.append(len(np.unique(id_ort_hum)))

n_associated_human_ort=np.array(n_associated_human_ort)



zfin_ort_common, ind_all, ind_zfin_to_human_ort=np.intersect1d(zfin_all, unique_zfin_ort, return_indices=True)

print('ort ENS:', len(zfin_ort_common))

final_genes_ort=ENS_all[ind_all]
n_associated_human_ort_final=n_associated_human_ort[ind_zfin_to_human_ort]


print('n human genes with ort:', np.sum(n_associated_human_ort_final))

final_human_ort_per_unique_zfin=[]
for i in range(len(ind_zfin_to_human_ort)):
    final_human_ort_per_unique_zfin.append(human_ort_per_unique_zfin[int(ind_zfin_to_human_ort[i])])

del zfin_ort_common, ind_all, ind_zfin_to_human_ort
del n_associated_human_ort

total_human_genes_ort=np.sum(n_associated_human_ort_final)



def orthologs(gene_analyse):
    common_genes_ort, ind_ort, ind_bulk=np.intersect1d(final_genes_ort, gene_analyse, return_indices=True)
    print(np.median(n_associated_human_ort_final[ind_ort]))
    return len(common_genes_ort)/len(gene_analyse)
    


#3.) We read the subclasses
label=['U_dev_U_tis', 'U_dev_S_tis', 'U_dev_hS_tis', 
       'S_dev_U_tis', 'S_dev_S_tis', 'S_dev_hS_tis',
       'hS_dev_U_tis', 'hS_dev_S_tis', 'hS_dev_hS_tis']
for k in (label):
    genes_analyze=np.loadtxt(path+"%s.txt" %k, dtype=str)
    genes_label=k
    globals()[f"{k}_frac_orthologs"]=orthologs(genes_analyze)
    

#one common figure
frac_genes_subclass=np.array([U_dev_U_tis_frac_orthologs, U_dev_S_tis_frac_orthologs, U_dev_hS_tis_frac_orthologs, 
                              S_dev_U_tis_frac_orthologs, S_dev_S_tis_frac_orthologs, S_dev_hS_tis_frac_orthologs, 
                              hS_dev_U_tis_frac_orthologs, hS_dev_S_tis_frac_orthologs, hS_dev_hS_tis_frac_orthologs])
groups=['U-Ut', 'U-St', 'U-hSt',
        'S-Ut', 'S-St', 'S-hSt', 
        'hS-Ut', 'hS-St', 'hS-hSt']
color=["#0D47A1", "#1976D2", "#5C6BC0", 
       "#7E57C2", "#AB47BC", "#EC407A", 
       "#F06292", "#F48FB1", "#F8BBD0"]
plt.figure(figsize=(4, 3), dpi=600)
plt.bar(groups, frac_genes_subclass, color=color)
plt.xticks(ticks=np.arange(len(groups)), labels=groups, fontsize=16, fontweight='bold', rotation=90)
plt.yticks(fontsize=14)
for i in range(len(groups)):
    plt.text(i-0.4, frac_genes_subclass[i], np.round(frac_genes_subclass[i], 2), fontsize=8, fontweight='bold', color='black')
plt.ylabel('% orthologs', fontsize=18, fontweight='bold')
plt.savefig(path+'frac_orthologs_per_gene_subclass.png', dpi=600, bbox_inches='tight')
plt.show()



