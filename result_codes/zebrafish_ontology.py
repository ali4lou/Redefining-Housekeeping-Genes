# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 14:22:31 2025

@author: logslab


Zebrafish ANATOMICAL ONTOLOGY
"""

#construcción de la matriz fenotípica 

import obonet
import numpy as np
import os
import networkx as nx
import matplotlib.pyplot as plt 
import pandas as pd
from collections import defaultdict
from scipy.stats import mstats, kstest, ttest_ind, fisher_exact
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



def find_non_common_elements(arr1, arr2):
    set1 = set(arr1)
    set2 = set(arr2)
    non_common = list(set1.symmetric_difference(set2))  # Elements in either set1 or set2 but not both
    return non_common

def find_ancestors(reverse_graph, term):
    ancestors = set()
    stack = [term]
    while stack:
        node = stack.pop()
        for parent in reverse_graph[node]:
            if parent not in ancestors:
                ancestors.add(parent)
                stack.append(parent)
    return ancestors


path='PATH_WHERE_YOU_KEEP_YOUR_DATASETS'


#0.) Zfin phen descrip
f=open(path+'anatomy_phen_descriptions.txt', 'r')
txt = f.read()
zfin_to_desc = txt.split('\n')
del txt, f
zfin_to_desc=np.delete(zfin_to_desc, len(zfin_to_desc)-1)
zfin_to_desc=np.array(zfin_to_desc)

zfin_desc_all= []
desc_all=[]
ini_stage=[]
final_stage=[]

for elemento in zfin_to_desc:
    partes = elemento.split("\t")
    zfin_desc_all.append(partes[0])
    desc_all.append(partes[1])
    ini_stage.append(partes[2])
    final_stage.append(partes[3])

del zfin_to_desc


ini_stage=np.array(ini_stage, dtype=str)
final_stage=np.array(final_stage, dtype=str)

unique_ini_stage, n_times_ini=np.unique(ini_stage, return_counts=True)



#0.5.) Stage ZFIN

f=open(path+'stage_ontology.txt', 'r')
txt = f.read()
stage_ontology = txt.split('\n')
del txt, f
stage_ontology=np.delete(stage_ontology, len(stage_ontology)-1)
stage_ontology=np.array(stage_ontology)

sorted_unique_stages=[]
desc_stage=[]
for elemento in stage_ontology:
    partes = elemento.split("\t")
    sorted_unique_stages.append(partes[0])
    desc_stage.append(partes[2])


desc_stage=np.delete(desc_stage, 0)
sorted_unique_stages=np.delete(sorted_unique_stages, 0)


n_times_desc_stage=np.zeros(len(desc_stage))
for i in range(len(desc_stage)):
    ind=np.where(unique_ini_stage==sorted_unique_stages[i])[0]
    print(ind)
    if len(ind>0):
        n_times_desc_stage[i]=n_times_ini[int(ind)]
        
n_zfin=0
total_zfin_per_stage=np.zeros(len(sorted_unique_stages))
final=int(len(sorted_unique_stages))
for i in range(len(zfin_desc_all)):
    ind_ini=np.where(sorted_unique_stages==ini_stage[i])[0]
    ind_final=np.where(sorted_unique_stages==final_stage[i])[0]
    if (len(ind_ini)>0) & (len(ind_final)>0):
        n_zfin=n_zfin+1
        for count in range(int(ind_ini), int(ind_final)):
            # print(count, i,int(ind_ini), int(ind_final))
            total_zfin_per_stage[count]=total_zfin_per_stage[count]+1
    else:
        if len(ind_final)>0:
            n_zfin=n_zfin+1
            total_zfin_per_stage[int(ind_final)]=total_zfin_per_stage[int(ind_final)]+1            
        if len(ind_ini)>0:
            n_zfin=n_zfin+1
            n_steps=final-int(ind_ini)
            for count in range(int(ind_ini), n_steps):
                total_zfin_per_stage[count]=total_zfin_per_stage[count]+1
                # print(count, int(ind_ini), n_steps)
        
    
#figure
fig, ax1 = plt.subplots(figsize=(8, 3), dpi=600)

ax1.plot(desc_stage, total_zfin_per_stage, 'mediumvioletred')  # Línea roja
ax1.set_ylabel("# total ZFAs", color="mediumvioletred", fontsize=18, fontweight='bold')
ax1.tick_params(axis="y", labelcolor="mediumvioletred", labelsize=15)

ax2 = ax1.twinx()
ax2.plot(desc_stage, n_times_desc_stage, 'slateblue')  # Línea azul
ax2.set_ylabel("# novel ZFAs", color="slateblue", fontsize=18, fontweight='bold')
ax2.tick_params(axis="y", labelcolor="slateblue", labelsize=15)

ax1.set_xticks(np.arange(len(desc_stage)))  # Posiciones en el eje X
ax1.set_xticklabels(desc_stage, rotation=90, fontsize=10)
plt.savefig(path+'ontology_along_stages.png', dpi=600, bbox_inches='tight')
plt.show()



#1.) Zfin gene related to ENS gene
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



#2.) gene-phen associations
with open(path+'phen_gene_data.txt', 'r') as file:
    data = file.read()
rows = data.split('\n')
rows = [row.strip() for row in rows if row.strip()]
del data, file

zfin_associated = []
phen_associated_id=[]
phen_associated_desc=[]
phen_thing=[]

for elemento in rows:
    partes = elemento.split("\t")
    
    for i in range(len(partes)):
        if partes[i].startswith("ZFA"):
            zfin_associated.append(partes[2])
            phen_associated_id.append(partes[i])
            phen_associated_desc.append(partes[int(i+1)])
            phen_thing.append(partes[10])
    
    

phen_associated_desc=np.array(phen_associated_desc)
phen_associated_id=np.array(phen_associated_id)
zfin_associated=np.array(zfin_associated)
phen_thing=np.array(phen_thing)

del rows, partes



#2.0.) We find the ENS genes associated with 'lethal'
lethal_ind=np.where(phen_thing=='lethal (sensu genetics)')[0]
print('Associations with lethal:', len(lethal_ind))
phen_id_associated_with_lethal=phen_associated_desc[lethal_ind]
phen_id_associated_with_lethal_unique, times_letal_phen=np.unique(phen_id_associated_with_lethal, return_counts=True)
lethal_zfin_genes=zfin_associated[lethal_ind]
common_g_lethal_zfin, ind_all_lethal, i=np.intersect1d(zfin_all, lethal_zfin_genes, return_indices=True)
lethal_genes=ENS_all[ind_all_lethal]


print('unique zfin genes:', len(np.unique(lethal_zfin_genes)))

#2.1.) We search the ENS genes associted
ENS_associated=[]
phen_associated_id_final=[]
phen_associated_desc_final=[]
phen_thing_final=[]
count_rep_genes=0
for i in range(len(zfin_associated)):
    ind_gene=np.where(zfin_all==zfin_associated[i])[0]
    if len(ind_gene)==1:
        phen_associated_id_final.append(phen_associated_id[i])
        phen_associated_desc_final.append(phen_associated_desc[i])
        phen_thing_final.append(phen_thing[i])
        ENS_associated.append(ENS_all[int(ind_gene)])
    if len(ind_gene)>1:
        count_rep_genes=count_rep_genes+1
        for j in range(len(ind_gene)):
            phen_associated_id_final.append(phen_associated_id[i])
            phen_associated_desc_final.append(phen_associated_desc[i])
            phen_thing_final.append(phen_thing[i])
            ENS_associated.append(ENS_all[int(ind_gene[j])])
    
ENS_associated=np.array(ENS_associated)
phen_associated_id_final=np.array(phen_associated_id_final)
phen_associated_desc_final=np.array(phen_associated_desc_final)
phen_thing_final=np.array(phen_thing_final)

del zfin_associated, phen_associated_desc, phen_associated_id, phen_thing


#3.) We search the ontology 
zfa_children=[]
zfa_parent=[]
with open(path+'anatomy_relationship.txt', 'r') as file:
    input_data = file.read()
rows = input_data.split('\n')
rows = [row.strip() for row in rows if row.strip()]
del file

for elemento in rows:
    partes = elemento.split("\t")
    if partes[2].startswith("is"):
        zfa_children.append(partes[1])
        zfa_parent.append(partes[0])
    if partes[2].startswith("part"):
        zfa_children.append(partes[1])
        zfa_parent.append(partes[0])

zfa_parent=np.array(zfa_parent)
zfa_children=np.array(zfa_children)


#exploratory analysis of the number of terms

unique_child=np.unique(zfa_children)
unique_parent=np.unique(zfa_parent)

child_inters_parent=np.intersect1d(unique_child, unique_parent)

specific_child=find_non_common_elements(unique_child, unique_parent)
real_specific_child=np.intersect1d(unique_child, specific_child)
top_ontology=np.intersect1d(unique_parent, specific_child)

n_phen_ont=len(real_specific_child)+len(unique_parent)
print('total zfa ontology', n_phen_ont)

print('We find %d specific pehnotypes' %len(real_specific_child))

zfa_unique=np.unique(phen_associated_id_final)
print('Total unique phenotypes associated:', len(zfa_unique))
print('Speficific phenotypes from that total:', len(np.intersect1d(real_specific_child, phen_associated_id_final)))
print('Parent phenotypes from that total:', len(np.intersect1d(unique_parent, phen_associated_id_final)))


# 3.1.) Build the graph and its reverse
graph = defaultdict(list)
reverse_graph = defaultdict(list)
for parent, child in zip(zfa_parent, zfa_children):
    graph[parent].append(child)
    reverse_graph[child].append(parent)

#3.2.) Find ancestors for a given term
results = {}
ancestors_list=[]
for term in zfa_unique:
    results[term] = find_ancestors(reverse_graph, term)
    ancestors = find_ancestors(reverse_graph, term)
    ancestors_list.append([term, list(ancestors)])

#3.2.1.) We find all the needed terms
all_redundant_terms=[]
for i in range(len(zfa_unique)):
    all_redundant_terms.append(ancestors_list[i][0])
    for j in range(len(ancestors_list[i][1])):
        all_redundant_terms.append(ancestors_list[i][1][j])

all_redundant_terms_unique=np.unique(np.array(all_redundant_terms))


#3.3.) We find the layers and the phenotypes
layers=[]
layers.append(top_ontology)
count_phen=1
for i in range(17):
    inner_terms=[]
    for j in range(len(layers[i])):
        ind_zfa=np.where(zfa_parent==layers[i][j])[0]
        for k in range(len(ind_zfa)):
            inner_terms.append(zfa_children[int(ind_zfa[k])])
            count_phen=count_phen+1
    layers.append(inner_terms)
    

        

#3.3.1.) We sort the associated phenotypes by layer
unique_genes=np.unique(ENS_associated)

zfa_layer=np.zeros(len(all_redundant_terms_unique))
for i in range(len(all_redundant_terms_unique)):
    for j in range(len(layers)):
        for k in range(len(layers[j])):
            if layers[j][k]==all_redundant_terms_unique[i]:
                if (zfa_layer[i]!=0): 
                    if j<zfa_layer[i]:
                        zfa_layer[i]=j
                else:
                    zfa_layer[i]=j


zfin_desc_all.append(top_ontology[0])
desc_all.append('zebrafish anatomical entity')

zfin_desc_all=np.array(zfin_desc_all)
desc_all=np.array(desc_all)

zfa_unique_desc=[]
for i in range(len(all_redundant_terms_unique)):
    ind=np.where(zfin_desc_all==all_redundant_terms_unique[i])[0]
    zfa_unique_desc.append(desc_all[int(ind[0])])
    
df_zfa=pd.DataFrame()
df_zfa['zfa']=all_redundant_terms_unique
df_zfa['desc']=zfa_unique_desc
df_zfa['layer']=zfa_layer

df_zfa=df_zfa.sort_values(by=['layer'])

zfa_unique_matrix=np.array(list(df_zfa['zfa']))
zfa_unique_desc_matrix=np.array(list(df_zfa['desc']))


#4.) We built the gene-phenotype association matrix
gene_phen_association_matrix=np.zeros((len(unique_genes), len(zfa_unique_matrix)))
for i in range(len(unique_genes)):
    ind_association=np.where(ENS_associated==unique_genes[i])
    searched_phen=phen_associated_id_final[ind_association]
    for j in range(len(searched_phen)):
        ind_phen=np.where(zfa_unique==searched_phen[j])[0]
        ind_matrix=np.where(zfa_unique_matrix==zfa_unique[int(ind_phen)])[0]
        gene_phen_association_matrix[i][int(ind_matrix)]=1
        for k in range(len(ancestors_list[int(ind_phen)][1])):
            ind_matrix=np.where(zfa_unique_matrix==ancestors_list[int(ind_phen)][1][k])[0]
            gene_phen_association_matrix[i][int(ind_matrix)]=1

np.savetxt(path+'gene_phen_association_matrix.txt', gene_phen_association_matrix)
np.savetxt(path+'genes_associated_phen.txt', unique_genes, fmt='%s', delimiter=',')
np.savetxt(path+'phen_anatomy.txt', zfa_unique_desc_matrix, fmt='%s', delimiter=',')



#5.) Pleiotropy
pleio=np.sum(gene_phen_association_matrix, axis=1)
plt.hist(pleio, bins=100)

        
        
#6.) We read the txt with the UE, the hS and S genes from bulk tissue

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


# #7.) We find the pleio associated with each gene
common_U_genes, ind_U_unique, ind_U_bulk =np.intersect1d(unique_genes, U_id_bulk, return_indices=True)
print('U genes:', len(U_id_bulk), ' - ', 'common U:', len(common_U_genes))
        
common_S_genes, ind_S_unique, ind_S_bulk=np.intersect1d(unique_genes, S_id_bulk, return_indices=True)
print('S genes:', len(S_id_bulk), ' - ', 'common S:', len(common_S_genes))
         
common_hS_genes, ind_hS_unique, ind_hS_bulk=np.intersect1d(unique_genes, hS_id_bulk, return_indices=True)
print('hS genes:', len(hS_id_bulk), ' - ', 'common hS:', len(common_hS_genes))
                
    
#7.0.) Lethal genes
check_eye=np.intersect1d(genes_bulk, lethal_genes)

U_genes_lethal=np.intersect1d(lethal_genes, common_U_genes)
S_genes_lethal=np.intersect1d(lethal_genes, common_S_genes)
hS_genes_lethal=np.intersect1d(lethal_genes, common_hS_genes)

print(len(U_genes_lethal)/len(lethal_genes), len(S_genes_lethal)/len(lethal_genes), len(hS_genes_lethal)/len(lethal_genes))
print(len(check_eye), len(U_genes_lethal)+len(S_genes_lethal)+len(hS_genes_lethal))

test_hipergeometrico(len(genes_bulk), len(check_eye), len(U_id_bulk), len(U_genes_lethal))
test_hipergeometrico(len(genes_bulk), len(check_eye), len(S_id_bulk), len(S_genes_lethal))
test_hipergeometrico(len(genes_bulk), len(check_eye), len(hS_id_bulk), len(hS_genes_lethal))



#figure
frac_genes=np.array([len(U_genes_lethal)/len(check_eye)*100, len(S_genes_lethal)/len(check_eye)*100, len(hS_genes_lethal)/len(check_eye)*100])
groups=['U', 'S', 'hS']
color=['royalblue', 'mediumturquoise', 'tomato']

plt.figure(figsize=(4, 3), dpi=600)
plt.bar(groups, frac_genes, color=color)
plt.xticks(ticks=np.arange(len(groups)), labels=groups, fontsize=16, fontweight='bold')
plt.yticks(fontsize=14)
for i in range(len(groups)):
    plt.text(i-0.25, frac_genes[i]-6, np.round(frac_genes[i], 2), fontsize=12, fontweight='bold', color='white')
if "U" in groups:
    idx_U = groups.index("U")  # posición de la barra 'S'
    plt.text(idx_U, frac_genes[idx_U]-3, "*",
             ha='center', va='bottom', fontsize=18, color='red')
plt.ylabel('% lethal genes', fontsize=18, fontweight='bold')
plt.savefig(path+'frac_lethal_genes_ZFA.png', dpi=600, bbox_inches='tight')
plt.show()

#figure of % expressed lethal genes in each stage per gene class

#index of bulk genes
g, i, ind_U_let=np.intersect1d(U_genes_lethal, genes_bulk, return_indices=True)
g, i, ind_S_let=np.intersect1d(S_genes_lethal, genes_bulk, return_indices=True)
g, i, ind_hS_let=np.intersect1d(hS_genes_lethal, genes_bulk, return_indices=True)
   
lethal_U_matrix=mean_bulk_matrix[ind_U_let, :]
lethal_S_matrix=mean_bulk_matrix[ind_S_let, :]
lethal_hS_matrix=mean_bulk_matrix[ind_hS_let, :]

act_stages_lethal_U=np.zeros(len(mean_bulk_matrix[0, :]))
act_stages_lethal_S=np.zeros(len(mean_bulk_matrix[0, :]))
act_stages_lethal_hS=np.zeros(len(mean_bulk_matrix[0, :]))

for i in range(len(mean_bulk_matrix[0, :])):
    act_stages_lethal_U[i]=len(np.where(lethal_U_matrix[:, i]>1)[0])/len(ind_U_let)
    # print(len(np.where(lethal_U_matrix[:, i]>1)[0]), len(ind_U_let))
    act_stages_lethal_S[i]=len(np.where(lethal_S_matrix[:, i]>1)[0])/len(ind_S_let)
    act_stages_lethal_hS[i]=len(np.where(lethal_hS_matrix[:, i]>1)[0])/len(ind_hS_let)
    print(len(np.where(lethal_hS_matrix[:, i]>1)[0]), len(ind_hS_let))


#mixed fig
time_unique=np.array(['0hpf', '0.75hpf', '2.25hpf', '3hpf', '4.3hpf', '5.25hpf', '6hpf', '8hpf', '10.3hpf', '16hpf', '19hpf', '24hpf', '30hpf', '36hpf', '2dpf', '3dpf', '4dpf', '5dpf'])

width = 0.25  # ancho de cada barra
x = np.arange(len(time_unique))
plt.figure(figsize=(11, 6), dpi=600)
plt.bar(x - width, act_stages_lethal_U * 100, width=width, label='U', color='royalblue')
plt.bar(x,         act_stages_lethal_S * 100, width=width, label='S', color='mediumturquoise')
plt.bar(x + width, act_stages_lethal_hS * 100, width=width, label='hS', color='tomato')

plt.xticks(ticks=x, labels=time_unique, rotation=90, fontsize=25)
plt.yticks(fontsize=25)
plt.xlabel('Stages', fontsize=30)
plt.ylabel('% lethal genes\n(<TPM> > 1)', fontsize=28, fontweight='bold')
plt.legend(fontsize=20)

plt.tight_layout()
plt.savefig(path + 'frac_lethal_genes_expression_per_stage_gene_class.png', dpi=600, bbox_inches='tight')
plt.show()

    

#7.1.) pleiotropy
gene_phen_matrix_U=gene_phen_association_matrix[ind_U_unique, :]
gene_phen_matrix_S=gene_phen_association_matrix[ind_S_unique, :]
gene_phen_matrix_hS=gene_phen_association_matrix[ind_hS_unique, :]

pleio_U=np.sum(gene_phen_matrix_U, axis=1)  
pleio_S=np.sum(gene_phen_matrix_S, axis=1)  
pleio_hS=np.sum(gene_phen_matrix_hS, axis=1)  
   
#figure
plt.figure(figsize=(4, 3), dpi=600)
plt.hist(pleio_U, color='royalblue', label='U', log=True)
plt.hist(pleio_S, color='mediumturquoise', label='S', log=True)
plt.hist(pleio_hS, color='tomato', label='hS', log=True)
plt.xlabel('Anatomical pleiotropy', fontsize=18, fontweight='bold')
plt.ylabel('# genes', fontsize=18)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.legend(fontsize=13)
plt.grid(False)
plt.tight_layout()
plt.savefig(path+'pleiotropy.png', dpi=600, bbox_inches='tight')
plt.show()


#figure boxplot
data=[pleio_U, pleio_S, pleio_hS]
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

plt.ylabel("Anatomical pleiotropy", fontsize=13, fontweight='bold')
# plt.xlabel("Gene classes", fontsize=15)
plt.xticks(fontsize=18, fontweight='bold')
plt.yticks(fontsize=13)
plt.yscale('log')

comparisons = [(0,1), (0,2), (1,2)]  
labels = [' ', '****', '****']   
y_max = max([max(d) for d in data])
h = 0.4 * y_max  
line_heights = [y_max + 1.0*y_max, y_max + 2.2*y_max, y_max + 3.4*y_max]
for i, ((x1, x2), label) in enumerate(zip(comparisons, labels)):
    y = line_heights[i]
    ax.plot([x1+1, x1+1, x2+1, x2+1], [y, y+h, y+h, y], lw=0.5, c='black')
    ax.text((x1 + x2 + 2) / 2, y + 0.1 * y_max, label, ha='center', fontsize=10, c='red')

plt.savefig(path+'pleiotropy_boxplot_log.png', dpi=600, bbox_inches='tight')

plt.show()


#7.2.) ks test between pleio distrib
from scipy import stats
ks_U_S, pvalue_U_S = stats.ks_2samp(pleio_U, pleio_S)
print('U vs S:', ks_U_S, pvalue_U_S)
ks_U_hS, pvalue_U_hS = stats.ks_2samp(pleio_U, pleio_hS)
print('U vs hS:', ks_U_hS, pvalue_U_hS)
ks_S_hS, pvalue_S_hS = stats.ks_2samp(pleio_S, pleio_hS)
print('S vs hS:', ks_S_hS, pvalue_S_hS)
print('significant difference between hS and S, U pleiotropy distrib')
    
print('median:', np.median(pleio_U), np.median(pleio_S), np.median(pleio_hS))

#7.3.) Pleiotropy and development in bulk
g, ind_bulk_hS, i=np.intersect1d(genes_bulk, common_hS_genes, return_indices=True)
g, ind_bulk_S, i=np.intersect1d(genes_bulk, common_S_genes, return_indices=True)
g, ind_bulk_U, i=np.intersect1d(genes_bulk, common_U_genes, return_indices=True)

bulk_U=mean_bulk_matrix[ind_bulk_U, :]
bulk_S=mean_bulk_matrix[ind_bulk_S, :]
bulk_hS=mean_bulk_matrix[ind_bulk_hS, :]

#7.3.1.) For each gene we find the stage of maximal exp 
stage_of_max_exp_U=np.zeros(len(common_U_genes))
for i in range(len(common_U_genes)):
    ind_max=np.argmax(bulk_U[i, :])
    stage_of_max_exp_U[i]=ind_max
    
stage_of_max_exp_S=np.zeros(len(common_S_genes))
for i in range(len(common_S_genes)):
    ind_max=np.argmax(bulk_S[i, :])
    stage_of_max_exp_S[i]=ind_max
    
stage_of_max_exp_hS=np.zeros(len(common_hS_genes))
for i in range(len(common_hS_genes)):
    ind_max=np.argmax(bulk_hS[i, :])
    stage_of_max_exp_hS[i]=ind_max

    
#7.3.2.) We perform a sliding window of its anatomical pleiotropy
sw=400
#U genes
stage_of_max_exp_sorted_U=np.sort(stage_of_max_exp_U)
index_stage_sort_U=np.argsort(stage_of_max_exp_U)
pleio_per_gen_sorted_U=np.zeros(len(stage_of_max_exp_U))
for i in range(len(index_stage_sort_U)):
    pleio_per_gen_sorted_U[i]=pleio_U[int(index_stage_sort_U[i])]

serie1U = pd.Series(stage_of_max_exp_sorted_U)
serie2U = pd.Series(pleio_per_gen_sorted_U)

sw_stage_U = serie1U.rolling(window=sw, center=False).mean()
sw_pleio_U = serie2U.rolling(window=sw, center=False).mean()

#S genes
stage_of_max_exp_sorted_S=np.sort(stage_of_max_exp_S)
index_stage_sort_S=np.argsort(stage_of_max_exp_S)
pleio_per_gen_sorted_S=np.zeros(len(stage_of_max_exp_S))
for i in range(len(index_stage_sort_S)):
    pleio_per_gen_sorted_S[i]=pleio_S[int(index_stage_sort_S[i])]

serie1S = pd.Series(stage_of_max_exp_sorted_S)
serie2S = pd.Series(pleio_per_gen_sorted_S)

sw_stage_S = serie1S.rolling(window=sw, center=False).mean()
sw_pleio_S = serie2S.rolling(window=sw, center=False).mean()

#hS genes
stage_of_max_exp_sorted_hS=np.sort(stage_of_max_exp_hS)
index_stage_sort_hS=np.argsort(stage_of_max_exp_hS)
pleio_per_gen_sorted_hS=np.zeros(len(stage_of_max_exp_hS))
for i in range(len(index_stage_sort_hS)):
    pleio_per_gen_sorted_hS[i]=pleio_hS[int(index_stage_sort_hS[i])]

serie1hS = pd.Series(stage_of_max_exp_sorted_hS)
serie2hS = pd.Series(pleio_per_gen_sorted_hS)

sw_stage_hS = serie1hS.rolling(window=sw, center=False).mean()
sw_pleio_hS = serie2hS.rolling(window=sw, center=False).mean()



# # figure
plt.figure(figsize=(5, 3), dpi=600)
plt.scatter(sw_stage_U, sw_pleio_U, s=1, color='royalblue', label='U')
plt.scatter(sw_stage_S, sw_pleio_S, s=1, color='mediumturquoise', label='S')
plt.scatter(sw_stage_hS, sw_pleio_hS, s=1, color='tomato', label='hS')
plt.xlabel('Stage max. expression', fontsize=13, fontweight='bold')
plt.ylabel('Anatomical pleiotropy', fontsize=13, fontweight='bold')
plt.xticks(np.arange(len(time_bulk)), time_bulk, fontsize=13, rotation=90)
plt.yticks(fontsize=13)
plt.legend(markerscale=5)
plt.savefig(path+'pleiotropy_vs_max_stage_exp_sw_%d.png' %sw, dpi=600, bbox_inches='tight')
plt.show()


#8.) Enrichment of phenotypes (Fisher)
# def enrichement_phen(submatrix_gene_phen, gene_array, phen, genes, phen_matrix):
phen=zfa_unique_desc_matrix

    
phen_n_times_matching=np.zeros(len(phen))
#Antes de calcular el enriquecimiento necesitamos saber cuántas veces sale el fenotipo
#en la matriz general de genes coincidentes
for i in range(len(phen)):
    phen_n_times_matching[i]=np.sum(gene_phen_association_matrix[:, i])
N_genes=len(unique_genes)    

    
#8.1.) ENRICHMENT U genes
odd_ratio_enriched=np.zeros(len(phen))
p_value_enriched=np.zeros(len(phen))
n_genes_subset=len(common_U_genes)
phen_enriched_fisher_U=[]
p_val_pleio=[]


n_genes_phen=[]
n_genes_phen_subset=[]
for fen in range(len(phen)):
    phen_n_times_subset=np.sum(gene_phen_matrix_U[:, fen])
    tabla=[[phen_n_times_subset, n_genes_subset-phen_n_times_subset],[phen_n_times_matching[fen], N_genes-phen_n_times_matching[fen]]]
    odd_ratio_enriched[fen], p_value_enriched[fen] = fisher_exact(tabla, alternative="greater") 
    if p_value_enriched[fen]<0.005:
        phen_enriched_fisher_U.append(phen[fen])
        p_val_pleio.append(p_value_enriched[fen])
        n_genes_phen.append(phen_n_times_matching[fen])
        n_genes_phen_subset.append(phen_n_times_subset)


p_val_pleio=np.array(p_val_pleio)
n_genes_phen=np.array(n_genes_phen)
n_genes_phen_subset=np.array(n_genes_phen_subset)
phen_enriched_fisher_U=np.array(phen_enriched_fisher_U)
ind_sorted=np.argsort(p_val_pleio)
p_val_pleio_sorted=p_val_pleio[ind_sorted]
phen_enriched_fisher_pleio_sorted=phen_enriched_fisher_U[ind_sorted]
n_genes_phen_sorted=n_genes_phen[ind_sorted]
n_genes_phen_subset_sorted=n_genes_phen_subset[ind_sorted]


df_pleio=pd.DataFrame()
df_pleio['phen']=phen_enriched_fisher_pleio_sorted
df_pleio['p-value']=p_val_pleio_sorted
df_pleio['n_genes_phen']=n_genes_phen_sorted
df_pleio['n_genes_subset_phen']=n_genes_phen_subset_sorted
df_pleio.to_csv(path+'enriched_phen_U_genes.csv', sep='\t')
    


#8.2.) ENRICHMENT S genes
odd_ratio_enriched=np.zeros(len(phen))
p_value_enriched=np.zeros(len(phen))
n_genes_subset=len(common_S_genes)
phen_enriched_fisher_S=[]
p_val_pleio=[]


n_genes_phen=[]
n_genes_phen_subset=[]
for fen in range(len(phen)):
    phen_n_times_subset=np.sum(gene_phen_matrix_S[:, fen])
    tabla=[[phen_n_times_subset, n_genes_subset-phen_n_times_subset],[phen_n_times_matching[fen], N_genes-phen_n_times_matching[fen]]]
    odd_ratio_enriched[fen], p_value_enriched[fen] = fisher_exact(tabla, alternative="greater") 
    if p_value_enriched[fen]<0.005:
        phen_enriched_fisher_S.append(phen[fen])
        p_val_pleio.append(p_value_enriched[fen])
        n_genes_phen.append(phen_n_times_matching[fen])
        n_genes_phen_subset.append(phen_n_times_subset)


p_val_pleio=np.array(p_val_pleio)
n_genes_phen=np.array(n_genes_phen)
n_genes_phen_subset=np.array(n_genes_phen_subset)
phen_enriched_fisher_S=np.array(phen_enriched_fisher_S)
ind_sorted=np.argsort(p_val_pleio)
p_val_pleio_sorted=p_val_pleio[ind_sorted]
phen_enriched_fisher_pleio_sorted=phen_enriched_fisher_S[ind_sorted]
n_genes_phen_sorted=n_genes_phen[ind_sorted]
n_genes_phen_subset_sorted=n_genes_phen_subset[ind_sorted]


df_pleio=pd.DataFrame()
df_pleio['phen']=phen_enriched_fisher_pleio_sorted
df_pleio['p-value']=p_val_pleio_sorted
df_pleio['n_genes_phen']=n_genes_phen_sorted
df_pleio['n_genes_subset_phen']=n_genes_phen_subset_sorted
df_pleio.to_csv(path+'enriched_phen_S_genes.csv', sep='\t')
  

#8.3.) ENRICHMENT hS genes
odd_ratio_enriched=np.zeros(len(phen))
p_value_enriched=np.zeros(len(phen))
n_genes_subset=len(common_hS_genes)
phen_enriched_fisher_hS=[]
p_val_pleio=[]


n_genes_phen=[]
n_genes_phen_subset=[]
for fen in range(len(phen)):
    phen_n_times_subset=np.sum(gene_phen_matrix_hS[:, fen])
    tabla=[[phen_n_times_subset, n_genes_subset-phen_n_times_subset],[phen_n_times_matching[fen], N_genes-phen_n_times_matching[fen]]]
    odd_ratio_enriched[fen], p_value_enriched[fen] = fisher_exact(tabla, alternative="greater") 
    if p_value_enriched[fen]<0.005:
        phen_enriched_fisher_hS.append(phen[fen])
        p_val_pleio.append(p_value_enriched[fen])
        n_genes_phen.append(phen_n_times_matching[fen])
        n_genes_phen_subset.append(phen_n_times_subset)


p_val_pleio=np.array(p_val_pleio)
n_genes_phen=np.array(n_genes_phen)
n_genes_phen_subset=np.array(n_genes_phen_subset)
phen_enriched_fisher_hS=np.array(phen_enriched_fisher_hS)
ind_sorted=np.argsort(p_val_pleio)
p_val_pleio_sorted=p_val_pleio[ind_sorted]
phen_enriched_fisher_pleio_sorted=phen_enriched_fisher_hS[ind_sorted]
n_genes_phen_sorted=n_genes_phen[ind_sorted]
n_genes_phen_subset_sorted=n_genes_phen_subset[ind_sorted]


df_pleio=pd.DataFrame()
df_pleio['phen']=phen_enriched_fisher_pleio_sorted
df_pleio['p-value']=p_val_pleio_sorted
df_pleio['n_genes_phen']=n_genes_phen_sorted
df_pleio['n_genes_subset_phen']=n_genes_phen_subset_sorted
df_pleio.to_csv(path+'enriched_phen_hS_genes.csv', sep='\t')
  





