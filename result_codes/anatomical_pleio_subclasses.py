"""
Created on Mon Jan 27 14:22:31 2025

@author: logslab


Zebrafish ANATOMICAL ONTOLOGY 
GENE SUBCLASSES
"""

#construcción de la matriz fenotípica 

import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd
from collections import defaultdict


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

path='/media/alicia/TOSHIBA EXT/zebrafish/files(Copy)/'

# path='PATH_WHERE_YOU_KEEP_YOUR_DATASETS'


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


#3.3.) We find the layers and the phentypes
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


#GENE SUBCLASSES
f=open(path+'genes_bulk.txt', 'r')
txt = f.read()
genes_bulk = txt.split('\n')
del txt, f
genes_bulk=np.delete(genes_bulk, len(genes_bulk)-1)
genes_bulk=np.array(genes_bulk)

#bulk genes and matrix
f=open(path+'time_bulk.txt', 'r')
txt = f.read()
time_bulk = txt.split('\n')
del txt, f
time_bulk=np.delete(time_bulk, len(time_bulk)-1)
time_bulk=np.array(time_bulk)

mean_bulk_matrix= np.genfromtxt(path+'mean_bulk_matrix.txt', delimiter=',', dtype=None, encoding=None)



#3.) We read the subclasses
def pleio_and_lethal(genes_analyze, label):
    
    # #7.) We find the pleio associated with each gene subclass
    common_genes, ind_unique, ind_bulk =np.intersect1d(unique_genes, genes_analyze, return_indices=True)
    print(label, 'genes:', len(genes_analyze), ' - ', 'common:', len(common_genes))
    
    frac_lethal_genes=len(np.intersect1d(lethal_genes, common_genes))/len(common_genes)
    n_lethal=len(np.intersect1d(lethal_genes, common_genes))
    let_genes=np.intersect1d(lethal_genes, common_genes)
    
    
    #index of bulk genes
    g, i, ind_let=np.intersect1d(let_genes, genes_bulk, return_indices=True)
    lethal_matrix=mean_bulk_matrix[ind_let, :]
    
    act_stages_lethal=np.zeros(len(mean_bulk_matrix[0, :]))
    
    if len(ind_let)>0:
        for i in range(len(mean_bulk_matrix[0, :])):
            act_stages_lethal[i]=len(np.where(lethal_matrix[:, i]>1)[0])/len(ind_let)
       
    pleio_matrix=gene_phen_association_matrix[ind_unique, :]
    pleio_genes_analyze=np.sum(pleio_matrix, axis=1)  
    
    return n_lethal, frac_lethal_genes, pleio_genes_analyze, act_stages_lethal




label=['U_dev_U_tis', 'U_dev_S_tis', 'U_dev_hS_tis', 
       'S_dev_U_tis', 'S_dev_S_tis', 'S_dev_hS_tis',
       'hS_dev_U_tis', 'hS_dev_S_tis', 'hS_dev_hS_tis']

for k in (label):
    genes_analyze=np.loadtxt(path+"%s.txt" %k, dtype=str)
    genes_label=k
    globals()[f"{k}_n_lethal"], globals()[f"{k}_frac_lethal"], globals()[f"{k}_pleio"], globals()[f"{k}_act_stages_lethal"]=pleio_and_lethal(genes_analyze, genes_label)
    

#figure lethal
frac_genes_subclass=np.array([U_dev_U_tis_frac_lethal, U_dev_S_tis_frac_lethal, U_dev_hS_tis_frac_lethal, 
                              S_dev_U_tis_frac_lethal, S_dev_S_tis_frac_lethal, S_dev_hS_tis_frac_lethal, 
                              hS_dev_U_tis_frac_lethal, hS_dev_S_tis_frac_lethal, hS_dev_hS_tis_frac_lethal])
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
plt.ylabel('% lethal genes', fontsize=18, fontweight='bold')
plt.savefig(path+'frac_lethal_genes_per_gene_subclass.png', dpi=600, bbox_inches='tight')
plt.show()

frac_genes_subclass=np.array([U_dev_U_tis_n_lethal, U_dev_S_tis_n_lethal, U_dev_hS_tis_n_lethal, 
                              S_dev_U_tis_n_lethal, S_dev_S_tis_n_lethal, S_dev_hS_tis_n_lethal, 
                              hS_dev_U_tis_n_lethal, hS_dev_S_tis_n_lethal, hS_dev_hS_tis_n_lethal])
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
plt.ylabel('# lethal genes', fontsize=18, fontweight='bold')
plt.savefig(path+'n_lethal_genes_per_gene_subclass.png', dpi=600, bbox_inches='tight')
plt.show()



#figures percentage lethal genes per gene subclass in dev

#mixed fig
time_unique=np.array(['0hpf', '0.75hpf', '2.25hpf', '3hpf', '4.3hpf', '5.25hpf', '6hpf', '8hpf', '10.3hpf', '16hpf', '19hpf', '24hpf', '30hpf', '36hpf', '2dpf', '3dpf', '4dpf', '5dpf'])

width = 0.25  # ancho de cada barra
x = np.arange(len(time_unique))
plt.figure(figsize=(11, 6), dpi=600)
plt.bar(x - width, S_dev_U_tis_act_stages_lethal * 100, width=width, label='S-Ut', color="#7E57C2")
plt.bar(x,         S_dev_S_tis_act_stages_lethal * 100, width=width, label='S-St', color= "#AB47BC")
plt.bar(x + width, S_dev_hS_tis_act_stages_lethal * 100, width=width, label='S-hSt', color="#EC407A")

plt.xticks(ticks=x, labels=time_unique, rotation=90, fontsize=25)
plt.yticks(fontsize=25)
plt.xlabel('Stages', fontsize=30)
plt.ylabel('% lethal genes\n(<TPM> > 1)', fontsize=28, fontweight='bold')
plt.legend(fontsize=24)
plt.title('S genes', fontsize=28)

plt.tight_layout()
plt.savefig(path + 'frac_lethal_genes_expression_per_stage_gene_subclass_S.png', dpi=600, bbox_inches='tight')
plt.show()

time_unique=np.array(['0hpf', '0.75hpf', '2.25hpf', '3hpf', '4.3hpf', '5.25hpf', '6hpf', '8hpf', '10.3hpf', '16hpf', '19hpf', '24hpf', '30hpf', '36hpf', '2dpf', '3dpf', '4dpf', '5dpf'])

width = 0.25  # ancho de cada barra
x = np.arange(len(time_unique))
plt.figure(figsize=(11, 6), dpi=600)
plt.bar(x - width, hS_dev_U_tis_act_stages_lethal * 100, width=width, label='hS-Ut', color="#F06292")
plt.bar(x,         hS_dev_S_tis_act_stages_lethal * 100, width=width, label='hS-St', color="#F48FB1")
plt.bar(x + width, hS_dev_hS_tis_act_stages_lethal * 100, width=width, label='hS-hSt', color="#F8BBD0")

plt.xticks(ticks=x, labels=time_unique, rotation=90, fontsize=25)
plt.yticks(fontsize=25)
plt.xlabel('Stages', fontsize=30)
plt.ylabel('% lethal genes\n(<TPM> > 1)', fontsize=28, fontweight='bold')
plt.legend(fontsize=24)
plt.title('hS genes', fontsize=28)

plt.tight_layout()
plt.savefig(path + 'frac_lethal_genes_expression_per_stage_gene_subclass_hS.png', dpi=600, bbox_inches='tight')
plt.show()


#figure of anatomical pleio distrib
data=[U_dev_U_tis_pleio, U_dev_S_tis_pleio, U_dev_hS_tis_pleio, 
      S_dev_U_tis_pleio, S_dev_S_tis_pleio, S_dev_hS_tis_pleio, 
      hS_dev_U_tis_pleio, hS_dev_S_tis_pleio, hS_dev_hS_tis_pleio]
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
plt.ylabel("Anatomical pleiotropy", fontsize=13, fontweight='bold')
# plt.xlabel("Gene classes", fontsize=15)
plt.xticks(fontsize=18, fontweight='bold', rotation=90)
plt.yticks(fontsize=13)

plt.yscale('log')

plt.savefig(path+'pleiotropy_boxplot_subclasses.png', dpi=600, bbox_inches='tight')

plt.show()


#KS test
from scipy import stats
ks_U_S, pvalue_U_S = stats.ks_2samp(U_dev_U_tis_pleio, U_dev_S_tis_pleio)
print('U vs S:', ks_U_S, pvalue_U_S)
ks_U_hS, pvalue_U_hS = stats.ks_2samp(U_dev_U_tis_pleio, U_dev_hS_tis_pleio)
print('U vs hS:', ks_U_hS, pvalue_U_hS)
ks_S_hS, pvalue_S_hS = stats.ks_2samp(U_dev_S_tis_pleio, U_dev_hS_tis_pleio)
print('S vs hS:', ks_S_hS, pvalue_S_hS)
print('No significant difference between all the age distrib')
    
ks_U_S, pvalue_U_S = stats.ks_2samp(S_dev_U_tis_pleio, S_dev_S_tis_pleio)
print('U vs S:', ks_U_S, pvalue_U_S)
ks_U_hS, pvalue_U_hS = stats.ks_2samp(S_dev_U_tis_pleio, S_dev_hS_tis_pleio)
print('U vs hS:', ks_U_hS, pvalue_U_hS)
ks_S_hS, pvalue_S_hS = stats.ks_2samp(S_dev_S_tis_pleio, S_dev_hS_tis_pleio)
print('S vs hS:', ks_S_hS, pvalue_S_hS)
print('No significant difference between all the age distrib')

ks_U_S, pvalue_U_S = stats.ks_2samp(hS_dev_U_tis_pleio, hS_dev_S_tis_pleio)
print('U vs S:', ks_U_S, pvalue_U_S)
ks_U_hS, pvalue_U_hS = stats.ks_2samp(hS_dev_U_tis_pleio, hS_dev_hS_tis_pleio)
print('U vs hS:', ks_U_hS, pvalue_U_hS)
ks_S_hS, pvalue_S_hS = stats.ks_2samp(hS_dev_S_tis_pleio, hS_dev_hS_tis_pleio)
print('S vs hS:', ks_S_hS, pvalue_S_hS)
print('No significant difference between all the age distrib')






