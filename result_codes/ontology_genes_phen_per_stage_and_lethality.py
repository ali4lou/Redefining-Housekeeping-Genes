# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 14:22:31 2025

@author: logslab


Zebrafish ANATOMICAL ONTOLOGY
"""


import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.patches as mpatches


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


#1.) Stage and phenotype
#We read commmon phenotypes
f=open(path+'phen_stage_analize.txt', 'r')
txt = f.read()
phen = txt.split('\n')
del txt, f
phen=np.delete(phen, len(phen)-1)
phen=np.array(phen)

final_stage=[]
ini_stage=[]
phen_stage=[]
for i in range(len(phen)):
    part = phen[i].split('\t')
    phen_stage.append(part[1])
    ini_stage.append(part[2])
    final_stage.append(part[3])
    
final_stage=np.array(final_stage)
phen_stage=np.array(phen_stage)
ini_stage=np.array(ini_stage)


gene_phen_association_matrix=np.loadtxt(path+'gene_phen_association_matrix.txt')

f=open(path+'genes_associated_phen.txt', 'r')
txt = f.read()
genes_associated_phen = txt.split('\n')
del txt, f
genes_associated_phen=np.delete(genes_associated_phen, len(genes_associated_phen)-1)
genes_associated_phen=np.array(genes_associated_phen)

f=open(path+'phen_anatomy.txt', 'r')
txt = f.read()
phen_anatomy = txt.split('\n')
del txt, f
phen_anatomy=np.delete(phen_anatomy, len(phen_anatomy)-1)
phen_anatomy=np.array(phen_anatomy)


#2.) gene-stage matrix
gene_stage_matrix=np.zeros((len(genes_associated_phen), len(desc_stage)))
for i in range(len(genes_associated_phen)):
    ind_phen_ass=np.where(gene_phen_association_matrix[i, :]>0)[0]
    phen_ass=phen_anatomy[ind_phen_ass]
    for j in phen_ass:
        ind_stage=np.where(phen_stage==j)[0]
        if len(ind_stage)>0:
            ind_ini_gene_stage=np.where(sorted_unique_stages==ini_stage[ind_stage])[0]
            ind_final_gene_stage=np.where(sorted_unique_stages==final_stage[ind_stage])[0]
            if (len(ind_ini_gene_stage)>0) & (len(ind_final_gene_stage)>0):
                for k in range(int(ind_ini_gene_stage), int(ind_final_gene_stage+1)):
                    gene_stage_matrix[i][k]=gene_stage_matrix[i][k]+1
    gene_stage_matrix[i, :]=gene_stage_matrix[i, :]/np.sum(gene_phen_association_matrix[i, :])
    print(i)



        
#3.) We read the txt with the UE, the hS and S genes from bulk tissue

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


#4.) Genes in the phen
common_U_genes, ind_U_unique, ind_U_bulk =np.intersect1d(genes_associated_phen, U_id_bulk, return_indices=True)
print('U genes:', len(U_id_bulk), ' - ', 'common U:', len(common_U_genes))
        
common_S_genes, ind_S_unique, ind_S_bulk=np.intersect1d(genes_associated_phen, S_id_bulk, return_indices=True)
print('S genes:', len(S_id_bulk), ' - ', 'common S:', len(common_S_genes))
         
common_hS_genes, ind_hS_unique, ind_hS_bulk=np.intersect1d(genes_associated_phen, hS_id_bulk, return_indices=True)
print('hS genes:', len(hS_id_bulk), ' - ', 'common hS:', len(common_hS_genes))
                

#7.1.) Gene - stage phenotype
gene_phen_matrix_U=gene_stage_matrix[ind_U_unique, :]
gene_phen_matrix_S=gene_stage_matrix[ind_S_unique, :]
gene_phen_matrix_hS=gene_stage_matrix[ind_hS_unique, :]


#figure
plt.figure(figsize=(18,6), dpi=600)

offset = 0.25
positions = np.arange(len(desc_stage))

# Propiedades de las cajas (sin borde)
box_props_U = dict(facecolor="royalblue", edgecolor="none")
box_props_S = dict(facecolor="mediumturquoise", edgecolor="none")
box_props_hS = dict(facecolor="tomato", edgecolor="none")

# Solo mantener la mediana visible
median_props = dict(color="black", linewidth=0.7)

# Dibujamos los boxplots
plt.boxplot([gene_phen_matrix_U[:, i] for i in range(len(desc_stage))],
             positions=positions - offset, widths=0.3,
             patch_artist=True, boxprops=box_props_U,
             medianprops=median_props, showfliers=False)

plt.boxplot([gene_phen_matrix_S[:, i] for i in range(len(desc_stage))],
             positions=positions, widths=0.3,
             patch_artist=True, boxprops=box_props_S,
             medianprops=median_props, showfliers=False)

plt.boxplot([gene_phen_matrix_hS[:, i] for i in range(len(desc_stage))],
             positions=positions + offset, widths=0.3,
             patch_artist=True, boxprops=box_props_hS,
             medianprops=median_props, showfliers=False)

# Configuración general
plt.xticks(np.arange(0, len(desc_stage)), desc_stage, rotation=90, fontsize=13)
plt.xlabel("Stage", fontsize=13)
plt.ylabel("% of associated phenotypes", fontsize=13)

# Leyenda con colores planos
legend_patches = [
    mpatches.Patch(color="royalblue", label="U"),
    mpatches.Patch(color="mediumturquoise", label="S"),
    mpatches.Patch(color="tomato", label="hS")
]
plt.legend(handles=legend_patches, loc="upper left")

plt.tight_layout()
plt.savefig(path+'gene_stage_phenotype.png', dpi=600, bbox_inches='tight')
plt.show()




#5.) Stage of max expresion and lethality
#5.1.) Stage of maximal expression
g, ind_bulk_hS, i=np.intersect1d(genes_bulk, common_hS_genes, return_indices=True)
g, ind_bulk_S, i=np.intersect1d(genes_bulk, common_S_genes, return_indices=True)
g, ind_bulk_U, i=np.intersect1d(genes_bulk, common_U_genes, return_indices=True)

bulk_U=mean_bulk_matrix[ind_bulk_U, :]
bulk_S=mean_bulk_matrix[ind_bulk_S, :]
bulk_hS=mean_bulk_matrix[ind_bulk_hS, :]

#5.2.) For each gene we find the stage of maximal exp 
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



#5.3.) Lethal genes
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


#5.3.1.) We find the ENS genes associated with 'lethal
lethal_ind=np.where(phen_thing=='lethal (sensu genetics)')[0]
print('Associations with lethal:', len(lethal_ind))
phen_id_associated_with_lethal=phen_associated_desc[lethal_ind]
phen_id_associated_with_lethal_unique, times_letal_phen=np.unique(phen_id_associated_with_lethal, return_counts=True)
lethal_zfin_genes=zfin_associated[lethal_ind]
common_g_lethal_zfin, ind_all_lethal, i=np.intersect1d(zfin_all, lethal_zfin_genes, return_indices=True)
lethal_genes=ENS_all[ind_all_lethal]

#5.4) Association of lethal genes and maximal expression
U_genes_lethal, i, ind_U_stage=np.intersect1d(lethal_genes, common_U_genes, return_indices=True)
S_genes_lethal, i, ind_S_stage=np.intersect1d(lethal_genes, common_S_genes, return_indices=True)
hS_genes_lethal, i, ind_hS_stage=np.intersect1d(lethal_genes, common_hS_genes, return_indices=True)

max_exp_U_lethal=stage_of_max_exp_U[ind_U_stage]
max_exp_S_lethal=stage_of_max_exp_S[ind_S_stage]
max_exp_hS_lethal=stage_of_max_exp_hS[ind_hS_stage]


#mixed fig
time_unique=np.array(['0hpf', '0.75hpf', '2.25hpf', '3hpf', '4.3hpf', '5.25hpf', '6hpf', '8hpf', '10.3hpf', '16hpf', '19hpf', '24hpf', '30hpf', '36hpf', '2dpf', '3dpf', '4dpf', '5dpf'])

act_stages_lethal_U=np.zeros(len(time_unique))
act_stages_lethal_S=np.zeros(len(time_unique))
act_stages_lethal_hS=np.zeros(len(time_unique))
for j in (max_exp_U_lethal):
    act_stages_lethal_U[int(j)]=act_stages_lethal_U[int(j)]+1
for j in (max_exp_S_lethal):
    act_stages_lethal_S[int(j)]=act_stages_lethal_S[int(j)]+1
for j in (max_exp_hS_lethal):
    act_stages_lethal_hS[int(j)]=act_stages_lethal_hS[int(j)]+1

act_stages_lethal_U=act_stages_lethal_U/len(U_genes_lethal)
act_stages_lethal_S=act_stages_lethal_S/len(S_genes_lethal)
act_stages_lethal_hS=act_stages_lethal_hS/len(hS_genes_lethal)

width = 0.25  # ancho de cada barra
x = np.arange(len(time_unique))
plt.figure(figsize=(11, 6), dpi=600)
plt.bar(x - width, act_stages_lethal_U * 100, width=width, label='U', color='royalblue')
plt.bar(x,         act_stages_lethal_S * 100, width=width, label='S', color='mediumturquoise')
plt.bar(x + width, act_stages_lethal_hS * 100, width=width, label='hS', color='tomato')

plt.xticks(ticks=x, labels=time_unique, rotation=90, fontsize=25)
plt.yticks(fontsize=25)
plt.xlabel('Stages', fontsize=30)
plt.ylabel('% lethal genes\nwith max exp', fontsize=28, fontweight='bold')
plt.legend(fontsize=20)

plt.tight_layout()
plt.savefig(path + 'frac_lethal_genes_stage_maximal_exp.png', dpi=600, bbox_inches='tight')
plt.show()


#5.5.) Percentage of genes (from those with max expression) taht are lethal
#lethal genes
U_genes_lethal, i, ind_U_stage=np.intersect1d(lethal_genes, common_U_genes, return_indices=True)
S_genes_lethal, i, ind_S_stage=np.intersect1d(lethal_genes, common_S_genes, return_indices=True)
hS_genes_lethal, i, ind_hS_stage=np.intersect1d(lethal_genes, common_hS_genes, return_indices=True)

max_exp_U_lethal=stage_of_max_exp_U[ind_U_stage]
max_exp_S_lethal=stage_of_max_exp_S[ind_S_stage]
max_exp_hS_lethal=stage_of_max_exp_hS[ind_hS_stage]


#mixed fig
time_unique=np.array(['0hpf', '0.75hpf', '2.25hpf', '3hpf', '4.3hpf', '5.25hpf', '6hpf', '8hpf', '10.3hpf', '16hpf', '19hpf', '24hpf', '30hpf', '36hpf', '2dpf', '3dpf', '4dpf', '5dpf'])

#all the developmental genes
act_stages_U=np.zeros(len(time_unique))
act_stages_S=np.zeros(len(time_unique))
act_stages_hS=np.zeros(len(time_unique))
for j in (stage_of_max_exp_U):
    act_stages_U[int(j)]=act_stages_U[int(j)]+1
for j in (stage_of_max_exp_S):
    act_stages_S[int(j)]=act_stages_S[int(j)]+1
for j in (stage_of_max_exp_hS):
    act_stages_hS[int(j)]=act_stages_hS[int(j)]+1

#lethal genes
act_stages_lethal_U=np.zeros(len(time_unique))
act_stages_lethal_S=np.zeros(len(time_unique))
act_stages_lethal_hS=np.zeros(len(time_unique))
for j in (max_exp_U_lethal):
    act_stages_lethal_U[int(j)]=act_stages_lethal_U[int(j)]+1
for j in (max_exp_S_lethal):
    act_stages_lethal_S[int(j)]=act_stages_lethal_S[int(j)]+1
for j in (max_exp_hS_lethal):
    act_stages_lethal_hS[int(j)]=act_stages_lethal_hS[int(j)]+1

act_stages_lethal_U=act_stages_lethal_U/act_stages_U
act_stages_lethal_S=act_stages_lethal_S/act_stages_S
act_stages_lethal_hS=act_stages_lethal_hS/act_stages_hS

width = 0.25  # ancho de cada barra
x = np.arange(len(time_unique))
plt.figure(figsize=(11, 6), dpi=600)
plt.bar(x - width, act_stages_lethal_U * 100, width=width, label='U', color='royalblue')
plt.bar(x,         act_stages_lethal_S * 100, width=width, label='S', color='mediumturquoise')
plt.bar(x + width, act_stages_lethal_hS * 100, width=width, label='hS', color='tomato')

plt.xticks(ticks=x, labels=time_unique, rotation=90, fontsize=25)
plt.yticks(fontsize=25)
plt.xlabel('Stages', fontsize=30)
plt.ylabel('from genes with\nmax exp, % lethal', fontsize=28, fontweight='bold')
plt.legend(fontsize=20)

plt.tight_layout()
plt.savefig(path + 'percentage_of_those_with_max_exp_are_lethal.png', dpi=600, bbox_inches='tight')
plt.show()





