
import obonet
import numpy as np
import os
import networkx as nx
import matplotlib.pyplot as plt 
import pandas as pd
import matplotlib as mpl
from scipy.cluster.hierarchy import dendrogram
# import gran_libreria_phenotypes as glp
from matplotlib import cm 
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.manifold import MDS
from scipy.stats import pearsonr
from collections import defaultdict
from scipy.stats import mstats, kstest, ttest_ind, fisher_exact
from scipy.stats import linregress
from collections import defaultdict
import matplotlib.pyplot as plt
from collections import Counter

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


path='/media/alicia/TOSHIBA EXT/zebrafish/ontology/'
path_sc='/media/alicia/TOSHIBA EXT/zebrafish/zebrafish_code/single_cell/'
path_bulk='/media/alicia/TOSHIBA EXT/zebrafish/zebrafish_code/bulk/'


#0.) Charge data
gene_phen_association_matrix = np.loadtxt(path + 'gene_phen_association_matrix.txt', dtype=float)
unique_genes_ini = np.loadtxt(path + 'genes_associated_phen.txt', dtype=str, delimiter=',')
zfa_unique_desc_matrix = np.loadtxt(path + 'phen_anatomy.txt', dtype=str, delimiter=',')


pleio = np.sum(gene_phen_association_matrix, axis=1)

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
            # n_steps=final-int(ind_final)
            # if n_steps>1:
            #     for count in range(int(ind_final), n_steps):
            #         total_zfin_per_stage[count]=total_zfin_per_stage[count]+1
            #         print(count, int(ind_ini), n_steps)
            # else:
            # print('hola', print(ind_ini))
        if len(ind_ini)>0:
            n_zfin=n_zfin+1
            n_steps=final-int(ind_ini)
            for count in range(int(ind_ini), n_steps):
                total_zfin_per_stage[count]=total_zfin_per_stage[count]+1
                # print(count, int(ind_ini), n_steps)
        
    

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



# 2.) gene-phen associations. also keep publications
with open(path+'phen_gene_data.txt', 'r') as file:
    data = file.read()
rows = data.split('\n')
rows = [row.strip() for row in rows if row.strip()]
del data, file


temp_pub_ids = []
temp_zfin_associated = []
temp_phen_associated_id = []
temp_phen_associated_desc = []
temp_phen_thing = []


zfin_to_ens_dict = defaultdict(list)
for z_id, e_id in zip(zfin_all, ENS_all):
    zfin_to_ens_dict[z_id].append(e_id)


print("Leyendo archivo de fenotipos...")
for elemento in rows:
    partes = elemento.split("\t")
    pub_id = next((p for p in partes if p.startswith("ZDB-PUB-")), None)
    
    for i in range(len(partes)):
        if partes[i].startswith("ZFA"):
            gene_zfin_id = partes[2]
            
            temp_pub_ids.append(pub_id)
            temp_zfin_associated.append(gene_zfin_id)
            temp_phen_associated_id.append(partes[i])
            temp_phen_associated_desc.append(partes[int(i+1)])
            temp_phen_thing.append(partes[10])
            break # Encontramos el ZFA, pasamos a la siguiente línea

del rows, partes

#plot
#=====================================================
ens_to_pubs_unfiltered = defaultdict(set)

for pub_id, zfin_id in zip(temp_pub_ids, temp_zfin_associated):
    if pub_id:
        ens_ids = zfin_to_ens_dict.get(zfin_id, [])
        for ens_id in ens_ids:
            ens_to_pubs_unfiltered[ens_id].add(pub_id)

pubs_per_gene_ini = np.array([len(ens_to_pubs_unfiltered[gen]) for gen in unique_genes_ini])

x = pleio
y = pubs_per_gene_ini

slope, intercept, r_value, p_value, std_err = linregress(x, y)
y_pred = slope * x + intercept
residuos = abs(y - y_pred)

umbral_99 = np.percentile(residuos, 99)

outliers = residuos >= umbral_99
datos_normales = ~(outliers)

plt.figure(figsize=(3.5, 2.5), dpi=600)
plt.scatter(x[datos_normales], y[datos_normales], alpha=0.5, color='dodgerblue', s=10)
plt.plot(x, y_pred, color='purple', linestyle='--', linewidth=1, label='Regression line')
plt.scatter(x[outliers], y[outliers], color='slateblue', s=10)
plt.xlabel('Anatomical pleiotropy', fontsize=12, fontweight='bold')
plt.ylabel('# publications', fontsize=12, fontweight='bold')
plt.legend(loc='upper left') # Movido a la izquierda para que no tape los hiper-estudiados
plt.grid(False)
plt.tight_layout()
plt.show()

#=====================================================

#3.Genes per pub
pub_to_ens_genes = defaultdict(set)

for pub_id, zfin_id in zip(temp_pub_ids, temp_zfin_associated):
    if pub_id:
        # Buscamos el ENS asociado a este ZFIN
        ens_ids = zfin_to_ens_dict.get(zfin_id, [])
        for ens_id in ens_ids:
            pub_to_ens_genes[pub_id].add(ens_id)

#We filter the publications by the number of genes associated per pub
pubs_filtradas = {pub for pub, genes_ens in pub_to_ens_genes.items() if (len(genes_ens) > 100) & (len(genes_ens) < 400)}

n_genes_per_pub = [len(genes) for genes in pub_to_ens_genes.values()]
plt.figure(figsize=(3, 2), dpi=600)
plt.hist(n_genes_per_pub, bins=100, log=True, color='red')
# plt.title("Distribución de Genes ENS únicos por Publicación")
plt.xlabel("# genes", fontsize=13, fontweight='bold')
plt.ylabel("# publications", fontsize=13, fontweight='bold')
plt.show()

np.percentile(n_genes_per_pub, 99.5)


ENS_associated = []
phen_associated_id_final = []
phen_associated_desc_final = []
phen_thing_final = []

count_rep_genes = 0

for i in range(len(temp_zfin_associated)):
    pub_id = temp_pub_ids[i]
    
    if pub_id not in pubs_filtradas:
        continue
        
    zfin_id = temp_zfin_associated[i]
    ens_ids = zfin_to_ens_dict.get(zfin_id, [])
    
    if len(ens_ids) == 1:
        phen_associated_id_final.append(temp_phen_associated_id[i])
        phen_associated_desc_final.append(temp_phen_associated_desc[i])
        phen_thing_final.append(temp_phen_thing[i])
        ENS_associated.append(ens_ids[0])
        
    elif len(ens_ids) > 1:
        count_rep_genes += 1
        for ens_id in ens_ids:
            phen_associated_id_final.append(temp_phen_associated_id[i])
            phen_associated_desc_final.append(temp_phen_associated_desc[i])
            phen_thing_final.append(temp_phen_thing[i])
            ENS_associated.append(ens_id)

ENS_associated = np.array(ENS_associated)
phen_associated_id_final = np.array(phen_associated_id_final)
phen_associated_desc_final = np.array(phen_associated_desc_final)
phen_thing_final = np.array(phen_thing_final)

del temp_pub_ids, temp_zfin_associated, temp_phen_associated_id, temp_phen_associated_desc, temp_phen_thing



#5.) COMPARISON OF PLEIO of filtered pub vs our pleio 
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


# 5.1.) Build the graph and its reverse
graph = defaultdict(list)
reverse_graph = defaultdict(list)
for parent, child in zip(zfa_parent, zfa_children):
    graph[parent].append(child)
    reverse_graph[child].append(parent)

#5.2.) Find ancestors for a given term
results = {}
ancestors_list=[]
for term in zfa_unique:
    results[term] = find_ancestors(reverse_graph, term)
    ancestors = find_ancestors(reverse_graph, term)
    ancestors_list.append([term, list(ancestors)])

#5.2.1.) We find all the needed terms
all_redundant_terms=[]
for i in range(len(zfa_unique)):
    all_redundant_terms.append(ancestors_list[i][0])
    for j in range(len(ancestors_list[i][1])):
        all_redundant_terms.append(ancestors_list[i][1][j])

all_redundant_terms_unique=np.unique(np.array(all_redundant_terms))


#5.3.) We find the layers and the phenotypes
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
    

        

#5.3.1.) We sort the associated phenotypes by layer
unique_genes_new=np.unique(ENS_associated)

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


#5.5) We built the gene-phenotype association matrix
gene_phen_association_matrix1=np.zeros((len(unique_genes_new), len(zfa_unique_matrix)))
for i in range(len(unique_genes_new)):
    ind_association=np.where(ENS_associated==unique_genes_new[i])
    searched_phen=phen_associated_id_final[ind_association]
    for j in range(len(searched_phen)):
        ind_phen=np.where(zfa_unique==searched_phen[j])[0]
        ind_matrix=np.where(zfa_unique_matrix==zfa_unique[int(ind_phen)])[0]
        gene_phen_association_matrix1[i][int(ind_matrix)]=1
        for k in range(len(ancestors_list[int(ind_phen)][1])):
            ind_matrix=np.where(zfa_unique_matrix==ancestors_list[int(ind_phen)][1][k])[0]
            gene_phen_association_matrix1[i][int(ind_matrix)]=1

    print(i)


#5.6) 
pleiotropy_scores_new = np.sum(gene_phen_association_matrix1, axis=1)

c, ind, ind_new=np.intersect1d(unique_genes_ini, unique_genes_new, return_indices=True)

print(pearsonr(pleio[ind], pleiotropy_scores_new[ind_new]))


print(pubs_filtradas)

plt.figure(figsize=(4, 3), dpi=600)
plt.scatter(pleio[ind], pleiotropy_scores_new[ind_new], s=7, color='forestgreen')
plt.xlabel("Pleiotropy (complete dataset)", fontsize=13, fontweight='bold')
plt.ylabel("Pleiotropy (large-scale\nmutagenesis screen)", fontsize=13, fontweight='bold')
plt.show()


#6.) Types of associations
archivo_datos = path+'phen_gene_data.txt'

conteo_mutaciones = Counter()
total_lineas_validas = 0


with open(archivo_datos, 'r') as f:
    for linea in f:
        if not linea.strip():
            continue
            
        partes = linea.strip().split('\t')
        
        fish_idx = next((i for i, p in enumerate(partes) if p.startswith("ZDB-FISH-")), None)
        
        if fish_idx is not None and (fish_idx + 1) < len(partes):
            genotipo = partes[fish_idx + 1]
            total_lineas_validas += 1
            
            es_mo = 'MO' in genotipo
            es_tg = 'Tg' in genotipo
            es_mutante = '<sup>' in genotipo
            
            if es_mo:
                conteo_mutaciones['Morpho (MO)'] += 1
            elif es_mutante and es_tg:
                conteo_mutaciones['Stable Mutant + Transgene'] += 1
            elif es_tg:
                conteo_mutaciones['Only Trasngene (Tg)'] += 1
            elif es_mutante:
                conteo_mutaciones['Pure mutation (Homo/Hetero)'] += 1
            else:
                conteo_mutaciones['Others / WT / Non class'] += 1

print("\n--- RESULTS ---")
print(f"Total associations: {total_lineas_validas}")

for categoria, cantidad in conteo_mutaciones.most_common():
    porcentaje = (cantidad / total_lineas_validas) * 100
    print(f"{categoria}: {cantidad} ({porcentaje:.2f}%)")

